import tables

from matplotlib import patches
import matplotlib.mlab as ml
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pickle as pickle
import os
from scipy import interpolate
import matplotlib.pyplot as plt
from PIL import Image
import astropy.io.fits as pyfits
from scipy.interpolate import griddata
import pyzdde.arraytrace as at
from collections import Counter

px = []
py = []
for i in range(-50, 51, 1):
    for j in range(-50, 51, 1):
        px.append(i / 50.)
        py.append(j / 50.)
px = np.array(px)
py = np.array(py)

idx = (px ** 2 + py ** 2) < 1


class Transformation(tables.IsDescription):
    wavelength = tables.Float32Col()
    shear = tables.Float32Col()
    rotation = tables.Float32Col()
    scale_x = tables.Float32Col()
    scale_y = tables.Float32Col()
    translation_x = tables.Float32Col()
    translation_y = tables.Float32Col()


def save_CCD_info_to_hdf(path, ccd):
    h5file = tables.open_file(path, "a")
    ccd_group = h5file.create_group(h5file.root, 'CCD', 'CCD information')
    ccd_group._v_attrs.Nx = ccd.Nx
    ccd_group._v_attrs.Ny = ccd.Ny
    ccd_group._v_attrs.pixelsize = ccd.pixelSize
    h5file.close()


def save_spectrograph_info_to_hdf(path, spec):
    h5file = tables.open_file(path, "w")
    spec_group = h5file.create_group(h5file.root, 'Spectrograph', "Spectrograph Information")
    spec_group._v_attrs.blaze = spec.blaze
    spec_group._v_attrs.gpmm = spec.grmm
    spec_group._v_attrs.name = spec.name
    h5file.close()


def save_transformation_to_hdf(path, res, fiber_number=1):
    h5file = tables.open_file(path, "a")
    gr = h5file.create_group(h5file.root, "fiber_" + str(fiber_number))
    gr._v_attrs.MatricesPerOrder = res['MatricesPerOrder']
    gr._v_attrs.norm_field = res['norm_field']
    gr._v_attrs.sampling_input_x = res['sampling_input_x']
    gr._v_attrs.field_with = res['field_width']
    gr._v_attrs.field_height = res['field_height']

    for order, r in res['matrices'].iteritems():
        tab = h5file.create_table("/fiber_" + str(fiber_number), 'order' + str(abs(order)), Transformation,
                                  "Affine Transformation", expectedrows=len(r), chunkshape=True)
        transf = tab.row
        for wl, pars in r.iteritems():
            transf['wavelength'] = wl
            transf['rotation'] = pars[0]
            transf['scale_x'] = pars[1]
            transf['scale_y'] = pars[2]
            transf['shear'] = pars[3]
            transf['translation_x'] = pars[4]
            transf['translation_y'] = pars[5]

            transf.append()
        tab.flush()
    h5file.close()


def save_psfs_to_hdf(path, res, fiber_number=1):
    h5file = tables.open_file(path, "a")
    if not (h5file.__contains__("/fiber_" + str(fiber_number))):
        gr = h5file.create_group(h5file.root, "fiber_" + str(fiber_number))
    else:
        gr = h5file.get_node(h5file.root, "fiber_" + str(fiber_number))

    for order, psfs in res.iteritems():
        if not (h5file.__contains__("/fiber_" + str(fiber_number) + "/psf_order_" + str(abs(order)))):
            gr = h5file.create_group("/fiber_" + str(fiber_number), "psf_order_" + str(abs(order)))
        else:
            gr = h5file.get_node("/fiber_" + str(fiber_number), "psf_order_" + str(abs(order)))

        for wl, data in psfs.iteritems():
            if not (
                    h5file.__contains__(
                        "/fiber_" + str(fiber_number) + "/psf_order_" + str(order) + "/wavelength_" + str(wl))):
                ar = h5file.create_array(gr, "wavelength_" + str(wl), np.array(data[1]))
                ar.attrs.wavelength = float(wl)
                ar.attrs.order = int(abs(order))
                for i, a in enumerate(data[0]._fields):
                    ar.set_attr(a, data[0][i])


def efficiency(scalingfactor, order, alpha, blaze, wl, n):
    bb = np.arcsin(-np.sin(alpha) + order * wl * n * 1E-6)
    return scalingfactor * np.sinc(order * (np.cos(alpha) / np.cos(alpha - blaze)) *
                                   (np.cos(blaze) - np.sin(blaze) / np.tan((alpha + bb) / 2.))) ** 2


class Spot(object):
    """ Class that describes a spot in a optical design

        It basically consists of a DDEArray
     """

    def __init__(self, wavelength, order, i, rd_in, rd_out, valid_only=True, circular_pupil=True):
        """
        Constructor
        :param wavelength: wavelength in microns
        :param order: order of diffraction of the echelle grating
        :param i: index of spot per order - makes it easier to create the spot_map but is probably redundant
        :param rd_in: DdeArray of input rays before raytracing
        :param rd_out: DdeArray of traced rays
        :param valid_only: if True, only rays within a circular aperture are traced (needed for spot diagrams)
                which are not vignetted
        :return:
        """

        self.wavelength = wavelength
        self.order = order
        self.i = i
        # number of rays
        self.Nrays = len(rd_in['z'][1:])

        # restrict rays to circular pupil or not
        if circular_pupil:
            px = rd_in['z'][1:]
            py = rd_in['l'][1:]
            idx = (px ** 2 + py ** 2) <= 1.
        else:
            idx = np.ones(self.Nrays)

        # restrict rays to non vignetted ones
        if valid_only:
            vig = rd_out['vigcode'][1:]
            err = rd_out['error'][1:]
            vig = np.logical_or(vig, err)
            index = np.logical_and(vig < 1, idx)
        else:
            index = idx

        self.hx = rd_in['x'][1:][index]
        self.hy = rd_in['y'][1:][index]
        self.x = rd_out['x'][1:][index]
        self.y = rd_out['y'][1:][index]

        self.px = rd_in['z'][1:][index]
        self.py = rd_in['l'][1:][index]

        self.barycenter = None
        self.xy_c = None
        self.rms = None
        self.rms_x = None
        self.rms_y = None

        self._calc_barycenter()
        self._calc_rms_radius()

    def _calc_barycenter(self):
        """
        calculate the barycenter of the spot
        """
        self.barycenter = {'x': np.average(self.x),
                           'y': np.average(self.y)}
        self.xy_c = {'x': self.x - self.barycenter['x'],
                     'y': self.y - self.barycenter['y']}

    def _calc_rms_radius(self):
        """
        calculate rms radius of the spot, radially, in x and y direction
        """
        self.rms = np.std(np.sqrt(self.xy_c['x'] ** 2 + self.xy_c['y'] ** 2))
        self.rms_x = np.std(np.sqrt(self.xy_c['x'] ** 2))
        self.rms_y = np.std(np.sqrt(self.xy_c['y'] ** 2))

    def EE_radius(self, EE=80., direction='r'):
        """
        Calculate encircled energy (EE) radius of the spot

        :param EE: encircled energy level in percent
        :param direction: direction in which EE is calculated (radial, x and y)
        :return:
        """
        n = len(self.xy_c['x'])
        if direction == 'r':
            return np.sort(np.sqrt(self.xy_c['x'] ** 2 + self.xy_c['y'] ** 2))[int(EE / 100. * n)] * 1000.
        if direction == 'x':
            return np.sort(np.sqrt(self.xy_c['x'] ** 2))[int(EE / 100. * n)] * 1000.
        if direction == 'y':
            return np.sort(np.sqrt(self.xy_c['y'] ** 2))[int(EE / 100. * n)] * 1000.

    def calc_weighted_barycenter(self, path_image=None, xy_c=None, radius=None, f=None, plot=False):
        """
        Calculates the barycenter of the spot weighted with an image.
        This can be used to calculate the spot barycenter weighted with a fiber far field (FF) illumination pattern.

        :param path_image: path to image that contains the weights
        :param xy_c: x and y coordinate of the center of the FF for interpolation, default is geometric image center
        :param radius: radius on the FF image that corresponds to p=radius, default is half image width
        :return: weighted barycenter
        """

        if isinstance(path_image, str):
            if path_image.lower().endswith('.fit') or path_image.lower().endswith('.fits'):
                weight_image = pyfits.open(path_image)[0].data[xy_c['y'] - np.ceil(radius):xy_c['y'] + np.ceil(radius),
                               xy_c['x'] - np.ceil(radius):xy_c['x'] + np.ceil(radius)]
            else:

                if xy_c == None:
                    xy_c = {}
                    dims = np.shape(np.array(Image.open(path_image).convert('L')))
                    xy_c['y'] = dims[0] / 2.
                    xy_c['x'] = dims[1] / 2.

                if radius == None:
                    radius = np.shape(np.array(Image.open(path_image).convert('L')))[0] / 2.

                # open image but only select relevant parts
                weight_image = np.array(Image.open(path_image).convert('L'))[
                               xy_c['y'] - np.ceil(radius):xy_c['y'] + np.ceil(radius),
                               xy_c['x'] - np.ceil(radius):xy_c['x'] + np.ceil(radius)]

            # normalized x and y coordinates (correspond to Px and Py in ZEMAX)
            xr = yr = np.arange(-1., 1., 1. / radius)

            # interpolation function
            f = interpolate.RectBivariateSpline(xr, yr, weight_image)

        w = f(self.px, self.py, grid=False)

        weighted_barycenter = {'x': np.average(self.x, weights=w),
                               'y': np.average(self.y, weights=w)}

        if plot:
            plt.figure()
            plt.scatter(self.px, self.py, c=w, linewidth=0., marker='o')
            plt.show()

        return weighted_barycenter


class Order(object):
    """ Class that describes an echelle order
    """

    def __init__(self, m, blazeWL, minWL, maxWL, minFSRwl, maxFSRwl):
        """
        Constructor

        :param m: order number
        :param blazeWL: blaze wavelength [micron]
        :param minWL: minimum wavelength that fits on chip [micron]
        :param maxWL: maximum wavelength that fits on chip [micron]
        :param minFSRwl: minimum FSR wavelength [micron]
        :param maxFSRwl: maximum FSR wavelength [micron]
        :return: None
        """
        self.m = m
        self.blazeWL = blazeWL
        self.minWL = minWL
        self.maxWL = maxWL
        self.minFSRwl = minFSRwl
        self.maxFSRwl = maxFSRwl

    def inFSR(self, wl):
        """
        checks if wavelength lies within FSR or not

        :param wl: wavelength [micron]
        :return: True/False
        """
        return self.maxFSRwl > wl > self.minFSRwl

    def inOrder(self, wl):
        """
        checks if wavelength lies in order (all chip) or not

        :param wl: wavelength [micron]
        :return: True/False
        """
        return self.maxWL > wl > self.minWL

    def info(self):
        print('Order ', self.m)
        print('FSR wavelength boundaries [microns]', self.minFSRwl, self.maxFSRwl)
        print('Chip wavelength boundaries [microns]', self.minWL, self.maxWL)


class CCD(object):
    """ CCD class, contains information about CCD detector """

    def __init__(self, Nx, Ny, pixelSize, dispersionDirection='x', name=''):
        """

        :param Nx: number of pixels in x - direction
        :param Ny: number of pixels in y - direction
        :param pixelSize: size of one pixel [micron]
        :param dispersionDirection: echelle dispersion direction
        :param name: name/identifier of the CCD detector
        :return:
        """
        self.Wx = Nx * pixelSize / 1000.
        self.Wy = Ny * pixelSize / 1000.
        self.Nx = Nx
        self.Ny = Ny
        self.pixelSize = pixelSize
        self.name = name

        self.xi = np.linspace(-Nx * pixelSize / 2000., Nx * pixelSize / 2000., Nx)
        self.yi = np.linspace(-Ny * pixelSize / 2000., Ny * pixelSize / 2000., Ny)

        self.extent = [-Nx * pixelSize / 2000.,
                       +Nx * pixelSize / 2000.,
                       -Ny * pixelSize / 2000.,
                       +Ny * pixelSize / 2000.]

        self.shape = [[-Nx * pixelSize / 2000., -Ny * pixelSize / 2000.],
                      [Nx * pixelSize / 2000., -Ny * pixelSize / 2000.],
                      [Nx * pixelSize / 2000., Ny * pixelSize / 2000.],
                      [-Nx * pixelSize / 2000., Ny * pixelSize / 2000.]
                      ]
        self.dispersionDirection = dispersionDirection


class Echelle():
    """
    class describing an echelle spectrograph
    """

    def __init__(self, ln=None, name=''):
        self.name = name
        self.savePath = 'PyEchelle_' + self.name
        if not os.path.exists(self.savePath):
            os.makedirs(self.savePath)
        # zemax surface number
        # self.ln= pyz.createLink()
        if ln is not None:
            import pyzdde.zdde as pyz
            import pyzdde.arraytrace as at  # Module for array ray tracing
            self.ln = ln
        self.zmx_nsurf = None
        # minimal/maximal order
        self.minord = None
        self.maxord = None
        # Blaze angle in degree
        self.blaze = None
        # gamma angle in degree
        self.gamma = None
        # groves per mm
        self.grmm = None
        # current order
        self.order = None
        self.theta = 0
        self.grp = None

        self.tracing = []

        self.x = []
        self.y = []
        self.orders = []

        self.file = None

        self.rays = []
        self.wls = []

        self.CCD = None
        self.Orders = {}
        self.spots = []

        self.order_configs = {}
        self.order_config_wave = {}

    def setCCD(self, CCD):
        self.CCD = CCD

    def saveOrders(self, filename='orders.pkl'):
        """
        Save Orders to file
        :param filename: filename
        :return: None
        """
        print('save orders')
        pickle.dump(self.Orders, open(self.savePath + '/' + filename, "wb"))

    def saveSpectrograph(self, filename=None):
        if filename == None:
            filename = self.name
        spec = {'blaze': self.blaze,
                'gamma': self.gamma,
                'theta': self.theta,
                'order': self.order,
                'name': self.name,
                'savePath': self.savePath,
                'minOrder': self.minord,
                'maxOrder': self.maxord,
                'grmm': self.grmm,
                'grp': self.grp,
                }

        pickle.dump(spec, open(self.savePath + '/' + filename + '.pkl', "wb"))

    def loadSpectrograph(self, filename=None):
        if filename == None:
            filename = self.name
        spec = pickle.load(open(self.savePath + '/' + filename + '.pkl'))
        self.blaze = spec['blaze']
        self.gamma = spec['gamma']
        self.theta = spec['theta']
        self.order = spec['order']
        self.minord = spec['minOrder']
        self.maxord = spec['maxOrder']
        self.grmm = spec['grmm']
        self.grp = spec['grp']

    def loadOrders(self, filename='orders.pkl'):
        """
        Load Orders from file
        :param filename: filename
        :return:
        """
        self.Orders = pickle.load(open(self.savePath + '/' + filename))

    def analyseZemaxFile(self, echellename='Echelle', thetaname='theta', blazename='blaze', gammaname='gamma'):
        """
        Analyses ZEMAX files and extract important parameters to specify Echelle Spectrograph.
        Looks for names in comment column of ZEMAX to detect specific surfaces.

        :param echellename: ZEMAX surface name of Echelle grating
        :param thetaname: ZEMAX surface name of theta angle
        :param blazename: ZEMAX surface name of blaze angle
        :param gammaname: ZEMAX surface name of gamma angle
        :return:
        """
        for i in range(self.ln.zGetNumSurf()):
            comm = self.ln.zGetComment(i)
            if comm == echellename:
                print('Echelle found ----------------------------')
                self.zmx_nsurf = i
                self.echelle_surface = i
                # grooves per mm
                self.grmm = self.ln.zGetSurfaceParameter(i, 1) * 1000.
                # current order
                self.order = int(self.ln.zGetSurfaceParameter(i, 2))
                print('Grooves per mm', self.grmm)
                print('Current order', self.order)
                print('Surface number', self.zmx_nsurf)
            elif comm == thetaname:
                print('Theta found ------------------------------')
                self.theta = float(self.ln.zGetSurfaceParameter(i, 4))
                print('theta', self.theta)

            elif comm == blazename:
                print('blaze found ------------------------------')
                b1 = abs(float(self.ln.zGetSurfaceParameter(i, 3)))
                b2 = abs(float(self.ln.zGetSurfaceParameter(i, 4)))
                b3 = abs(float(self.ln.zGetSurfaceParameter(i, 5)))
                self.blaze = max((b1, b2, b3))
                print('blaze', self.blaze)

            elif comm == gammaname:
                print('gamma found ------------------------------')
                b1 = abs(float(self.ln.zGetSurfaceParameter(i, 3)))
                b2 = abs(float(self.ln.zGetSurfaceParameter(i, 4)))
                self.gamma = max((b1, b2))
                print('gamma', self.gamma)

    def trace(self, wave=1, hx=0, hy=0, N=101, intensity=1.):
        self.ln.zGetUpdate()
        self.ln.zPushLens()
        Nx = Ny = int(np.sqrt(N))
        rd_in = at.getRayDataArray(Nx * Ny, tType=0, mode=0)
        rd_out = at.getRayDataArray(Nx * Ny, tType=0, mode=0)
        k = 0
        for i in np.linspace(-1., 1., Nx):
            for j in np.linspace(-1., 1., Ny):
                k += 1
                rd_out[k].x = hx
                rd_out[k].y = hy
                rd_out[k].z = i  # px
                rd_out[k].l = j  # py
                rd_out[k].intensity = intensity
                rd_out[k].wave = wave

                rd_in[k].x = hx
                rd_in[k].y = hy
                rd_in[k].z = i  # px
                rd_in[k].l = j  # py
                rd_in[k].intensity = intensity
                rd_in[k].wave = wave

        ret = at.zArrayTrace(rd_out, timeout=5000)

        return np.array(rd_in, dtype=at.DdeArrayData), np.array(rd_out, dtype=at.DdeArrayData)

    def trace_rays(self, wave, FIELD):
        self.ln.zGetUpdate()
        self.ln.zPushLens()
        numRays = 10201
        rd = at.getRayDataArray(numRays, tType=0, mode=0)
        # Fill the rest of the ray data array
        k = 0
        for i in range(-50, 51, 1):
            for j in range(-50, 51, 1):
                k += 1
                rd[k].y = FIELD
                rd[k].z = i / 50.  # px
                rd[k].l = j / 50.  # py
                rd[k].intensity = 1.0
                rd[k].wave = wave

        # Trace the rays
        ret = at.zArrayTrace(rd, timeout=5000)
        return rd

    def order_to_config(self, order):
        return self.order_configs[order]

    def clear_configs(self):
        c, cc, rc = self.ln.zGetConfig()
        for i in range(cc):
            self.ln.zDeleteConfig(1)
        self.ln.zPushLens()
        for i in range(rc):
            self.ln.zDeleteMCO(1)

    def clear_merit_function(self):
        while (self.ln.zDeleteMFO(1) > 1):
            self.ln.zDeleteMFO(1)

    def set_config_and_wavelength(self, wavelength_per_order=7):
        self.clear_configs()

        self.ln.zSetMulticon(0, 1, 'PAR2', self.echelle_surface, 0, 0)
        self.ln.zInsertMCO(2)
        self.ln.zSetMulticon(0, 2, 'WAVE', 0, 0, 0)

        self.order_configs = {}
        for i, o in enumerate(self.Orders):
            self.ln.zInsertConfig(i + 1)
            self.ln.zSetMulticon(i + 1, 1, self.Orders[o].m, 0, 0, 0, 1, 0)
            self.ln.zSetMulticon(i + 1, 2, self.Orders[o].blazeWL, 0, 0, 0, 1, 0)
            self.order_configs[o] = i + 1
            # self.ln.zInsertMFO(i + 1)
            # self.ln.zSetOperandRow(i + 1, 'CONF', i+1)

        c, cc, rc = self.ln.zGetConfig()
        self.ln.zDeleteConfig(cc)
        self.ln.zPushLens()

    def clear_wavelength(self):
        n = self.ln.zGetNumWave()

    def set_config_and_wavelength_from_list(self, orders, wavelength, posx, posy):
        self.clear_configs()
        self.clear_merit_function()

        self.ln.zSetMulticon(0, 1, 'PAR2', self.echelle_surface, 0, 0)

        # add unique orders to multi config file
        unique_orders = np.unique(np.array(orders))
        self.order_configs = dict(zip(unique_orders, range(len(unique_orders))))
        for i, o in enumerate(unique_orders):
            self.ln.zInsertConfig(i + 1)
            self.ln.zSetMulticon(i + 1, 1, o, 0, 0, 0, 1, 0)
            self.order_configs[o] = i + 1

        self.ln.zPushLens()
        c, cc, rc = self.ln.zGetConfig()
        self.ln.zDeleteConfig(cc)
        self.ln.zPushLens()

        # # add as many rows needed for the order with the most wavelength
        n_wavelength = Counter(orders).most_common(1)[0][1]
        self.ln.zSetWave(0, 1, n_wavelength)
        self.ln.zGetUpdate()
        #

        for n in range(n_wavelength):
            self.ln.zInsertMCO(n + 2)
            self.ln.zSetPrimaryWave(n + 1)
            self.ln.zSetMulticon(0, n + 2, 'WAVE', n + 1, n + 1, n + 1)
            for i in unique_orders:
                self.ln.zSetMulticon(self.order_to_config(i), n + 2, self.Orders[i].blazeWL, 0, 0, 0, 1, 0)
        row_count = {}
        for uo in unique_orders:
            row_count[uo] = 2

        for o, wl, px, py in zip(orders, wavelength, posx, posy):
            config = self.order_to_config(o)
            rc = row_count[o]
            self.ln.zSetMulticon(config, rc, wl, 0, 0, 0, 1, 0)
            self.set_merit_function(o, rc - 1, px, py)
            row_count[o] += 1

        self.ln.zPushLens()

    def set_merit_function(self, order, wave, posx, posy, clear=False):
        if clear:
            self.clear_merit_function()
        n = self.ln.zGetNumSurf()
        cf = self.order_to_config(order)
        self.ln.zInsertMFO(1)
        self.ln.zSetOperandRow(1, 'REAY', n, wave, 0, 0, 0, 0, tgt=posy)

        self.ln.zInsertMFO(1)
        self.ln.zSetOperandRow(1, 'REAX', n, wave, 0, 0, 0, 0, tgt=posx)

        self.ln.zInsertMFO(1)
        self.ln.zSetOperandRow(1, 'CONF', cf)
        self.ln.zPushLens()

    def read_merit_function_position_difference(self, n):
        REAX = []
        REAY = []
        dx = []
        dy = []
        for i in range(n):
            data = self.ln.zGetOperandRow(i)
            if data[0] == 'REAX':
                dx.append((data[9] - data[11]) * data[10])
                REAX.append(data[11])
            if data[0] == 'REAY':
                dy.append((data[9] - data[11]) * data[10])
                REAY.append(data[11])
        print("Median deviation XY: ", np.median(np.array(dx)) * 1000., np.median(np.array(dy)) * 1000.)
        plt.figure()
        plt.plot()
        plt.axis('equal')

        for x, y, dxx, dyy in zip(REAX, REAY, dx, dy):
            plt.scatter(x, y)
            plt.arrow(x, y, dxx * 100, dyy * 100)
        plt.show()

    def do_spectral_format(self, nPerOrder=7, FSRonly=True, hx=0, hy=0):
        s = []
        for o in list(self.Orders.values()):
            print('Trace order', o.m)
            self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, o.m)
            # self.ln.zPushLens()
            if FSRonly:
                wl = np.linspace(o.minFSRwl, o.maxFSRwl, nPerOrder)
            else:
                wl = np.linspace(o.minWL, o.maxWL, nPerOrder)
            for i, w in enumerate(wl):
                self.ln.zSetWave(1, w, 1.)
                self.ln.zGetUpdate()
                # self.ln.zPushLens()
                rayTraceData = self.ln.zGetTrace(1, 0, -1, hx, hy, 0, 0)
                error, vig, x, y, z, l, m, n, l2, m2, n2, intensity = rayTraceData
                s.append([o.m, w, x, y])
        return s

    def get_psfs(self, nPerOrder=1, fieldnumber=3, fieldposition=[0., 0.]):
        psfs = {}
        old_field = self.ln.zGetField(fieldnumber)
        self.ln.zSetField(fieldnumber, fieldposition[0], fieldposition[1])
        for o in list(self.Orders.values()):
            print('Trace order', o.m)
            self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, o.m)
            wl = np.linspace(o.minWL, o.maxWL, nPerOrder)
            psfs[o.m] = {}
            for w in wl:
                self.ln.zSetWave(1, w, 1.)
                psf = self.ln.zGetPSF(which='huygens')
                print(psf)
                psfs[o.m][w] = psf
        # restore field
        self.ln.zSetField(fieldnumber, old_field.xf, old_field.yf)
        return psfs

    def do_affine_transformation_calculation(self, nPerOrder=10,
                                             norm_field=[[-1, 1], [-1, -1], [1, -1], [1, 1], [0, 0]], fw=None, fh=None):
        """
        Calculates Affine Matrices that describe spectrograph

        The spectrograph can be described by affine transformations from the input slit to the focal plane.
        an affine transofmration can be described by a 3x3 matrix.
        this function calculates the 3x3 matrix per wavelength and order that matches the input slit to the focal plane

        :param nPerOrder: number of wavelength steps across one order
        :param norm_field: corner points and center point in normalized coordinates
        :param fw: fiber/slit width [microns]
        :param fh: fiber/slit height [microns]
        :return:
        """
        from skimage import transform as tf
        sampling_input_x = int(fw)
        res = {'MatricesPerOrder': nPerOrder,
               'norm_field': norm_field,
               'sampling_input_x': sampling_input_x}

        # find field dimensions in ZEMAX
        field_info = self.ln.zGetField(0)

        # TODO: raise error
        if field_info[0] is not 1:
            print('Field coordinates have the wrong format')

        zmx_fields = []
        for ii in range(1, field_info[1] + 1):
            field = self.ln.zGetField(ii)
            zmx_fields.append([field[0], field[1]])
        zmx_fields = np.array(zmx_fields)
        norm_field = np.array(norm_field)

        if fw is None:
            fw = (np.max(zmx_fields[:, 0]) - np.min(zmx_fields[:, 0])) * 1000.  # slit width in microns
            fw *= (np.max(norm_field[:, 0]) - np.min(norm_field[:, 0])) / 2.
        if fh is None:
            fh = (np.max(zmx_fields[:, 1]) - np.min(zmx_fields[:, 1])) * 1000.  # slit height in microns
            fh *= (np.max(norm_field[:, 1]) - np.min(norm_field[:, 1]))
        print('Field width: ' + str(fw))
        print('Field height: ' + str(fh))

        res['field_width'] = fw
        res['field_height'] = fh
        sampling_x = sampling_input_x
        sampling_y = sampling_input_x * fh / fw

        src = np.array(norm_field, dtype=float)
        src[:, 0] -= np.min(src[:, 0])
        src[:, 1] -= np.min(src[:, 1])

        src[:, 0] /= np.max(src[:, 0])
        src[:, 1] /= np.max(src[:, 1])

        # src[:, 0] *= sampling_x
        # src[:, 1] *= sampling_y

        ppp = []
        dst_x = []
        dst_y = []
        orders = []
        wavelength = []
        for o in list(self.Orders.values()):
            print('trace order ' + str(o.m))
            wl = np.linspace(o.minWL, o.maxWL, nPerOrder)
            self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, o.m)
            # print(wl, o.m)
            for w in wl:
                self.ln.zSetWave(1, w, 1.)
                self.ln.zGetUpdate()
                for f in norm_field:
                    rayTraceData = self.ln.zGetTrace(1, 0, -1, f[0], f[1], 0, 0)
                    error, vig, x, y, z, l, m, n, l2, m2, n2, intensity = rayTraceData
                    dst_x.append(x)
                    dst_y.append(y)
                    orders.append(o.m)
                    wavelength.append(w)
                # plt.figure()
                # plt.scatter(dst_x, dst_y)
                # plt.show()

            # ppp.append(np.array(self.do_spectral_format(nPerOrder=nPerOrder, FSRonly=False, hx=f[0], hy=f[1])))
        # ppp = np.array(ppp)
        dst_x = np.array(dst_x)
        dst_y = np.array(dst_y)
        dst = np.vstack((dst_x, dst_y))
        dst /= ((self.CCD.pixelSize) / 1000.)
        dst += self.CCD.Nx / 2
        dst = dst.reshape(2, len(dst[0]) / len(norm_field), len(norm_field)).transpose((1, 2, 0))
        orders = np.array(orders)
        wavelength = np.array(wavelength)

        orders = orders.reshape((len(orders) / len(norm_field), len(norm_field)))
        wavelength = wavelength.reshape((len(wavelength) / len(norm_field), len(norm_field)))

        affine_matrices = {}
        transformations = {}

        for order, wavel, p in zip(orders, wavelength, dst):
            params = tf.estimate_transform('affine', src, p)
            if affine_matrices.has_key(order[0]):
                affine_matrices[order[0]].update({wavel[0]: np.array(
                    [params.rotation, params.scale[0], params.scale[1], params.shear, params.translation[0],
                     params.translation[1]])})
            else:
                affine_matrices[order[0]] = {wavel[0]: np.array(
                    [params.rotation, params.scale[0], params.scale[1], params.shear, params.translation[0],
                     params.translation[1]])}

        res['matrices'] = affine_matrices
        return res

    def walk_trough_configs(self, nWl=7, nPerSpot=5001, hx=0., hy=0.):
        actC, nC, operandC = self.ln.zGetConfig()
        for i in range(1, nC + 1):
            for j in range(1, nWl + 1):
                self.ln.zSetConfig(i)
                wl = self.ln.zGetWave(j).wavelength
                print(wl)

                rd_in, rd_out = self.trace(j, hx=hx, hy=hy, N=nPerSpot)
                o = self.ln.zGetSurfaceParameter(self.zmx_nsurf, 2)

                self.spots.append(Spot(wl, o, i - 1, rd_in, rd_out))

    def do_spots(self, nPerOrder=5, nOrders=5, FSRonly=True, nPerSpot=5001, hx=0, hy=0, everyNthOrder=5):
        n = everyNthOrder
        for o in list(self.Orders.values()):
            if n < everyNthOrder:
                n += 1
            else:
                print('Trace order', o.m)
                self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, o.m)
                if FSRonly:
                    wl = np.linspace(o.minFSRwl, o.maxFSRwl, nPerOrder)
                else:
                    wl = np.linspace(o.minWL, o.maxWL, nPerOrder)
                for i, w in enumerate(wl):
                    self.ln.zSetWave(1, w, 1.)
                    rd_in, rd_out = self.trace(hx=hx, hy=hy, N=nPerSpot)
                    self.spots.append(Spot(w, o.m, i, rd_in, rd_out))
                n -= everyNthOrder

    def do_spot_diagrams(self, order='all', nPerOrder=5, field=0):
        if order == 'all':
            for o in self.tracing:
                if o[0] <= self.maxord and o[0] >= self.minord:
                    print(("Trace order...", o[0]))
                    self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, o[0])
                    wl = np.linspace(o[1], o[2], nPerOrder)
                    for i, w in enumerate(wl):
                        self.ln.zSetWave(1, w, 1.)
                        asdf = self.trace_rays(1, field)
                        a = np.array(asdf, dtype=at.DdeArrayData)
                        wl = self.ln.zGetWave(self.ln.zGetPrimaryWave()).wavelength

                        vig = a['vigcode'][1:]
                        err = a['error'][1:]
                        vig = np.logical_and(vig, err)
                        index = np.logical_and(vig < 1, idx)
                        if np.max(index) > 0:
                            self.rays.append([a['x'][index], a['y'][index]])
                            self.wls.append(wl)

    def saveSpots(self, filename='spots.pkl'):
        print('save spots')
        pickle.dump(self.spots, open(self.savePath + filename, "wb"))

    def loadSpots(self, filename='spots.pkl'):
        self.spots = pickle.load(open(self.savePath + filename))

    def do_tracing(self, order='all', n=1000):
        if order == 'all':
            for o in self.tracing:
                print(("Trace order...", o[0]))
                self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, o[0])
                array = self.file.create_array(self.file.root, 'order' + str(o[0]), atom=np.array([3.]),
                                               shape=(2 * 4 * n,))
                wlarray = self.file.create_array(self.file.root, 'wl_order' + str(o[0]), atom=np.array([3.]),
                                                 shape=(n,))
                wl = np.linspace(o[1], o[2], n)

                for i, w in enumerate(wl):
                    self.ln.zSetWave(1, w, 1.)
                    xy = self.ln.zGetTrace(1, 0, -1, -1, -1, 0, 0)
                    array[i * 4 * 2] = xy[2]
                    array[i * 4 * 2 + 1] = xy[3]

                    xy = self.ln.zGetTrace(1, 0, -1, 1, -1, 0, 0)
                    array[i * 4 * 2 + 2] = xy[2]
                    array[i * 4 * 2 + 3] = xy[3]

                    xy = self.ln.zGetTrace(1, 0, -1, 1, 1, 0, 0)
                    array[i * 4 * 2 + 4] = xy[2]
                    array[i * 4 * 2 + 5] = xy[3]

                    xy = self.ln.zGetTrace(1, 0, -1, -1, 1, 0, 0)
                    array[i * 4 * 2 + 6] = xy[2]
                    array[i * 4 * 2 + 7] = xy[3]

                    wlarray[i] = w
                self.file.flush()
            self.file.close()

        else:
            self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, self.tracing[0][0])

            array = self.file.create_array(self.file.root, 'order' + str(self.tracing[0][0]), atom=np.array([3.]),
                                           shape=(2 * 4 * n,))
            wl = np.linspace(self.tracing[0][1], self.tracing[0][2], n)

            for i, w in enumerate(wl):
                self.ln.zSetWave(1, w, 1.)
                xy = self.ln.zGetTrace(1, 0, -1, -1, -1, 0, 0)
                array[i * 4 * 2] = xy[2]
                array[i * 4 * 2 + 1] = xy[3]

                xy = self.ln.zGetTrace(1, 0, -1, 1, -1, 0, 0)
                array[i * 4 * 2 + 2] = xy[2]
                array[i * 4 * 2 + 3] = xy[3]

                xy = self.ln.zGetTrace(1, 0, -1, 1, 1, 0, 0)
                array[i * 4 * 2 + 4] = xy[2]
                array[i * 4 * 2 + 5] = xy[3]

                xy = self.ln.zGetTrace(1, 0, -1, -1, 1, 0, 0)
                array[i * 4 * 2 + 6] = xy[2]
                array[i * 4 * 2 + 7] = xy[3]

            self.file.close()

    def setFile(self, name='MaroonXblue.h5', mode='w'):
        self.file = tables.open_file(name, mode=mode)

    def wavelength_to_order(self, wl):
        """
        Returns the order in which the wavelength appears.
        Returns empty list if wavelength is outside the spectral range.
        Returns a list of tuples, with the order number and a string indicating whether it is within FSR or not.
        :param wl: wavelength [micron]
        :return: list of tuples (order number, 'FSR'/'CCD')
        """
        res = []
        for o in list(self.Orders.values()):
            if o.inFSR(wl):
                res.append((o.m, 'FSR'))
            elif o.inOrder(wl):
                res.append((o.m, 'CCD'))
        return res

    def calc_wl(self):
        print('Calc wavelength')

        def find_lmin(order, dwl=0.0001):
            wl = self.ln.zGetWave(1)[0]
            vig = False
            wlmin = wl
            while vig < 1:
                wl = wl - dwl
                self.ln.zSetWave(1, wl, 1.)
                xy = self.ln.zGetTrace(1, 0, -1, 0, 0, 0, 0)
                vig = np.logical_or(xy[1], xy[0])
            else:
                print('vignetting at surface ', xy[1], self.order, wl)
                wlmin = wl
                xmin = xy[2]
                ymin = xy[3]
            self.x.append(xmin)
            self.y.append(ymin)
            return wlmin, xmin, ymin

        def find_lmax(order, dwl=0.0001):
            wl = self.ln.zGetWave(1)[0]
            vig = False
            wlmin = wl
            while vig < 1:
                wl = wl + dwl
                self.ln.zSetWave(1, wl, 1.)
                xy = self.ln.zGetTrace(1, 0, -1, 0, 0, 0, 0)
                vig = np.logical_or(xy[1], xy[0])
            else:
                print('vignetting at surface ', xy[1], self.order, wl)
                wlmin = wl
                xmin = xy[2]
                ymin = xy[3]
            self.x.append(xmin)
            self.y.append(ymin)
            return wlmin, xmin, ymin

        gamma_rad = np.deg2rad(self.gamma)
        blaze_rad = np.deg2rad(self.blaze)
        theta_rad = np.deg2rad(self.theta)
        self.grp = 1000. / self.grmm
        alpha = blaze_rad + theta_rad
        beta = blaze_rad - theta_rad
        c0 = self.grp * np.cos(gamma_rad)
        c1 = c0 * (np.sin(alpha) + np.sin(beta))
        c2 = c0 * np.cos(beta)
        c3 = self.grp * np.cos(blaze_rad) * (1. - np.tan(self.theta) * np.tan(blaze_rad))
        print(self.order + 1, c1 / (self.order + 1))
        self.ln.zSetWave(0, 1, 1)
        self.ln.zPushLens()

        vig = False
        # find max order
        o_working = self.order
        print('find max order --------------------')
        while vig < 1 and abs(self.order) < abs(self.maxord):
            if self.order > 0:
                self.order += 1
            else:
                self.order -= 1

            blazeWL = abs(c1 / self.order)
            print('Order: ', self.order, 'Blaze wl: ', blazeWL)
            self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, self.order)
            self.ln.zSetWave(1, blazeWL, 1.)
            self.ln.zGetUpdate()
            self.ln.zPushLens()
            xy = self.ln.zGetTrace(1, 0, -1, 0, 0, 0, 0)
            vig = np.logical_or(xy[1], xy[0])
            if vig < 1:
                self.x.append(xy[2])
                self.y.append(xy[3])
                self.orders.append(self.order)
                self.ln.zSetWave(1, blazeWL, 1.)
                self.ln.zPushLens()
                wmax = find_lmax(self.order)[0]

                self.ln.zSetWave(1, blazeWL, 1.)
                self.ln.zPushLens()

                wmin = find_lmin(self.order)[0]
                print("Order added ", self.order, wmin, wmax, blazeWL)
                self.Orders[self.order] = Order(self.order, blazeWL, wmin, wmax,
                                                blazeWL - blazeWL / self.order / 2.,
                                                blazeWL + blazeWL / self.order / 2.)
        # find min order
        vig = False
        self.order = o_working + 1
        print('find min order')
        while vig < 1 and abs(self.order) > abs(self.minord):
            print('test order', self.order, self.minord)
            if self.order > 0:
                self.order -= 1
            else:
                self.order += 1

            blazeWL = abs(c1 / self.order)
            self.ln.zSetSurfaceParameter(self.zmx_nsurf, 2, self.order)
            self.ln.zSetWave(1, blazeWL, 1.)
            self.ln.zPushLens()
            xy = self.ln.zGetTrace(1, 0, -1, 0, 0, 0, 0)
            vig = np.logical_or(xy[1], xy[0])
            if vig < 1:
                print('ok')
                self.orders.append(self.order)
                self.x.append(xy[2])
                self.y.append(xy[3])
                self.ln.zSetWave(1, blazeWL, 1.)
                self.ln.zPushLens()
                wmax = find_lmax(self.order)[0]

                self.ln.zSetWave(1, blazeWL, 1.)
                self.ln.zPushLens()
                wmin = find_lmin(self.order)[0]

                self.Orders[self.order] = Order(self.order, blazeWL, wmin, wmax,
                                                blazeWL - blazeWL / self.order / 2.,
                                                blazeWL + blazeWL / self.order / 2.)

    def spots_on_CCD(self):
        plt.figure()
        for s in self.spots:
            plt.scatter(s.x, s.y)
        plt.show()

    def EE_map(self, direction='r', plotSpots=True, zoom=150, save='', vmax=15., vmin=0., hx=0, hy=0, showplot=False,
               EE_ratio=80., additional_spots=[]):
        """
        generates encircled energy map from traced spots.

        :param direction: 'r', 'x' or 'y'
        :param plotSpots: plots spot diagramms as an overlay
        :param zoom: zoom of the individual spot diagrams
        :return:
        """
        print('EE map')

        fig, ax = plt.subplots()

        X = []
        Y = []
        R = []

        for s in self.spots:
            if np.mean(s.hx) == hx:
                if np.mean(s.hy) == hy:
                    X.append(s.barycenter['x'])
                    Y.append(s.barycenter['y'])
                    R.append(s.EE_radius(direction=direction, EE=EE_ratio))

            if plotSpots:
                if np.mean(s.hx) == hx:
                    if np.mean(s.hy) == hy:
                        ax.scatter(s.barycenter['x'] + zoom * s.xy_c['x'], -s.barycenter['y'] + zoom * s.xy_c['y'],
                                   s=.2, facecolor='black', lw=0)

        X = np.array(X)
        Y = np.array(Y)
        R = np.array(R)

        xi = np.linspace(-self.CCD.Wx / 2., self.CCD.Wx / 2., 101)
        yi = np.linspace(-self.CCD.Wy / 2., self.CCD.Wy / 2., 101)
        zi = griddata((X, Y), R, (xi[None, :], yi[:, None]), method='linear')

        ax.set_xlim((np.min(xi), np.max(xi)))
        ax.set_ylim((np.min(yi), np.max(yi)))
        ax.set_xlabel('Detector x [mm]')
        ax.set_ylabel('Detector y [mm]')

        im = ax.imshow(zi, interpolation='nearest', extent=[np.min(xi), np.max(xi), np.min(yi), np.max(yi)], vmin=vmin,
                       vmax=vmax)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cb = plt.colorbar(im, cax=cax)
        cb.set_label('EE' + str(EE_ratio) + ' radius [micron]')

        # for s in additional_spots:
        # ax.plot(additional_spots[:,0], additional_spots[:,1], 'w+')
        plt.tight_layout()
        if not save == '':
            plt.savefig(save, dpi=300)
        if showplot:
            plt.show()

    def spot_map(self):
        n_spots = len(self.spots)
        orders = []
        for s in self.spots:
            orders.append(s.order)
        unique_orders = np.unique(np.array(orders))
        n_orders = len(unique_orders)

        print('number of spots:', n_spots)
        print('number of orders:', n_orders)
        print('spot per order', n_spots / n_orders)

        fig, axarr = plt.subplots(n_orders, n_spots / n_orders, figsize=(n_spots / n_orders, n_orders + 1))
        for s in self.spots:
            j = np.where(unique_orders == s.order)[0][0]
            axarr[j, s.i].scatter(s.x, s.y, s=1, alpha=0.05, facecolor='blue', lw=0)
            axarr[j, s.i].set_xticklabels([])
            axarr[j, s.i].set_yticklabels([])
            axarr[j, s.i].axis('equal')
            axarr[j, s.i].xaxis.set_visible(False)
            axarr[j, s.i].yaxis.set_visible(False)
            axarr[j, s.i].axis('off')

            axarr[j, s.i].annotate('{:.4f}'.format(s.wavelength), xy=(
                s.barycenter['x'] - self.CCD.pixelSize / 2000., s.barycenter['y'] + self.CCD.pixelSize / 2000.),
                                   fontsize=8)
            axarr[j, s.i].add_patch(patches.Polygon(
                [[s.barycenter['x'] - self.CCD.pixelSize / 1000., s.barycenter['y'] + self.CCD.pixelSize / 1000.],
                 [s.barycenter['x'] + self.CCD.pixelSize / 1000., s.barycenter['y'] + self.CCD.pixelSize / 1000.],
                 [s.barycenter['x'] + self.CCD.pixelSize / 1000., s.barycenter['y'] - self.CCD.pixelSize / 1000.],
                 [s.barycenter['x'] - self.CCD.pixelSize / 1000., s.barycenter['y'] - self.CCD.pixelSize / 1000.]],
                fill=False))
            axarr[j, s.i].set_xlim((s.barycenter['x'] - self.CCD.pixelSize / 2000. * 1.15,
                                    s.barycenter['x'] + self.CCD.pixelSize / 2000. * 1.15))
            axarr[j, s.i].set_ylim((s.barycenter['y'] - self.CCD.pixelSize / 2000. * 1.15,
                                    s.barycenter['y'] + self.CCD.pixelSize / 2000. * 1.15))

        plt.show()

    def differential_FF_effects(self, path_image1, path_image2, xy1=None, xy2=None, r1=None, r2=None, k=1, plot=False,
                                saveplot=False):

        if path_image1.lower().endswith('.fit') or path_image1.lower().endswith('.fits'):
            weight_image1 = pyfits.open(path_image1)[0].data
        else:
            weight_image1 = np.array(Image.open(path_image1).convert('L'))

        if xy1 is None:
            xy1 = {}
            dims = np.shape(weight_image1)
            xy1['y'] = dims[0] / 2.
            xy1['x'] = dims[1] / 2.

        if r1 is None:
            r1 = np.shape(weight_image1)[0] / 2.

        weight_image1 = weight_image1[xy1['y'] - np.ceil(r1):xy1['y'] + np.ceil(r1),
                        xy1['x'] - np.ceil(r1):xy1['x'] + np.ceil(r1)]
        # normalized x and y coordinates (correspond to Px and Py in ZEMAX)
        xr = yr = np.arange(-1., 1., 1. / r1)

        # interpolation function
        f1 = interpolate.RectBivariateSpline(xr, yr, weight_image1)

        if path_image2.lower().endswith('.fit') or path_image2.lower().endswith('.fits'):
            weight_image2 = pyfits.open(path_image2)[0].data
        else:
            weight_image2 = np.array(Image.open(path_image2).convert('L'))

        if xy2 is None:
            xy2 = {}
            dims = np.shape(weight_image2)
            xy2['y'] = dims[0] / 2.
            xy2['x'] = dims[1] / 2.

        if r2 is None:
            r2 = np.shape(weight_image2)[0] / 2.

        weight_image2 = weight_image2[xy2['y'] - np.ceil(r2):xy2['y'] + np.ceil(r2),
                        xy2['x'] - np.ceil(r2):xy2['x'] + np.ceil(r2)]
        # normalized x and y coordinates (correspond to Px and Py in ZEMAX)
        xr2 = yr2 = np.arange(-1., 1., 1. / r2)

        # interpolation function
        f2 = interpolate.RectBivariateSpline(xr2, yr2, weight_image2)

        if plot or saveplot:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        xMean = []
        yMean = []
        diffxMean = []
        diffyMean = []
        for s in self.spots:
            b1 = s.calc_weighted_barycenter(f=f1)
            xMean.append(b1['x'])
            yMean.append(b1['y'])

            b2 = s.calc_weighted_barycenter(f=f2)

            diff = {'x': (b1['x'] - b2['x']) * 1000000. * 0.106,
                    'y': (b1['y'] - b2['y']) * 1000000. * 0.106}

            diffxMean.append(diff['x'])
            diffyMean.append(diff['y'])

            if plot or saveplot:
                # plt.scatter(b1['x'], b1['y'])
                arr = patches.FancyArrow(b1['x'], b1['y'], k * diff['x'], k * diff['y'], head_width=0.25, width=0.05,
                                         fc='k', aa=True)
                ax.add_patch(arr)
                # plt.scatter(b2['x'], b2['y'])

        xMean = np.array(xMean)
        yMean = np.array(yMean)
        diffxMean = np.array(diffxMean)
        diffyMean = np.array(diffyMean)

        if plot or saveplot:
            if self.CCD.dispersionDirection == 'x':
                zi = ml.griddata(xMean, yMean, diffyMean, self.CCD.xi, self.CCD.yi, interp='linear')
            else:
                zi = ml.griddata(xMean, yMean, diffxMean, self.CCD.xi, self.CCD.yi, interp='linear')

            img = ax.imshow(zi, interpolation='none', extent=self.CCD.extent, aspect='equal', origin='lower')
            cbar = fig.colorbar(img)
            cbar.ax.set_ylabel('shift in dispersion direction [m/s]')
        if plot:
            plt.show()
        if saveplot:
            plt.savefig(path_image2[:-3] + '_FF_rv.png')

        return diffxMean, diffyMean

    def calc_echellogram(self, nPerOrder=15, FSRonly=False):
        s00 = np.array(self.do_spectral_format(nPerOrder=nPerOrder, FSRonly=FSRonly))
        plt.figure()
        plt.plot(s00[:, 2], s00[:, 3], 'g+')
        plt.show()


if __name__ == "__main__":
    pass

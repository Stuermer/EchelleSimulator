# this script shows how to generate a spectrograph model file from the Zemax model of a spectrograph using PyZDDE
import PyEchelle
import pyzdde.zdde as pyz

# Create pyzdde link to ZEMAX
ln = pyz.createLink()

# Create echelle spectrograph
spectrograph = PyEchelle.Echelle(ln, 'Spectrograph')

# Analyze ZEMAX file to extract e.g. blaze and gamma angles
spectrograph.analyseZemaxFile(echellename='Echelle', blazename='Blaze', gammaname='Gamma')

# Define minimum and maximum order
spectrograph.minord = 77
spectrograph.maxord = 159

# Define CCD including pixel size
spectrograph.setCCD(PyEchelle.CCD(4096, 4096, 15, name='CCD'))

# Calculate wavelength bounds for all orders
spectrograph.calc_wl()

# Optional: Save/Load orders to avoid recalculating wavelength bounds next time
# spectrograph.saveOrders()
# spectrograph.loadOrders()

# Generate spectrograph model file for Echelle++
filename = 'model_spectrograph.hdf'
PyEchelle.save_spectrograph_info_to_hdf(filename, spectrograph)
PyEchelle.save_CCD_info_to_hdf(filename, spectrograph.CCD)

# Calculate affine transformation matrices across each order for fiber 1
att = spectrograph.do_affine_transformation_calculation(50, fw=50., fh=50.)
PyEchelle.save_transformation_to_hdf(filename, att, 1)

# Calculate PSFs across each order for fiber 1
psf = spectrograph.get_psfs(15, 1, [0, 0])
PyEchelle.save_psfs_to_hdf(filename, psf, 1)

# Close ZEMAX link
ln.close()

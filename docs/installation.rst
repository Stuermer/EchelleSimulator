Installation
============

There are two ways to install Echelle++: Building from source or using the docker image.

In principle we recommend building the source code, for optimal performance. However, we happily provide a docker image for a hassle-free installation.

Prerequisites
^^^^^^^^^^^^^
Echelle++ requires some 3rd party packages listed in below. Please install them for your platform.

 * GCC>4.9, capable of handling C++11 syntax (MSVC is untested, but should work)
 * `HDF5 <https://www.hdfgroup.org/hdf5/>`_
 * `CFITSIO <https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html>`_
 * `CCFits <http://heasarc.gsfc.nasa.gov/fitsio/ccfits/>`_
 * `Curl <https://curl.haxx.se/libcurl/>`_


Installation on ubuntu 18.04 (or similar)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install required libraries::
    
    sudo apt-get install libhdf5-dev libccfits-dev libcurl4-openssl-dev git git-lfs cmake build-essential

Build Echelle++:

.. code-block:: bash
    
    git clone https://github.com/Stuermer/EchelleSimulator.git
    cd EchelleSimulator
    git-lfs pull
    mkdir build
    cd build
    cmake ../
    make

This will compile the **Debug** built of Echelle++. This is good for developing / testing. But if you want full speed and multiprocessing change the 2nd last line to:

.. code-block:: bash

    cmake ../ -DCMAKE_BUILD_TYPE=Release

Other available build types are *RelNoParallel* and *RelWithDebInfo*. See *CMakeList.txt* for further infos.

After Echelle++ is successfully built, use the *echellesimulator* executable and start simulating spectra !

Using docker
^^^^^^^^^^^^

Install docker for your platform.

Run:

.. code-block:: bash

    docker run -v /path/to/output_directory:/home/simulations stuermer/echellesimulator

to download latest version of Echelle++ and run it with the desired command line arguments.

Note: **/path/to/output_directory** is the absolute path to your local directory where the resulting simulations are saved.
**/home/simulations** is the mount point of that directory inside the docker container. Please don't change that.

Further notes:
On Unix and MacOS **/path/to/output_directory** looks just like a regular path, e.g.:

.. code-block:: bash

    /path/to/output_directory = /home/USERNAME/simulations

For Windows, the path requires some special formatting and should look like this:

.. code-block:: bash

    /path/to/output_directory = ///c/Users/USERNAME/simulations

You might also want to check whether the docker-machine has enough memory and CPU's assigned to it.
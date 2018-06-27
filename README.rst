Installation
============

Scafacos
--------
To install and build the Scafacos libraries with dipole support, the following steps need to be taken (This will install to $HOME/lib and $HOME/include):::
    
    git clone git://github.com/scafacos/scafacos --branch dipoles
    cd scafacos
    git submodule init
    git submodule update
    ./bootstrap
    ./configure --with-internal-pfft --with-internal-pnfft --enable-fcs-solvers=direct,pnfft,p2nfft,p3m --prefix=$HOME --enable-shared --disable-fcs-fortran --enable-fcs-dipoles
    make -j8
    make install

Espresso
--------
Espresso is built as described at http://espressomd.org/html/doc/installation.html
However, two additional steps need to be tkaen.
* Before running CMake, the PKG_CONFIG_PATH environment variable needs to be adapted to include the path to Scafacos'  configuration:::

    export PKG_CONFIG_PATH=$HOME/lib/pkg_config:PKG_CONFIG_PATH

* Before running Espresso, one needs to ensure that Scafacos is found in the library path:::
    
    export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH

Files in this distribution
==========================

Files independent of a specific simulation project:

  * add-reference-forces-torques.py: Adds reference forces and torques obtained by direct summation to a particle configuration
  * tune.py: Tunes the P2nfft method parameters to obtain the fastest computation providing a given accuracy. The tuning is based on a particle configuration with reference forces and torques
  * get-accuracy.py: Prints out the accuracy provided by a P2nfft parameter set for a given particle configuration which includes reference forces and torques
  * p2nfft_common.py: Some shared routines

Files pertaining to the demonstration simulation:

  * model.py: The base class for the simulation model
  * gen-system.py: Generates and thermalizes an initial simulation system
  * run.py: Runs the simulation
  * run-test.sh: Shell script that runs a brief test simulation (aprox. 5 mins)
  * run-full.sh: Runs the production simulation as shown in the article




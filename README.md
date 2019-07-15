# r-ncdf4-build-opendap-windows
A windows build of the ncdf4 R package, that has opendap capability

The R package ncdf4 is a fantastic tool to read and handle netcdf files. But... on Windows you may have experienced that it is not possible to read files made available on the web as opendap files (e.g. thredds servers). It does work on Linux and Mac, though, but compilation of these functions for Windows was a hard nut to crack.
This build has resolved the problem. The package is available as a zip file for installation.

The current version is compatible with R 3.5.1.

For the user who wants to install and use ncdf4 on opendap files in Windows, only a single file is of importance: ncd4_1.15.zip. This contains the compiled package. All other files document the procedure of compilation. More details are given in procedure.rtf. There should normally be no need to go through the compilation process.

Instructions: download the ncdf4_1.15.zip file. Install the package ncdf4, choose installing from local archive, select the downloaded zip file, install, use.

Enjoy!

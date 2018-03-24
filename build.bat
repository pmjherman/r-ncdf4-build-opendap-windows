path = D:\herman\Documents\R\R-3.3.0\bin;D:\herman\Rtools\bin;D:\herman\Rtools\mingw_32\bin;C:\msys64\mingw64\bin;C:\Program Files\Common Files\Microsoft Shared\Microsoft Online Services;C:\Program Files (x86)\Common Files\Microsoft Shared\Microsoft Online Services;C:\WINDOWS\system32;C:\WINDOWS;C:\WINDOWS\System32\Wbem;C:\WINDOWS\System32\WindowsPowerShell\v1.0\;C:\Program Files (x86)\Microsoft Application Virtualization Client;C:\Program Files (x86)\Skype\Phone\;C:\Program Files\MiKTeX 2.9\miktex\bin\x64\;C:\Program Files\TortoiseSVN\bin;C:\Program Files\PostgreSQL\9.4\bin;C:\Program Files\MATLAB\MATLAB Compiler Runtime\v82\runtime\win64;D:\herman\DflowFMlibs\current
R CMD REMOVE ncdf4
tar cvf ncdf4_1.16.tar ncdf4 
gzip -v ncdf4_1.16.tar
set BINPREF="d:/herman/Rtools/mingw_64/bin/"
set NCDF_DIR="D:/herman/Rncdf4/"
R CMD INSTALL --build --no-multiarch D:/herman/Rncdf4/ncdf4_1.16.tar.gz
path = C:\Rtools\bin;"c:\Program Files\R\R-3.5.3\bin\";C:\WINDOWS\system32;C:\WINDOWS;C:\WINDOWS\System32\Wbem;C:\WINDOWS\System32\WindowsPowerShell\v1.0\;C:\WINDOWS\System32\OpenSSH\;C:\Program Files\MiKTeX 2.9\miktex\bin\x64\;C:\Program Files\Git\cmd;C:\Program Files\TortoiseSVN\bin;C:\Program Files (x86)\Calibre2\;C:\Users\herman\AppData\Local\Microsoft\WindowsApps;C:\Users\herman\AppData\Local\Programs\Microsoft VS Code\bin;C:\Users\herman\AppData\Local\Continuum\anaconda2\;C:\Users\herman\AppData\Local\Continuum\anaconda2\Scripts\;
R CMD REMOVE ncdf4
tar cvf ncdf4_1.16.tar ncdf4 
gzip -v ncdf4_1.16.tar
set BINPREF="C:/Rtools/mingw_64/bin/"
set NCDF_DIR=C:/Rncdf4
R CMD INSTALL --build --no-multiarch C:/Rncdf4/ncdf4_1.16.tar.gz
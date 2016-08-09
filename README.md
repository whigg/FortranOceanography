# FortranOceanography
Utilities for oceanography

* _convert_SMAP_v1_SSS_L3_netcdf.f90_: convert the binary files obtained from [ftp://ftp.remss.com/smap/L3/v1.0/8day_running/](the REMSS FTP) into [http://www.unidata.ucar.edu/software/netcdf/](netCDF). 
To compile:
```bash
 gfortran -I/usr/include convert_SMAP_v1_SSS_L3_netcdf.f90 -L/usr/lib/ -lnetcdf -lnetcdff -o convert_SMAP_v1_SSS_L3_netcdf.a
```

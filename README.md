# FortranOceanography
Utilities for oceanography

* _convert_SMAP_v1_SSS_L3_netcdf.f90_: convert the binary files obtained from [the REMSS FTP](ftp://ftp.remss.com/smap/L3/v1.0/8day_running) into [netCDF](http://www.unidata.ucar.edu/software/netcdf/). 
To compile:
```bash
 gfortran -I/usr/include convert_SMAP_v1_SSS_L3_netcdf.f90 -L/usr/lib/ -lnetcdf -lnetcdff -o convert_SMAP_v1_SSS_L3_netcdf.a
```

<img src="https://cloud.githubusercontent.com/assets/11868914/17514539/6f048808-5e33-11e6-9f56-9ffa9b944d8c.png" width="500">

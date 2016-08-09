! FORTRAN 90 reader for RSS SMAP Salinity Version 1.0 (BETA) Level 3 files
! to netCDF (CF-compliant) format.
!
!
! Compilation: you need to have the netCDF library installed
!
!  gfortran -I/usr/include convert_SMAP_v1_SSS_L3_netcdf.f90 -L/usr/lib/  -lnetcdf -lnetcdff -o convert_SMAP_v1_SSS_L3_netcdf.a
! ICTS SOCIB, 2016
!
! -----------------------------------------------------------------------------------------

program convert_SMAP_v1_SSS_L3_netcdf
  use netcdf
  implicit none
  

  character(len=250)                  ::  filename='sss_smap_8day_running_2016_182_v1.0.dat' ! to be specified by user
  character(len=50)                   ::  outfile='sss_smap_8day_running_2016_182_v1.0.nc'
  integer(4), parameter               ::  mlon=1440, mlat=720, iu=3
  real(4), dimension(mlon,mlat)       ::  map_sss
  integer(4), dimension(mlon,mlat)    ::  map_num 
  integer, parameter                  ::  ndims=2
  real(4), dimension(mlon)            ::  lon
  real(4), dimension(mlat)            ::  lat
 
  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  integer :: ncid, varid_lon, varid_lat, varid_sss, varid_ssscount, dimids(ndims)
  integer :: lon_dimid, lat_dimid
  integer :: ilon, ilat  
  integer :: nf90_put_att  

   ! We also create some parameters for the variable attributes
  character(len=50), parameter        ::  MISSING_VALUE = 'missing_value'
  character(len=50), parameter        ::  STANDARD_NAME = 'standard_name'
  character(len=50), parameter        ::  LONG_NAME = 'long_name'
  character(len=50), parameter        ::  UNITS = 'units'
 
  ! Create longitude and latitude vectors
  do ilon=1,mlon
    lon(ilon)=0.0+(ilon-1)*360./mlon
  end do

  do ilat=1,mlat
    lat(ilat)=-90.+(ilat-1)*180./mlat
  end do 

  ! Read the fortran binary file
  open(unit=iu,form='unformatted',file=filename,action='read',access='stream')
  read(iu) map_num
  read(iu) map_sss
  close(iu)

  ! Create the netCDf file
  call check( nf90_create(outfile, NF90_CLOBBER, ncid) )

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  call check( nf90_def_dim(ncid, "lat", mlat, lat_dimid) )
  call check( nf90_def_dim(ncid, "lon", mlon, lon_dimid) )
  
  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  !dimids =  (/ y_dimid, x_dimid /)
  dimids(1) = lon_dimid
  dimids(2) = lat_dimid

  ! Define the variable. The type of the variable in this case is
  ! NF90_INT (4-byte integer).
  call check( nf90_def_var(ncid, "lat", NF90_REAL, lat_dimid, varid_lat) )
  call check( nf90_def_var(ncid, "lon", NF90_REAL, lon_dimid, varid_lon) )
  call check( nf90_def_var(ncid, "l4_sss", NF90_REAL, dimids, varid_sss) )
  call check( nf90_def_var(ncid, "l4_sss_count", NF90_INT, dimids, varid_ssscount) )
 
  ! Add attributes to the variables
  call check( nf90_put_att(ncid, varid_lon, STANDARD_NAME, 'longitude') )
  call check( nf90_put_att(ncid, varid_lat, STANDARD_NAME, 'latitude') )
  call check( nf90_put_att(ncid, varid_lon, LONG_NAME, 'Longitude coordinate') )
  call check( nf90_put_att(ncid, varid_lat, LONG_NAME, 'Latitude coordinate') )

  call check( nf90_put_att(ncid, varid_sss, "missing_value", real(-9999.0)) )
  call check( nf90_put_att(ncid, varid_sss, STANDARD_NAME, 'sea_water_salinity') )
  call check( nf90_put_att(ncid, varid_sss, LONG_NAME, 'Sea water salinity') )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )

  ! Write the pretend data to the file. Although netCDF supports
  ! reading and writing subsets of data, in this case we write all the
  ! data in one operation.
  call check( nf90_put_var(ncid, varid_sss, map_sss) )
  call check( nf90_put_var(ncid, varid_ssscount, map_num) )
  call check( nf90_put_var(ncid, varid_lon, lon) )
  call check( nf90_put_var(ncid, varid_lat, lat) )  
  
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )


contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end program convert_SMAP_v1_SSS_L3_netcdf

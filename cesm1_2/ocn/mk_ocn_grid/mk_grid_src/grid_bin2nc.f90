!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !ROUTINE: grid_bin2nc
! !INTERFACE:

program grid_bin2nc

! !DESCRIPTION:
! This program takes information from standard POP files and writes a
! netCDF file.  The netCDF file contains all the information needed by
! g3 (= gCubed, a.k.a. Global Grid Generator framework).
!
! The POP binary files use radians and centimeters.  The netCDF file
! uses degrees and meters.
!
! !REVISION HISTORY:

! !USES:


implicit none
include 'netcdf.inc'

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer, parameter ::         &
     char_len  = 80                     ,&
     log_kind  = kind(.true.)           ,&
     int_kind  = kind(1)                ,&
     i4        = selected_int_kind(6)   ,&
     i8        = selected_int_kind(13)  ,&
     r4        = selected_real_kind(6)  ,&
     r8        = selected_real_kind(13)

integer ::      &
      att_ival,     &! netCDF data type
      att_num,      &! the volatile number currently assigned to an attribute
      dimid, dimid1, dimid2,        &! netCDF dimension id
!D       dimids(NF90_MAX_VAR_DIMS), &
      dim1, dim2,   &! length of longitude, latitude dimensions
      iostatus,     &! status flag
      io_record_length, &
      itype,        &! netCDF data type
      length,       &! number of elements in a dimension
      i, j, k, n,   &
      ncid,         &! netCDF file id
      nsize,        &! size parameter returned by inquire function (BAD; inquire doesn't wk)
      num_atts,     &! number of global attributes
      num_dims,     &! number of dimensions
      num_vars,     &! number of variables
      varid,        &! netCDF variable id
      xtype

   logical(log_kind) :: &
      att_lval           ! temp space for logical attribute

   integer(i4), allocatable :: kmt(:,:), cellId_array(:,:), &
		inc_array(:,:)

   real(r4) :: att_rval       ! temp space for real attribute
   real(r4) :: a, b, c        ! temp variables

   real(r8), allocatable :: nc_array(:,:), binary_array(:,:), &
                 plot_array(:,:), elevation(:,:), dz(:), z(:)

   real(r8) ::     &
      att_dval       ! temp space for double attribute

   logical(log_kind) :: &
      attrib_error        ! error flag for reading attributes
   character(len=80) :: filename
!cas
    integer,dimension(2)::ostart,ocount


! Parameters(coded as variables)
real(r8) :: PI, DegToRad, MToCm, RadToDeg, CmToM
!D character(len=NF90_MAX_NAME) :: name

PI = 4.0 * atan(1.0)
DegToRad = PI/180.0
RadToDeg = 180.0/PI
MToCm = 100.0
CmToM = 0.01



! Begin executable code
print *, 'Enter length of longitude dimension'
read *, dim1
print *, 'Longitude entered ', dim1

print *, 'Enter length of latitude dimension'
read *, dim2
print *, 'Latitude entered ', dim2

iostatus = nf_noerr
print *, 'Enter netCDF output filename'
read *, filename
iostatus = nf_create(trim(filename), nf_clobber, ncid)
call check(iostatus)
print *, 'All output data will go to netCDF file ', trim(filename)
print *, 'NCID ', ncid

iostatus = nf_def_dim(ncid, 'longitude', dim1, dimid1)
call check(iostatus)
print *, 'DIMID 1 ', dimid1

iostatus = nf_def_dim(ncid, 'latitude', dim2+1, dimid2)
call check(iostatus)
print *, 'DIMID 2 ', dimid2

!D iostatus = nf90_Inquire(ncid=ncid, nDimensions = num_dims)
!D call check(iostatus)
!D print *, 'NUMBER OF DIMENSIONS ', num_dims
!D if(num_dims > 0) then
!D     do dimid = 1, num_dims
!D         iostatus = nf90_Inquire_Dimension(ncid=ncid, dimid=dimid, &
!D                                           name=name, len=length)
!D         print *, '      DIM ', dimid, ' NAME ', name, ' LENGTH ', length
!D     end do
!D end if

allocate(nc_array     (0:dim1-1, 0:dim2) )
allocate(inc_array     (0:dim1-1, 0:dim2) )
allocate(cellId_array (0:dim1-1, 0:dim2) )
allocate(binary_array (1:dim1,   1:dim2) )
allocate(plot_array   (1:dim1,   0:dim2) )
allocate(elevation (1:dim1,   1:dim2) )
allocate(kmt (1:dim1,   1:dim2) )

ostart =(/   1,  1/)
ocount =(/dim1, dim2+1/)

print *, 'Dimensions ', dimid1, ' ', dimid2
print *, dim1
print *, dim2

! Cell ID is defined only for backwards compatibility.  This variable
! will be removed someday.  No realistic data need be stored for it.
iostatus = nf_def_var(ncid, 'cellId', nf_int,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
print *,' cellID defined'

! Elevation is negative below sea level
iostatus = nf_def_var(ncid, 'elevation',nf_float, 2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
iostatus = nf_put_att_text(ncid, varid, 'units', len_trim('meters'),'meters')
call check(iostatus)
iostatus = nf_put_att_text(ncid, varid, 'long_name', len_trim('Topography'),'Topography')
call check(iostatus)
print *,' elevation defined'

! KMT
iostatus = nf_def_var(ncid, 'kmt', nf_float,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
print *,' kmt defined'

! These variables are described in the POP manuals
iostatus = nf_def_var(ncid, 'ULAT', nf_float,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
iostatus = nf_def_var(ncid, 'ULON', nf_float,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
iostatus = nf_def_var(ncid, 'HTN', nf_float,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
iostatus = nf_def_var(ncid, 'HTE', nf_float,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
iostatus = nf_def_var(ncid, 'HUS', nf_float,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
iostatus = nf_def_var(ncid, 'HUW', nf_float,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
iostatus = nf_def_var(ncid, 'ANGLE', nf_float,2, (/ dimid1, dimid2 /), varid)
call check(iostatus)
print *,' cell lengths and angles defined'

iostatus = nf_enddef(ncid)
call check(iostatus)

! cellId
cellId_array = 1
iostatus = nf_put_var_int(ncid, 1, cellId_array)
call check(iostatus)

! Depth profile normally in [centimeter]
print *, 'Enter depth profile filename'
read *, filename
print *, 'First Depth profile will be read from text file ', trim(filename)
open(unit=9, access='sequential', action='read', file=trim(filename),  &
                status='old')
print *, 'Depth profile will be read from text file ', trim(filename)

! The depth file is always small, so easier to count the layers than to
! prompt for them.
n = 0
do while(iostatus==0)
   read(unit=9, fmt=*, iostat=iostatus) a
   n = n + 1
end do
n = n - 1
!D print *, 'Number of layers ', n

! After counting the layers, go back to the top of the file and read
! the values.
rewind(unit=9)
allocate(dz(1:n))
allocate(z(1:n))
a = 0.0
do i = 1, n
   read(unit=9, fmt=*, iostat=iostatus) dz(i)
   a = a - 0.01*dz(i) ! Convert to meters, negative below sea level
   z(i) = a
!D    print *, 'Z(', i, ') = ', z(i)
end do

close(unit=9)

!inquire(iolength=io_record_length) kmt(1:dim1, 1:dim2)
!D print *, ' Record length for kmt file ', io_record_length
! inquire is not working anymore. 4apr13 nanr
io_record_length = (dim1*dim2*8)
print *, ' Calculated length for kmt record ', io_record_length

print *, 'Enter KMT filename'
read *, filename
print *, 'KMT filename ', trim(filename)
open(unit=10, access='direct', action='read', file=trim(filename), &
                recl=io_record_length, status='old')
print *, 'KMT will be read from binary file ', trim(filename)

! elevation
nc_array = 1.0
read(10, rec=1) kmt

close(unit=10)

! print *, 'kmt = ', kmt

! Set an arbitrarily large positive elevation to represent land.
! This value yields a reasonable spread of colors for the ocean.
elevation = 7000.0
do i = 1, n
   where(kmt == i) elevation = z(i)
end do


nc_array(0,        1:dim2) = elevation (dim1,     1:dim2)
nc_array(1:dim1-1, 1:dim2) = elevation (1:dim1-1, 1:dim2)
nc_array(0:dim1-1,      0) = elevation (1:dim1,        1)
! iostatus = nf_put_vara_real(ncid, 2, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 2, ostart,ocount,nc_array)
call check(iostatus)
! print *,' nc_array elevation', nc_array
! print *,'  elevation', elevation 
nc_array(0:dim1-1, 1:dim2) = real(kmt (1:dim1, 1:dim2))
nc_array(0:dim1-1,      0) = real(kmt (1:dim1,        1))
! iostatus = nf_put_vara_real(ncid, 3, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 3, ostart,ocount,nc_array)
call check(iostatus)
! print *,' nc_array kmt', nc_array 

! Horizontal grid points
! To display the horizontal grid correctly we need the coordinates
! of the points at the 'northeast' corner of every T-cell.  This is
! sufficient for the east-west plotting because the coordinates
! wrap around the globe.  Along the 'southern' edge, however, we
! need an extra row of coordinates.  These could be estimated from
! the HTE/HTN values of the cells, but we read them from a pre-existing
! binary file produced for the purpose of plotting.
! inquire(iolength=io_record_length) plot_array(1:dim1, 0:dim2)
!D print *, ' Record length for plot file ', io_record_length
print *, ' Record length for plot file ', io_record_length
! inquire is not working anymore. 4apr13 nanr
io_record_length = (dim1*(dim2+1)*8)
print *, ' Calc   length for plot file ', io_record_length

print *, 'Enter binary plot filename '
read *, filename
open(unit=10, access='direct', action='read', file=trim(filename), &
                recl=io_record_length, status='old')
print *, 'Grid points will be read from binary file ', trim(filename)

! ULAT
read(10, rec=1) plot_array
nc_array(0,          0:dim2) = RadToDeg * plot_array (dim1,     0:dim2)
nc_array(1:dim1-1,   0:dim2) = RadToDeg * plot_array (1:dim1-1, 0:dim2)
! iostatus = nf_put_vara_real(ncid, 4, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 4, ostart,ocount,nc_array)
call check(iostatus)
print *,' nc_array ulat', nc_array 

! ULON
read(10, rec=2) plot_array
nc_array(0,          0:dim2) = RadToDeg * plot_array (dim1,     0:dim2)
nc_array(1:dim1-1,   0:dim2) = RadToDeg * plot_array (1:dim1-1, 0:dim2)
! iostatus = nf_put_vara_real(ncid, 5, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 5, ostart,ocount,nc_array)
call check(iostatus)
print *,' nc_array ulon', nc_array 
close(unit=10)

!! inquire(iolength=io_record_length) binary_array(1:dim1, 1:dim2)
!D print *, ' Record length for binary input ', io_record_length
! inquire is not working anymore. 4apr13 nanr
io_record_length = (dim1*(dim2)*8)
print *, ' Record length for binary input ', io_record_length

print *, 'Enter binary grid filename'
read *, filename
open(unit=10, access='direct', action='read', file=trim(filename), &
                recl=io_record_length, status='old')
print *, 'Cell data will be read from binary file ', trim(filename)

! The cell data will all be padded with an extra 'southern' row
! (nc_array(0,:).  The data values used for padding are chosen to
! be inconspicuous when the raw data are displayed.
! HTN
read(10, rec=3) binary_array
nc_array(0,        1:dim2) = CmToM * binary_array (dim1,     1:dim2)
nc_array(1:dim1-1, 1:dim2) = CmToM * binary_array (1:dim1-1, 1:dim2)
nc_array(0:dim1-1,      0) = CmToM * binary_array (1:dim1,        1)
! iostatus = nf_put_vara_real(ncid, 6, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 6, ostart,ocount,nc_array)
call check(iostatus)
print *,' nc_array htn', nc_array 

! HTE
read(10, rec=4) binary_array
nc_array(0,        1:dim2) = CmToM * binary_array (dim1,     1:dim2)
nc_array(1:dim1-1, 1:dim2) = CmToM * binary_array (1:dim1-1, 1:dim2)
nc_array(0:dim1-1,      0) = CmToM * binary_array (1:dim1,        1)
! iostatus = nf_put_vara_real(ncid, 7, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 7, ostart,ocount,nc_array)
call check(iostatus)
print *,' nc_array hte', nc_array 

! HUS
read(10, rec=5) binary_array
nc_array(0,        1:dim2) = CmToM * binary_array (dim1,     1:dim2)
nc_array(1:dim1-1, 1:dim2) = CmToM * binary_array (1:dim1-1, 1:dim2)
nc_array(0:dim1-1,      0) = CmToM * binary_array (1:dim1,        1)
! iostatus = nf_put_vara_real(ncid, 8, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 8, ostart,ocount,nc_array)
call check(iostatus)
print *,' nc_array hus', nc_array 

! HUW
read(10, rec=6) binary_array
nc_array(0,        1:dim2) = CmToM * binary_array (dim1,     1:dim2)
nc_array(1:dim1-1, 1:dim2) = CmToM * binary_array (1:dim1-1, 1:dim2)
nc_array(0:dim1-1,      0) = CmToM * binary_array (1:dim1,        1)
! iostatus = nf_put_vara_real(ncid, 9, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 9, ostart,ocount,nc_array)
call check(iostatus)
print *,' nc_array huw', nc_array 

! ANGLE
read(10, rec=7) binary_array
nc_array(0,        1:dim2) = RadToDeg * binary_array (dim1,     1:dim2)
nc_array(1:dim1-1, 1:dim2) = RadToDeg * binary_array (1:dim1-1, 1:dim2)
nc_array(0:dim1-1,      0) = RadToDeg * binary_array (1:dim1,        1)
! iostatus = nf_put_vara_real(ncid, 10, ostart,ocount,nc_array)
iostatus = nf_put_vara_double(ncid, 10, ostart,ocount,nc_array)
call check(iostatus)
print *,' nc_array angle', nc_array 

close(unit=10)


!D iostatus = nf90_Inquire(ncid=ncid, nVariables = num_vars)
!D call check(iostatus)
!D print *, 'NUMBER OF VARIABLES ', num_vars
!D if(num_vars > 0) then
!D     do varid = 1, num_vars
!D         iostatus = nf90_Inquire_Variable(ncid=ncid, varid=varid, &
!D                                          name=name, xtype=xtype, &
!D                                          ndims=num_dims,            &
!D                                          dimids=dimids(:num_dims),  &
!D                                          nAtts=num_Atts)
!D         print *, '      VAR ', varid, ' NAME ', name, ' TYPE ', xtype, &
!D              ' N DIMS ', num_dims, ' N ATTS ', num_Atts
!D         do n = 1, num_dims
!D             print *, '                 DIMENSION ID ', n
!D         end do
!D         do att_num = 1, num_Atts
!D             iostatus = nf90_inq_attname(ncid=ncid, varid=varid, &
!D                                          attnum=att_num, name=name)
!D             print *, '                 ATTRIBUTE ', att_num, ' NAMED ', name
!D         end do
!D     end do
!D end if


iostatus = nf_close(ncid)
call check(iostatus)
!D print *, 'CLOSE', ncid

print *, 'DONE'
end program grid_bin2nc

!-----------------------------------------------------------------------
!EOC

 subroutine check(status)
 include 'netcdf.inc'
    integer, intent(in) :: status
    if(status /= nf_noerr) then
       print *, nf_strerror(status)
       stop 'ERROR in netcdf'
    end if
 end subroutine check

program gridkmt_nc2bin

!-------------------------------------------------------------------------------
! PURPOSE:
! o Reads the kmt field from a netCDF file generated using grid2nc
!   and writes it out as a standard POP binary kmt file.
!-------------------------------------------------------------------------------

use netcdf

implicit none

!--- local ---
   integer         :: io_record_length,nx,ny
   integer         :: did,fid,vid,rcode
   integer, allocatable :: kmt(:,:),kmtin(:,:)
   character(len=80) :: filein, fileout,dname

   print *, 'Enter netCDF filename containing kmt array:'
   read *, filein
   print *, 'Enter POP binary output file name:'
   read *, fileout

   rcode = nf90_open(filein,NF90_NOWRITE,fid)
   rcode = nf90_inq_dimid(fid, 'longitude' , did)
   rcode = nf90_Inquire_Dimension(fid, did, dname, nx  )
   rcode = nf90_inq_dimid(fid, 'latitude' , did)
   rcode = nf90_Inquire_Dimension(fid, did, dname, ny  )

   allocate(kmtin(nx,ny)) 
   allocate(kmt(nx,ny-1)) 

   rcode = nf90_inq_varid(fid,'kmt'  ,vid)
   rcode = nf90_get_var(fid,vid,kmtin)
   rcode = nf90_close(fid)

   kmt = kmtin(1:nx,2:ny)

   inquire(iolength=io_record_length) kmt(1:nx,1:ny-1)
   open(10, access='direct', file=trim(fileout), &
                recl=io_record_length, status='unknown')
   write(10, rec=1) kmt
   close(10)


end program gridkmt_nc2bin


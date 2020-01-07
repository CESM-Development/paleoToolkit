program topography_driver

   implicit none
   integer     :: nlat,nlng,nlev

   print *, " "
   print *, "enter grid dimensions: nlng,nlat,nlev "
   read(*,*) nlng,nlat,nlev
   print *, ">>>>>",nlng,nlat,nlev

   call topography(nlng,nlat,nlev)

end program topography_driver

!---------------------------------------------------------------------

subroutine topography(nlng,nlat,nlev)

   implicit none

   integer  :: nlng,nlat,nlev,n,iocean,jocean
   real*8     ::  ulat(nlng,0:nlat),ulng(nlng,0:nlat)
   real*8     ::  h(nlng,0:nlat)
   integer  ::  kmt(nlng,nlat)
   integer  ::  nxt, nyt
   real*8     ::  dz(nlev),zt(nlev)
   character*80  :: filename

   print *, " "
   print *, "enter name of direct access file containing"
   print *, "double-precision fields:"
   print *, "ulat(1:nlng,0:nlat),ulng(1:nlng,0:nlat)"
   read(*,'(a80)') filename
   print *, ">>>>>",filename

   open (53, file=filename, form='unformatted', status='old', & 
     & access='direct',action='read',recl=(nlat+1)*nlng*8)
   read (53,rec=1) ulat
   read (53,rec=2) ulng
   close(53)
   print *, "file read: ",filename
   !print *, "ulat = ", ulat
   !print *, "ulon = ", ulng

   !print *," before zcalc "
   call zcalc(nlev,dz,zt)

   !print *," after zcalc "

   call add_topo(nlng,nlat,nlev,ulng,ulat,h)

   print *, " "
   print *, "enter: "
   print *, "  0 to continue "
   print *, "  1 to output depth field H(i,j) to direct access file"
   read(*,*) n
   print *, ">>>>>",n
   if(n.eq.1) then
     print *, " "
     print *, "enter filename"
     read(*,'(a80)') filename
     print *, ">>>>>",filename
     open (54, file=filename, form='unformatted', status='unknown', &
  &   access='direct',action='write',recl=nlat*nlng*8)
     write(54,rec=1) h(:,1:nlat)
     close(54)
     print *, "file written: ",filename
   endif

   call depth_levels(nlng,nlat,nlev,dz,zt,h,kmt)

   print *, " "
   print *, "enter: "
   print *, "  0 to continue "
   print *, "  1 to remove isolated points at all levels from KMT"
   read(*,*) n
   print *, ">>>>>",n
   if (n.eq.1) call remove_isolated_points(nlng,nlat,nlev,kmt)

   print *, " "
   print *, "enter: "
   print *, "  0 to continue "
   print *, "  1 to remove disconnected bays that are not"
   print *, "    connected to a selected point"
   read(*,*) n
   print *, ">>>>>",n
   if (n.eq.1) then
2       print *, " "
        print *, " enter indices (iocean, jocean) "
        print *, " for a select point in the ocean"
        read(*,*) iocean,jocean
        print *, ">>>>>",iocean,jocean
        if (kmt(iocean,jocean).eq.0) then
          print *, " "
          print *, " that is a land point, try again "
          go to 2
        endif
        call remove_disconnected_bays(iocean,jocean,nlng,nlat,kmt)
   endif

   print *, " "
   print *, "enter: "
   print *, "  0,     to continue "
   print *, &
     & "  n > 0, to force KMT to have a minimum value of n at ocean pts"
   read(*,*) n
   print *, ">>>>>",n
   if(n.gt.0) then
     where(kmt.gt.0) kmt = max(kmt,n)
   endif

   print *, " "
   print *, "enter: "
   print *, "  0 to continue "
   print *, "  1 to output KMT field to direct access file"
   read(*,*) n
   print *, ">>>>>",n
   if(n.eq.1) then
     print *, " "
     print *, "enter filename"
     read(*,'(a80)') filename
     print *, ">>>>>",filename
     open (55, file=filename, form='unformatted', status='unknown', &
  &    access='direct',action='write',recl=nlat*nlng*4)
     write(55,rec=1) kmt
     close(55)
     print *, "file written: ",filename
   endif
    print *,kmt
    
end subroutine topography

!----------------------------------------------------------------------
 
subroutine zcalc(nlev,dz,zt)
 
!----------------------------------------------------------------------
!     Reads depth profile for number of levels in model grid
!     calculates depths from the profile
!----------------------------------------------------------------------
!      input:
!        nlev = No. of depths in model
!      output:
!        zt = depth to center of each level
!        dz = depth profile
!      external requirement:
!        'in_depths.dat' = ascii file of dz for nlev levels
!----------------------------------------------------------------------
 
   implicit none

   integer :: nlev,k
   real*8    :: dz(nlev),zt(nlev)
   character*80 :: filename
   character(len=*),parameter :: F00 = "(x,i3,x,f12.3,x,f12.3)"

!  read dz depth profile from file
   print *, " "
   print *, " enter name of ascii grid file with depth profile in cm"
   read(*,'(a80)') filename
   print *, ">>>>>",filename

   open(54,file=filename,status='old',form='formatted')
   do k = 1,nlev
      read (54,*) dz(k)
   enddo
!     print *,k,dz
   close (54)
 
   print *, 'file read: ',filename
!     calculate z from dz
 
   zt(0) = 0.
   zt(1)=.5*dz(1)
   do k=2,nlev
      zt(k)=zt(k-1) + 0.5*(dz(k)+dz(k-1))
   enddo
   print *, " "
   print *, "  k       dz(m)       zt(m)"
   do k=0,nlev
      write(*,F00) k,dz(k)*0.01,zt(k)*0.01
   enddo
 
end
 
!---------------------------------------------------------------------
!     subroutine linking the mesh and the topography
!---------------------------------------------------------------------
 
subroutine add_topo(nlng,nlat,nlev,ulng,ulat,h)

   implicit none
   include 'netcdf.inc'


   integer :: nlng,nlat,nlev,nxt,nyt
   real*8    :: ulat(nlng,0:nlat),ulng(nlng,0:nlat)
   real*8    :: h(nlng,0:nlat)
   integer, allocatable :: itopo(:,:)
   real, allocatable    :: topo(:,:)
   real*8    :: x(3,5),ullng(4),ullat(4),pi
   real*8    :: real_lat,real_lng,hh,hland,hocean,xp,yp,zp
   real*8    :: xv1,xv2,xpv,yv1,yv2,ypv,zv1,zv2,zpv,fac
   integer :: i,j,k,i1,i2,im1,j1,l,kp1,ln,la
   integer :: lngmax,lngmin,latmax,latmin,nbland,nbpocean,nbpland
   integer :: did,fid,vid,rcode,vtype,vdims,vnatts,vdimids(2)
   logical :: incell, land
   character*60 :: fname,vname

   character(len=*),parameter :: F00 = "(20i6)"

   pi = 3.141592653589793d0
   fac = nxt/360.

   print *,'nxt = ',nxt
 
!-----------------------------------------------------------------
!     reading topography...
!     The  integer values (depth at each point) are in meters, positive
!     above sealevel, neg. below.
!     The resolution should be at least 1/2degx1/2deg.
!     The data is assumed to have a top latitude row at 90N, and a first
!     longitude column at 0E.
!-----------------------------------------------------------------

   print *, " "
   print *, " reading topo data..."
   print *, " enter name of netcdf file with <topo> field:"
   read(*,'(a60)') fname

   rcode = nf_open(trim(fname),nf_nowrite,fid)
   rcode = nf_inq_varid(fid,'topo',vid)
   rcode = nf_inq_var(fid,vid,vname,vtype,vdims,vdimids,vnatts)
   rcode = nf_inq_dim(fid,vdimids(1),vname,nxt)
   rcode = nf_inq_dim(fid,vdimids(2),vname,nyt)
   print *, " dimension of <topo> field: ",nxt,nyt
   allocate(topo(nxt,nyt))
   allocate(itopo(nxt,nyt))
   rcode = nf_get_var_real(fid,vid,topo)
   rcode = nf_close(fid)
 
   do j=1,nyt
      itopo(1:nxt,nyt+1-j) = int(topo(1:nxt,j))
   enddo
 
   fac = nxt/360.
   print *, " "
   print *, " calculating depth field..."
 
!----------------------------------------------------------
!     is there land at T-points?
!----------------------------------------------------------
 
   nbland = 0
!  loop over mesh grid
   do j=1,nlat
      do i=1,nlng
         if (i.eq.1) then
            im1 = nlng
         else
            im1 = i-1
         endif

!     what are the coordonates of the 4 edges of the cell?
         do k=1,4
            if (k.eq.1.or.k.eq.2) then
               j1 = j-1
            else
               j1 = j
            endif
            if (k.eq.1.or.k.eq.4) then
               i1 =  im1
            else
               i1 = i
            endif
            x (1,k) = cos(ulng(i1,j1))*cos(ulat(i1,j1))
            x (2,k) = sin(ulng(i1,j1))*cos(ulat(i1,j1))
            x (3,k) = sin(ulat(i1,j1))
            ullng (k) = ulng(i1,j1)
            ullat (k) = ulat(i1,j1)
         enddo
 
 
!----------------------------------------------------------
!     what are the (long-lat) boundaries of the cell ?
!----------------------------------------------------------
 
         lngmax = int(max(ulng(i,j),ulng(im1,j), &
     &           ulng(i,j-1),ulng(im1,j-1))*180.*fac/pi)+1
         lngmin = int(min(ulng(i,j),ulng(im1,j), &
     &           ulng(i,j-1),ulng(im1,j-1))*180.*fac/pi)+1
         latmin = (90*fac - int(max(ulat(i,j),ulat(i,j-1), &
     &           ulat(im1,j),ulat(im1,j-1))*fac*180./pi)) + 1
         latmax = (90*fac - int(min(ulat(i,j),ulat(i,j-1), &
     &           ulat(im1,j),ulat(im1,j-1))*fac*180./pi)) + 1
         if (float(lngmax-lngmin)/float(nxt).gt.0.25) then
            l = lngmin
            lngmin = lngmax
            lngmax = l+nxt
         endif
 
!----------------------------------------------------------
!     calculation of depth
!----------------------------------------------------------
 
         hland = 0.
         hocean = 0.
         incell = .true.
         land = .false.
         nbpocean =0
         nbpland =0
         do j1=latmin,latmax
            do i2=lngmin,lngmax
 
               i1 = mod(i2,nxt)
               incell = .true.
               real_lng = float(i1-1)*2.*pi/float(nxt)
               real_lat = pi/2.-float(j1-1)*pi/float(nyt)
               xp = cos(real_lng)*cos(real_lat)
               yp = sin(real_lng)*cos(real_lat)
               zp = sin(real_lat)
 
               do k=1,4
                  kp1 = mod(k,4)+1
                  xv1 = x(1,k) - xp
                  yv1 = x(2,k) - yp
                  zv1 = x(3,k) - zp
                  xv2 = x(1,kp1) - xp
                  yv2 = x(2,kp1) - yp
                  zv2 = x(3,kp1) - zp
                  xpv = yv1*zv2 - zv1*yv2
                  ypv =-xv1*zv2 + zv1*xv2
                  zpv = xv1*yv2 - yv1*xv2
                  incell = incell.and.(xpv*xp+ypv*yp+zpv*zp.ge.0)
               enddo
 
               if (incell) goto 20
               goto 10
20               continue
               if (itopo(i1,j1).ge.0.0) then
                  nbpland  = nbpland +1
!     hland = hland + float(itopo(i1,j1))
               else
                  hocean = hocean + float(itopo(i1,j1))
                  nbpocean  = nbpocean +1
               endif
10               continue
            enddo
         enddo

	!print *,'Got this far'
 
         if (nbpland.gt.nbpocean) then
            hh =  0.0
         else
            if (nbpocean+nbpland.eq.0) then
               hh=0.0
               do k=1,4
                  ln = int(ullng(k)*180.*fac/pi)+1
                  la = (90*fac - int(ullat(k)*fac*180./pi)) + 1
                  hh = hh+ itopo(ln,la)
                  nbpocean = 1
               enddo
               hh = hh/4.
!                 print *,i,j,tlng(i,j)*180./pi,tlat(i,j)*180./pi,hh
            else
               hh = (hocean+hland)/float(nbpland+nbpocean)
            endif
         endif
         if (hh.ge.0.0) then
            nbland = nbland +1
         endif
         h(i,j) = -hh*100.0
         enddo
!        print *,'nb de points :',nbland,' pour j ',j
   enddo
 
end subroutine add_topo
 
!----------------------------------------------------------------------
 
subroutine depth_levels(nlng,nlat,nlev,dz,zt,h,kmt)
 
!-----------------------------------------------------------------------
!     find KMT field given H field
!
!     input:
!       H = depth field (cm)
!       Z = model depths (cm)
!       nlng, nlat, nlev = No. of longitudes, latitudes, depths in model
!
!     output:
!       KMT field
!-----------------------------------------------------------------------

   implicit none

   integer :: i,j,k,nlat,nlng,nlev
   real*8    :: h(nlng,0:nlat),dz(nlev),zt(nlev)
   integer :: kmt(nlng,nlat)
   real*8    :: hmin
 
!***********************************************************************

   print *, " "
   print *, " computing kmt field ..."
   print *, " "
   print *, "enter shallow water cuttoff (m)' "
   print *, "  (points with depths shallower than this cuttoff"
   print *, "   will be automatically set to land)"
   read(*,*) hmin
   print *, ">>>>>",hmin
!  hmin=5.
   hmin = hmin*100.   ! convert to cm
 
   do j = 1,nlat
      do i = 1,nlng
         if(h(i,j) .lt. hmin) h(i,j) = 0.0
         if((h(i,j) .ge. hmin) .and. (h(i,j) .lt. zt(1))) &
     &         h(i,j) = zt(1)
         kmt(i,j) = 0
         if(h(i,j) .lt. zt(1))  kmt(i,j) = 0
         if(h(i,j) .ge. zt(nlev)) kmt(i,j) = nlev
      enddo
   enddo
 
   do k = 1,nlev-1
      do j = 1,nlat
         do i = 1,nlng
            if((h(i,j) .ge. zt(k)) .and. (h(i,j) .lt. zt(k+1))) &
     &            kmt(i,j) = k
         enddo
      enddo
   enddo
 
   do i = 1,nlng
      kmt(i,1)    = 0
      kmt(i,nlat) = 0
   enddo
 
end subroutine depth_levels
 
!---------------------------------------------------------------------

subroutine remove_isolated_points(nlng,nlat,nlev,kmt)

!---------------------------------------------------------------------
!     removes isolated points at all levels.
!     (points sandwiched between land on lateral boundaries:
!           land to N and S, or, land to E and W)
!---------------------------------------------------------------------

   implicit none

   integer :: i,j,k,nlng,nlat,nlev
   integer :: im1,ip1,jm1,jp1,nremoved
   integer :: kmt(nlng,nlat), kmtx(nlng,nlat)
   logical :: land,lande,landw,landn,lands

   print *, " "
   print *, "removing isolated points... "

   kmtx = kmt   ! copy kmt into temp array

   do k = nlev,1,-1

10      nremoved = 0

     do i=1,nlng

       if(i.eq.1) then
         im1 = nlng
       else
         im1 = i-1
       endif
       if(i.eq.nlng) then
         ip1 = 1
       else
         ip1 = i+1
       endif

     do j=2,nlat-1    ! j=1 and j=nlat are assumed to be land

       jm1 = j-1 
       jp1 = j+1 

!         ...note: the following must be modified for tripole grids

       land   = k.gt.kmtx(i,j)
       lande  = k.gt.kmtx(ip1,j)
       landw  = k.gt.kmtx(im1,j)
       landn  = k.gt.kmtx(i,jp1)
       lands  = k.gt.kmtx(i,jm1)


       if (.not.land.and.((lande .and. landw).or.(landn.and.lands))) then
         kmtx(i,j) = kmtx(i,j) - 1   ! remove cell at level k
         nremoved = nremoved + 1
       endif

     enddo  ! j
     enddo  ! i 

     print *, 'points removed this pass: ', k, nremoved
     if (nremoved .ne. 0) go to 10

   enddo    ! k

   kmt = kmtx    ! overwrite kmt 

end subroutine remove_isolated_points

!---------------------------------------------------------------------

subroutine remove_disconnected_bays(iocean,jocean,imt,jmt,KMTX)

!-----------------------------------------------------------------------
!     remove isolated lakes, bays, single points from grid (KMT field)
!     note: this routine assumes land along northern and southern
!     boundaries: KMTX(:,1) = KMTX(:,jmt) = 0
!-----------------------------------------------------------------------

   implicit none

   integer :: i,j,iocean,jocean,imax,imt,jmt
   logical, dimension(imt,jmt) :: MASKC
   integer, dimension(imt,jmt) :: KMTX,KMUX,KINT
   integer, dimension(imt,jmt) :: CFN,CFS,CFE,CFW,I1,I2
   real*4, dimension(imt,jmt) :: CN,CS,CE,CW,C0

   print *, " "
   print *, " removing disconnected bays..."

   do i=1,imt
   do j=1,jmt
     I1(i,j) = i
     I2(i,j) = j
   enddo
   enddo
   KINT = I1 + (I2-1)*imt

   KMUX = min(KMTX,cshift(KMTX,+1,1))
   KMUX = min(KMUX,cshift(KMUX,+1,2))

   CFE = KMUX + cshift(KMUX,-1,2)
   CFW = cshift(CFE,-1,1)
   CFN = KMUX + cshift(KMUX,-1,1)
   CFS = cshift(CFN,-1,2)

   where (CFE.gt.0)
      CFE = 1
   elsewhere
      CFE = imt*jmt + 1
   endwhere

   where (CFW.gt.0)
      CFW = 1
   elsewhere
      CFW = imt*jmt + 1
   endwhere

   where (CFN.gt.0)
      CFN = 1
   elsewhere
      CFN = imt*jmt + 1
   endwhere

   where (CFS.gt.0)
      CFS = 1
   elsewhere
      CFS = imt*jmt + 1
   endwhere

   MASKC = .true.
   where (KMTX.eq.0)
      MASKC = .false.
      KINT = imt*jmt + 1
   endwhere

   imax = (imt+jmt)*3
   do i = 1,imax
      CE = cshift(KINT,+1,1)*(1.0*CFE)
      CW = cshift(KINT,-1,1)*(1.0*CFW)
      CN = cshift(KINT,+1,2)*(1.0*CFN)
      CS = cshift(KINT,-1,2)*(1.0*CFS)
      C0 = KINT*1.0
      where (MASKC)
         KINT = nint(min(CE,CW,CN,CS,C0))
      endwhere
   enddo

!     print *, 'KMTX(iocean,jocean) =', KMTX(iocean,jocean)
!     print *, 'KINT(iocean,jocean) =', KINT(iocean,jocean)

   where (KINT.ne.KINT(iocean,jocean))
      KINT = 0
   elsewhere
      KINT = KMTX
   endwhere

   print *, ' '
   print *, 'surface count =',count(KINT.ne.0)

   KMTX = KINT

end subroutine remove_disconnected_bays

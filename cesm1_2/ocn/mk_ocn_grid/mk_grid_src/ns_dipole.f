
      program ns_dipole

      implicit none

      integer nlng,nlat,nlat_n,nlat_s,ndelta,jcon,n
      integer i,j,jm1,jp1,im1,ip1
      real*8 phi_pole_n,lambda_pole_n,phi_pole_s,lambda_pole_s
      real*8 pi,dx,delta_pole,delta
!!    real*8 deq,dsig,pi,dx,delta_pole,delta
      real*8, allocatable :: ulat_n(:,:),ulng_n(:,:),utan_n(:,:)
      real*8, allocatable :: htn_n(:,:),hte_n(:,:)
      real*8, allocatable :: ulat_s(:,:),ulng_s(:,:),utan_s(:,:)
      real*8, allocatable :: htn_s(:,:),hte_s(:,:)
      real*8, allocatable :: ulat(:,:),ulng(:,:),utan(:,:)
      real*8, allocatable :: htn(:,:),hte(:,:)
      real*8, allocatable :: hus(:,:),huw(:,:)
      double precision deq,dsig
      character*80 filename

      pi = 3.141592653589793d0

      print *, " "
      print *, "nlng:    # zonal grid points"
      print *, "nlat_n:  # meridional grid points in northern hemishere"
      print *, "nlat_s:  # meridional grid points in southern hemishere"
      print *, " "
      print *, "enter nlng, nlat_n, nlat_s"
      read(*,*) nlng,nlat_n,nlat_s
      print *, ">>>>>",nlng,nlat_n,nlat_s

      nlat = nlat_s + nlat_n 

      print *, " "
      print *, "global array dimensions are: ", nlng, nlat

      print *, " "
      print *, " lamda_pole_n : longitude of northern grid pole"
      print *, "   phi_pole_n :  latitude of northern grid pole"
      print *, " lamda_pole_s : longitude of southern grid pole"
      print *, "   phi_pole_s :  latitude of southern grid pole"
200   print *, " "
      print *, "enter lambda_pole_n,phi_pole_n (degrees)"
      read(*,*) lambda_pole_n,phi_pole_n
      print *, ">>>>>",lambda_pole_n,phi_pole_n
      if(phi_pole_n.le.0.d0) then
        print *, " "
        print *, " north grid pole latitude must be > 0"
        go to 200
      endif
      if(abs(lambda_pole_n).gt.180.d0) then
        print *, " "
        print *, 
     &  "enter a lambda_pole_n between -180 (180W) and +180 (180E)"
        go to 200
      endif

300   print *, " "
      print *, "enter lambda_pole_s,phi_pole_s (degrees)"
      read(*,*) lambda_pole_s,phi_pole_s
      print *, ">>>>>",lambda_pole_s,phi_pole_s
      if(phi_pole_s.ge.0.d0) then
        print *, " "
        print *, " south grid pole latitude must be < 0"
        go to 300
      endif
      if(abs(lambda_pole_s).gt.180.d0) then
        print *, " "
        print *, 
     &  "enter a lambda_pole_s between -180 (180W) and +180 (180E)"
        go to 300
      endif

      dx = 360.d0/float(nlng)
      delta_pole = lambda_pole_n - lambda_pole_s
      delta = delta_pole/dx
      ndelta = nint(delta)
      if (delta.ne.float(ndelta)) then
        delta_pole = float(ndelta)*dx
        lambda_pole_s = lambda_pole_n - delta_pole
        print *, " "
        print *, " Longitudes of north and south grid poles"
        print *, " must differ by an integer times the "
        print *, " Equatorial grid spacing. Southern"
        print *, " grid pole latitude has been changed"
        print *, " to the nearest allowed value:"
        print *, "   lambda_pole_s = ",lambda_pole_s
        print *, "   delta_pole    = ",delta_pole
      endif

      lambda_pole_n=lambda_pole_n*pi/180.d0
      phi_pole_n=phi_pole_n*pi/180.d0
      lambda_pole_s=lambda_pole_s*pi/180.d0
      phi_pole_s=phi_pole_s*pi/180.d0
      delta_pole = delta_pole*pi/180.d0

      print *, " "
      print *, "deq : the desired meridional spacing along the Equator"
      print *, "dsig: gaussian half width for meridional profile for"
      print *, "      enhanced meridional spacing along the Equator:"
      print *, "jcon: if not 0, forces (nearly) constant meridional"
      print *, "      grid spacing in the last jcon rows of cells next"
      print *, "      to the grid poles (used to move grid boundaries "
      print *, "      closer to the grid poles without increasing nlat"
      print *, " "
      print *, "enter deq, dsig (degrees), jcon"
      read(*,*) deq,dsig,jcon
      print *, ">>>>>",deq,dsig,jcon

c     ...allocate northern grid

      allocate (ulat_n(nlng,0:nlat_n))
      allocate (ulng_n(nlng,0:nlat_n))
      allocate (utan_n(nlng,0:nlat_n))
      allocate (hte_n(nlng,0:nlat_n))
      allocate (htn_n(nlng,0:nlat_n))

c     ...calculate northern grid

      print *, " "
      print *, " NORTHERN GRID:"

      call northern_dipole(nlng,nlat_n            ! input
     &   ,phi_pole_n,lambda_pole_n,deq,dsig,jcon  ! input
     &   ,ulat_n,ulng_n,htn_n,hte_n,utan_n)       ! output

c     ...allocate southern grid

      allocate (ulat_s(nlng,0:nlat_s))
      allocate (ulng_s(nlng,0:nlat_s))
      allocate (utan_s(nlng,0:nlat_s))
      allocate (hte_s(nlng,0:nlat_s))
      allocate (htn_s(nlng,0:nlat_s))

c     ...calculate southern grid

      print *, " "
      print *, " SOUTHERN GRID:"

      call northern_dipole(nlng,nlat_s             ! input
     &   ,-phi_pole_s,lambda_pole_s,deq,dsig,jcon  ! input
     &   ,ulat_s,ulng_s,htn_s,hte_s,utan_s)        ! output

c     ...shift southern grid arrays to align with northern grid

      ulat_s = cshift(ulat_s,ndelta,1)
      ulng_s = cshift(ulng_s,ndelta,1)
      utan_s = cshift(utan_s,ndelta,1)
      htn_s  = cshift(htn_s ,ndelta,1)
      hte_s  = cshift(hte_s ,ndelta,1)

c     ...allocate global grid

      allocate (ulat(nlng,0:nlat))
      allocate (ulng(nlng,0:nlat))
      allocate (utan(nlng,0:nlat))
      allocate (hte(nlng,0:nlat))
      allocate (htn(nlng,0:nlat))
      allocate (hus(nlng,0:nlat))
      allocate (huw(nlng,0:nlat))

c     ...write southern grid into global arrays

      ulng(:,0:nlat_s) =  ulng_s(:,nlat_s:0:-1)
      ulat(:,0:nlat_s) = -ulat_s(:,nlat_s:0:-1)
      utan(:,0:nlat_s) = -utan_s(:,nlat_s:0:-1)
      htn (:,0:nlat_s) =  htn_s (:,nlat_s:0:-1)
      hte (:,1:nlat_s) =  hte_s (:,nlat_s:1:-1)

c     ...write northern grid into global arrays

      ulng(:,nlat_s:nlat) =  ulng_n(:,0:nlat_n)
      ulat(:,nlat_s:nlat) =  ulat_n(:,0:nlat_n)
      utan(:,nlat_s:nlat) =  utan_n(:,0:nlat_n)
      htn (:,nlat_s:nlat) =  htn_n (:,0:nlat_n)
      hte (:,nlat_s+1:nlat) =  hte_n (:,1:nlat_n)

c     ...calculate hus,huw

      do j = 1,nlat
         jm1 = j-1
         jp1 = j+1
         if (j.eq.1) jm1 = 1
         if (j.eq.nlat) jp1 = nlat
         do i = 1,nlng
            im1 = i-1
            ip1 = i+1
            if (i.eq.1) im1 = nlng
            if (i.eq.nlng) ip1 = 1
            hus(i,j) = 0.25d0*(htn(i  ,j) + htn(i  ,jm1)
     &                       + htn(ip1,j) + htn(ip1,jm1))
            huw(i,j) = 0.25d0*(hte(i,  j) + hte(im1,j  )
     &                       + hte(i,jp1) + hte(im1,jp1))
         enddo
      enddo

c999  continue
c     print *, " "
c     print *, "print column of ulngs, ulats, enter i "
c     print *, "(-1 to continue)"
c     read(*,*)i
c     if(i.eq.-1)go to 1000
c     do j=0,nlat
c      print *,j,ulng(i,j),ulat(i,j),htn(i,j),hte(i,j)
c      print *,j,htn(i,j),hte(i,j),hus(i,j),huw(i,j)
c     enddo 
c     go to 999
c1000 continue

c     ...write direct access grid file grid.plot.da for plotting

      print *, " "
      print *, "enter: "
      print *, "  0 to continue "
      print *, "  1 to output ulat, ulng to direct access file"
      read(*,*) n
      print *, ">>>>>",n

      if(n.eq.1) then

      print *, "enter name of direct access file containing"
      print *, "double-precision fields:"
      print *, "ulat(1:nlng,0:nlat),ulng(1:nlng,0:nlat)"
      read(*,'(a80)') filename
      print *, ">>>>>",filename

      open (53, file=filename, form='unformatted', status='replace',
     & access='direct',action='write',recl=(nlat+1)*nlng*8)
      write (53,rec=1) ulat
      write (53,rec=2) ulng
      close(53)
      print *, " "
      print *, " file written: ", filename

      endif

c     ...write direct access grid file grid.pop.da for POP

      print *, " "
      print *, "enter: "
      print *, "  0 to continue "
      print *, "  1 to direct access grid file for POP"
      read(*,*) n
      print *, ">>>>>",n

      if(n.eq.1) then

        print *, "enter name of direct access POP grid file"
        print *, "containing double-precision fields:"
        print *, " ulat(1:nlng,1:nlat)"
        print *, " ulng(1:nlng,1:nlat)"
        print *, "  htn(1:nlng,1:nlat)"
        print *, "  hte(1:nlng,1:nlat)"
        print *, "  hus(1:nlng,1:nlat)"
        print *, "  huw(1:nlng,1:nlat)"
        print *, " angle(1:nlng,1:nlat)"
        read(*,'(a80)') filename
      print *, ">>>>>",filename

        open (53, file=filename, form='unformatted', status='replace',
     &   access='direct',action='write',recl=nlat*nlng*8)
        write (53,rec=1) ulat(:,1:nlat)
        write (53,rec=2) ulng(:,1:nlat)
        write (53,rec=3) htn(:,1:nlat)
        write (53,rec=4) hte(:,1:nlat)
        write (53,rec=5) hus(:,1:nlat)
        write (53,rec=6) huw(:,1:nlat)
        write (53,rec=7) utan(:,1:nlat)
        close(53)
        print *, " "
        i = nlat*nlng*8
        print *, " file written: ", filename, " Record length = ", i
        i = nlat*nlng*7*8
        print *, " file length = ", i
      endif

c     ...write ascii grid file grid.ascii for VTK

c     print *, " "
c     print *, " enter (1/0) to write a ascii grid file for VTK"
c     read(*,*) i

c     if (i.eq.1) then
c       filename = './grid.ascii'
c       open(53,file=filename,form='formatted',status='unknown')
c       write (53,*)
c    &  'DIPOLE GRID:  ULAT, ULNG, I,    J,     
c    &  HTN,     HTE,     HUS,     HUW,     UTAN'
c       write (53,'(6i8)') nlng,                    ! Num of I-values
c    &                  1 + nlat,                   ! Num of J-values
c    &                  nlng/2,                        ! Crossover I index
c    &                  nlng-1,                 ! (==East I index)
c    &                  nlat/2,                          ! Crossover J index
c    &                  nlat                      ! North J index
c       do j = 0, nlat
c         write (53,*) 'ROW J = ', j
c           ! At a roughly constant latitude J, loop over the longitude
c           ! points and write records LON, LAT, I, J.
c         do i =  1,nlng
c             write (53,'(2f15.7,2i9,x,f15.7
c    &         ,x,f15.7,x,f15.7,x,f15.7,x,f15.7)') 
c    &        -180.0+ulng(i,j)*180.d0/pi, ulat(i,j)*180./pi, i-1, j 
c    &      , htn(i,j)*0.01, hte(i,j)*0.01                
c    &      , hus(i,j)*0.01, huw(i,j)*0.01                
c    &      , utan(i,j)*180./pi
c         end do
c       end do
c       close(53)
c       print *, " "
c       print *, " file written: ", filename
c     endif

      end
 
c*********************************************************************
c*********************************************************************

      subroutine northern_dipole(nlng,nlat    ! input
     &   ,phi_pole,lambda_pole,deq,dsig,jcon  ! input
     &   ,ulat,ulng,htn,hte,utan)             ! output

      implicit none

      integer nlng,nlat,jcon,nlast,i,j,j_match
      double precision phi_pole,lambda_pole,fi,pi,phi_match
      double precision ulat(nlng,0:nlat),ulng(nlng,0:nlat)
      double precision utan(nlng,0:nlat)
      double precision hte(nlng,0:nlat),htn(nlng,0:nlat)
      double precision deq,dsig,xsquare
      double precision xdeg(nlng),ydeg(0:nlat)
      double precision x0(1:nlng),y0(1:nlng),z0(1:nlng)

c---------------------------------------------------------------------
c     set some constants, initialize output arrays
c---------------------------------------------------------------------

      pi = 3.141592653589793d0 
      xsquare = 0.5d0*pi
      nlast = 10000
c     jcon = 0
      j_match = 0
      phi_match = 0.d0

      ulat = 0.d0
      ulng = 0.d0
      utan = 0.d0
      htn = 0.d0
      htn = 0.d0
      hte = 0.d0

c---------------------------------------------------------------------
c     calculation of x-spacing at equator
c---------------------------------------------------------------------

c     ... constant x-spacing assumed

      fi = 2.0d0*pi/float(nlng)
      do i = 1,nlng
        xdeg(i) = i*fi
        ulat(i,j_match) = phi_match
        ulng(i,j_match) = xdeg(i)*cos(phi_match)
      enddo
 
c---------------------------------------------------------------------
c     calculation of the grid
c---------------------------------------------------------------------
 
      call gridgen(nlng,nlat                           ! input
     &            ,phi_pole,lambda_pole,phi_match,xdeg ! input
     &            ,j_match,jcon,nlast,deq,dsig,xsquare ! input
     &            ,ulat,ulng                           ! input, output
     &            ,htn,hte,utan)               ! output 
 
      end

c*********************************************************************
c*********************************************************************

      subroutine gridgen(nlng,nlat                         !input
     &            ,phi_pole,lambda_pole,phi_match,xdeg     ! input
     &            ,j_match,jcon,nlast,deq,dsig,xsquare     ! input
     &            ,ulat,ulng                           ! input, output
     &            ,htn,hte,utan)                   ! output)
 
c---------------------------------------------------------------------
c     this routine generates grid variables for a composite
c     mesh consisting of a polar mesh below phi = phi_match
c     and a distorted mesh with a displace pole for
c     phi > phi_match
c---------------------------------------------------------------------
c     input variables (all angles in radians):
c
c       phi_pole     latitude  of true displaced pole 
c       lambda_pole       longitude of true displaced pole
c       phi_match    latitude  where grids match together 
c       ydeg(0:nlat)      latitude spacing at lambda_pole + xsquare
c---------------------------------------------------------------------
c     output variables:
c
c       ulat,ulng,utan,hte,htn
c---------------------------------------------------------------------
 
      implicit none

      integer nlng,nlat,j_match,jcon,nlast
      double precision p25
     &  , del, dmin, dmax, alpha0, phi0, fy, dxavg, phiold
     &  , phi_match, dy, yy, aspect, dxmin, dymin
     &  , deq, dsig, xsquare
      double precision radius,delta,c0,d0,c00,d00,yl,yu,gg
      double precision lambda_pole,phi_pole,smean,srms
      double precision x0(nlng),y0(nlng),z0(nlng)
      double precision x00(nlng),y00(nlng),z00(nlng)
      double precision ulat(nlng,0:nlat),ulng(nlng,0:nlat)
      double precision utan(nlng,0:nlat)
      double precision hte(nlng,0:nlat),htn(nlng,0:nlat)
      double precision xdeg(nlng),ydeg(0:nlat)
      double precision htep(nlng),htnp(nlng)
      integer i,j,im1,ip1,jm1,jp1,n,ii,jj,imax,jmax,imin,jmin,nnn,iii
 
      parameter (p25 = 0.25d0)
c     parameter (radius = 6370.0d5)
      parameter (radius = 6371.22d5)
      double precision pi
      parameter (pi = 3.141592653589793d0)

c*********************************************************************

      dxavg = 360.d0/float(nlng)

      dmin = pi
      do i = 1,nlng/2
        delta = abs(xdeg(i)-xsquare)
        if (delta.lt.dmin) then
          dmin = delta
          ii = i
        endif
      enddo

      j = j_match
      delta = cos(phi_match)
      do i = 1,nlng
        ulng(i,j) = xdeg(i)
        ulat(i,j) = phi_match
        if (i.gt.1) then
          htn(i,j) = delta*(ulng(i,j)-ulng(i-1,j))*radius
        else
          htn(i,j) = delta*ulng(i,j)*radius
        endif
      enddo
      ydeg(j) = ulat(nlng/2,j)

      do j = j_match-1,0,-1
	print *,'in AAAA'
        fy = dsig*pi/180.d0
        fy = 1.0d0 - (1.0d0-deq/dxavg)*exp(-0.5d0*(ulat(ii,j+1)/fy)**2)
c       fy = 1.0d0
        jj = j
        if(j.lt.jcon)jj=jcon
        del = (xdeg(ii) - xdeg(ii-1))*cos(ulat(ii,jj+1))*fy
        do i = 1,nlng
          ulng(i,j) = xdeg(i)
          ulat(i,j) = ulat(i,j+1) - del
          hte(i,j+1) = abs(ulat(i,j) - ulat(i,j+1))*radius
          delta = cos(ulat(i,j))
          if (i.gt.1) then
            htn(i,j) = delta*(ulng(i,j)-ulng(i-1,j))*radius
          else
            htn(i,j) = delta*ulng(i,j)*radius
          endif
        enddo
        ydeg(j) = ulat(nlng/2,j)
      enddo

      do j = 0,j_match-1
	print *,'in BBB'
        fy = dsig*pi/180.d0
        fy = 1.0d0 - (1.0d0-deq/dxavg)*exp(-0.5d0*(ulat(ii,j+1)/fy)**2)
c       fy = 1.0d0
c       print *,j,ydeg(j)*180.d0/pi,fy
      enddo

c   ... initialize x0,y0,z0,c0,d0

      c0 = 0.
      d0 = sin(ulat(1,j_match))

      do i = 1,nlng
        x0(i) = cos(ulat(i,j_match))*cos(ulng(i,j_match))
        y0(i) = cos(ulat(i,j_match))*sin(ulng(i,j_match))
        z0(i) = sin(ulat(i,j_match))
      enddo

c   ... displaced-pole grid above j_match
        
      yy = 0.0
      phiold = phi_match + yy
      dy = 0.01d0*pi/180.  ! 100th of a degree 

      do j = j_match,nlat-1

        do i = 1,nlng
          hte(i,j+1) = 0.0
          htn(i,j+1) = 0.0
        enddo

        jj = j
        if(j.gt.nlat-jcon)jj=nlat-jcon

        fy = dsig*pi/180.d0
        fy = 1.0d0 - (1.0d0-deq/dxavg)*exp(-0.5*(phiold/fy)**2)

        do n = 1,nlast

c         ...set alpha0, phi0

          yy = yy + dy
          delta = yy/(pi-phi_pole-phi_match)
          delta = pi*delta*0.5d0
          delta = sin(delta)
          alpha0 = pi*0.5d0 - (pi*0.5d0-phi_pole)*delta
          phi0 = phi_match + yy

          c00 = c0
          d00 = d0
          do i = 1,nlng
            x00(i) = x0(i)
            y00(i) = y0(i)
            z00(i) = z0(i)
          enddo

          call jcircle(j,nlng,nlat,alpha0,phi0    ! input
     &                ,x0,y0,z0,c0,d0             ! input, output
     &                ,htep,htnp)                 ! output

c         ...increment sum for hte

          if (hte(ii,j+1)+htep(ii).ge.fy*htn(ii,jj)) go to 1000

          do i = 1,nlng
            hte(i,j+1) = hte(i,j+1) + htep(i)
          enddo

        enddo

1000    continue

        yu = yy
        yl = yy - dy

        do n = 1,nlast

          yy = 0.5d0*(yu + yl)
          delta = yy/(pi-phi_pole-phi_match)
          delta = pi*delta*0.5d0
          delta = sin(delta)
          alpha0 = pi*0.5d0 - (pi*0.5d0-phi_pole)*delta
          phi0 = phi_match + yy
          c0 = c00
          d0 = d00
          do i = 1,nlng
            x0(i) = x00(i)
            y0(i) = y00(i)
            z0(i) = z00(i)
          enddo

          call jcircle(j,nlng,nlat,alpha0,phi0    ! input
     &                ,x0,y0,z0,c0,d0             ! input, output
     &                ,htep,htnp)                 ! output

          gg = fy*htn(ii,jj) - (hte(ii,j+1)+htep(ii))

c         if (abs(gg).le.0.1d0) go to 2000
          if (abs(gg).le.1.0d0) go to 2000
          if(gg.gt.0.0d0) then
            yl = yy
          else
            yu = yy
          endif

        enddo

        print *, "no converge: ", fy*htn(ii,jj),hte(ii,j+1)+htep(ii)
2000    continue

        do i = 1,nlng
          hte(i,j+1) = hte(i,j+1) + htep(i)
        enddo
      
        phiold = phi0
        ydeg(j+1) = phi0
c       print *, j+1,phi0*180.d0/pi,n,fy

c       ...find ulat, ulng, and angle utan
c          also copy htn to 2d array

        do i = 1,nlng

          htn(i,j+1) =  htnp(i)

          ulat(i,j+1) = asin(z0(i))
          if (x0(i).ne.0.0d0.or.y0(i).ne.0.0d0) 
     &      ulng(i,j+1) = atan2(y0(i),x0(i))

          delta = sqrt((y0(i)**2*(1.d0+c0**2)+(c0*z0(i)-x0(i))**2)
     &          *(1.d0-z0(i)**2))
          delta = (z0(i)*d0-1.d0)/delta
          if (abs(delta).gt.1.d0) delta = sign(1.d0,delta)
          utan(i,j+1) = sign(acos(-delta),c0*y0(i))

        enddo

c     ...symmetrize right and left hemispheres

        do i = 1,nlng/2

          iii = nlng - i
          htn(iii+1,j+1) = htn(i,j+1)
          hte(iii,j+1)   = hte(i,j+1)
          ulng(iii,j+1)  = 2.d0*pi - ulng(i,j+1)
          ulat(iii,j+1)  = ulat(i,j+1)
          utan(iii,j+1)  = -utan(i,j+1)

        enddo

      enddo

C     print *, " output hte(1,:) to file? (1-yes, 0-no)"
C     read(*,*) n
C     if(n.eq.1)then
C       open(77,file='hte.dat',form='formatted',status='unknown')
C       do j=1,nlat
C         write(77,*) j,hte(1,j)
C       enddo
C       close(77)
C     endif

C     print *, " "
C     do j = 0,nlat-1
C       if(j.eq.0)then
C         print *, j,ydeg(j)
C       else
C         print *, j,ydeg(j)*180./pi,(ulat(ii,j)-ulat(ii,j-1))*180./pi
C       endif
C     enddo

c   ... print min,center,max latitudes of last circle

      print *, " "
      print *, " latitudes of highest circle"
      print *, "  low point:",
     &  (180.d0/pi)*acos((c0*d0 + sqrt(1.d0+c0**2-d0**2))/(1.d0+c0**2))
      print *, "     center:",
     &  (180.d0/pi)*acos((c0*d0)/(1.d0+c0**2))
      print *, " high point:",
     &  (180.d0/pi)*acos((c0*d0 - sqrt(1.d0+c0**2-d0**2))/(1.d0+c0**2))
 
c   ... adjust longitudes 
      do j=0,nlat
         do i=1,nlng
	    ulng(i,j) = mod(ulng(i,j)+lambda_pole+6*pi,2*pi)
         enddo
      enddo
 
c   ...print min cell widths

      dxmin = 1.0d+10
      dymin = 1.0d+10
      do i = 1,nlng
      do j = 1,nlat
        delta = htn(i,j)
        if((delta.lt.dxmin) .and. (delta.ne.0.d0)) then
          imax = i
          jmax = j
          dxmin = delta
        endif
        delta = hte(i,j)
        if((delta.lt.dymin) .and. (delta.ne.0d0)) then
          imin = i
          jmin = j
          dymin = delta
        endif
      enddo
      enddo
      print *, " "
      print *, " min dx :"
      print *, imax,jmax,htn(imax,jmax),htn(nlng/2,nlat)
     &                                 ,htn(1,1)
      print *, " long,lat: "
      print *, ulng(imax,jmax)*180.d0/pi, ulat(imax,jmax)*180.d0/pi
      print *, " "
      print *, " min dy:"
      print *, imin,jmin,hte(imin,jmin)
      print *, " lat,long: "
      print *, ulng(imin,jmin)*180.d0/pi, ulat(imin,jmin)*180.d0/pi
      print *, " "

c   ...print max, min aspect ratios

      dmin = 1.0d0
      dmax = 1.0d0
      smean = 0.0d0
      srms = 0.0d0
      do i = 1,nlng
      do j = 1,nlat
        aspect = htn(i,j)/hte(i,j)
        delta = abs(aspect-dmax)
        if(delta.gt.dmax)then
          imax = i
          jmax = j
          dmax = delta
        endif
        delta = abs(aspect-dmin)
        if(delta.lt.dmin)then
          imin = i
          jmin = j
          dmin = delta
        endif
c       if(aspect.lt.1.0d0) aspect = 1.0d0/aspect
        smean = smean + aspect
        srms = srms + aspect**2
      enddo
      enddo
      print *, " "
      print *, " max aspect ratio:"
      print *, imax,jmax,htn(imax,jmax)/hte(imax,jmax) 
      print *, " long,lat: "
      print *, ulng(imax,jmax)*180.d0/pi, ulat(imax,jmax)*180.d0/pi
      print *, " "
      print *, " min aspect ratio:"
      print *, imin,jmin,htn(imin,jmin)/hte(imin,jmin) 
      print *, " lat,long: "
      print *, ulng(imin,jmin)*180.d0/pi, ulat(imin,jmin)*180.d0/pi
      print *, " "
      print *, " mean aspect ratio:"
      print *, smean/float(nlng*nlat)
      print *, " rms  aspect ratio:"
      print *, sqrt(srms/float(nlng*nlat))
      print *, " "

c   ...output xdeg,ydeg to files

c     open(17,file='xdeg.out',status='unknown')
c     do i=1,nlng
c       write(17,*) xdeg(i)*180.d0/pi, fx(i)*180.d0/pi
c     enddo
c     close(17)

c     open(18,file='ydeg.out',status='unknown')
c     do j=1,nlat
c       write(18,*) ydeg(j)*180.d0/pi, (ydeg(j)-ydeg(j-1))*180.d0/pi
c     enddo
c     close(18)
 
c   ...output circumpherences of j-circles and i-circles

c     open(11,file='jcircle.out',status='unknown')
c     do j = 1,nlat
c       hjcircle(j) = 0.0d0
c       do i = 1,nlng
c         hjcircle(j) = hjcircle(j) + htn(i,j)
c       enddo
c       write(11,*) j, hjcircle(j)*1.d-5  ! (km)
c     enddo
c     close(11)

c     open(12,file='icircle.out',status='unknown')
c     del = 0.0d0
c     do i = 1,nlng
c       hicircle(i) = 0.0d0
c       do j = j_match,nlat
c         hicircle(i) = hicircle(i) + hte(i,j)
c       enddo
c       del = del + hicircle(i)
c     enddo
c     do i = 1,nlng
c       write(12,*) i, (hicircle(i)/del)*360.d0  ! (degrees)
c     enddo
c     close(12)

      end
 
c*********************************************************************

      subroutine jcircle(j,nlng,nlat,alpha0,phi0    ! input
     &                  ,x0,y0,z0,c0,d0             ! input, output
     &                  ,htep,htnp)                 ! output

c     ...

      implicit none
      integer nlng,nlat
      integer i,j,im1,n
      double precision a,b,c0,d0,c,d,c1,d1,c2,atmp,btmp,ctmp,ar,ba,bb
      double precision alph,alpha,phi,xq,zq,sq,radius
      double precision alpha0,phi0,rr
      double precision x,y,z,x1,y1,z1,r1tmp,r2tmp,r3tmp
      double precision p,psi,eps,amu,anu,q,r,delta,rnorm,r1norm 
      double precision x0(1:nlng),y0(1:nlng),z0(1:nlng)
      double precision xi(nlng),yi(nlng),zi(nlng)
      double precision htep(nlng),htnp(nlng)
      data eps,radius/1.d-8,6371.22d5/
      double precision pi
      parameter (pi = 3.141592653589793d0)
      
c---------------------------------------------------------------------

      alph = (sqrt(1.d0+c0**2-d0**2) + c0*d0)/(1.d0+c0**2)
      alph = acos(alph)
c     print *, j,' alph =',alph*180.d0/pi,alpha0*180.d0/pi
      alph = max(alph,alpha0)
c     alph = alpha0
      d1 = 1./sin(alph)
      c1 = (d1 - sin(alph))/cos(alph)
      if (abs(alph-pi/2).lt.eps) then
        psi = 0.0d0
        xq = 0.0d0
        zq = 0.0d0
      else
        xq = (d1 - d0)/(c1 - c0)
        zq = d0 - c0*xq
        psi = atan(zq/xq)
      endif
      phi = phi0 + psi         ! phi(j+1)   in system j
      alpha = alph - psi            ! alpha(j) in system j
      c = sin(phi)/(1.d0/cos(alpha) + cos(phi)) ! c(j+1) " "
      d = c/cos(alpha)                        ! d(j+1) " "
      anu = atan(c)
      if(abs(cos(alpha)).gt.eps) then
         r3tmp = sin(anu)/cos(alpha)
      else
         r3tmp = 0.d0
      endif
      r2tmp = sqrt(1.d0 - r3tmp*r3tmp)

      do i = 1,nlng

c       ... rotate to find x(j),y(j),z(j) in system j

        x1 =  x0(i)*cos(psi) + z0(i)*sin(psi)
        z1 = -x0(i)*sin(psi) + z0(i)*cos(psi)
        y1 =  y0(i)

c       ...solving px^2 + 2qx + r = 0

        if (abs(y1).gt.abs(x1-cos(alpha))) then
          ar = (cos(alpha)-x1)/y1
          ba = cos(alpha)
          bb = sign(1.0d0,ar)*cos(alpha)
          p = 1.d0 + (1.d0+c**2)*ar**2
          q = -(ba + c*d*ar**2)
          r = ba**2 + (d**2 - 1.d0)*ar**2
          sq = q*q - r*p
          delta = sqrt(max(sq,0.d0))
          if (i.eq.nlng .or. bb*(nlng/2-i).lt.0.d0) then
            x = (-q+delta)/p
          else
            x = (-q-delta)/p
          endif
          z = d - c*x
          y = sqrt(max(1.d0 - z**2 - x**2,0.d0))
          if(i.gt.nlng/2) y = -y
          rr = sqrt(x**2+y**2+z**2)
          x=x/rr
          y=y/rr
          z=z/rr
          if(abs(ar).gt.eps) then
            amu = atan(1.d0/ar)
          else
            amu = pi/2.d0
          endif
        else
          a = y1/(cos(alpha)-x1)
          b = a*cos(alpha)
          p = 1.d0 + a**2 + c**2
          q = -(c*d+a*b)
          r = d**2 + b**2 - 1.d0
          sq = q*q - r*p
          delta = sqrt(max(sq,0.d0))
          if (i.eq.nlng .or. b*(nlng/2-i).lt.0.d0) then
            x = (-q+delta)/p
          else
            x = (-q-delta)/p
          endif
          y = b - a*x
          z = d - c*x
          r = sqrt(x**2+y**2+z**2)
          x=x/r
          y=y/r
          z=z/r
          amu = atan(a)
        endif

c       ...rotation of angle -psi around the y-axis

        y0(i)=y
        x0(i)=cos(-psi)*x+sin(-psi)*z
        z0(i)=-sin(-psi)*x+cos(-psi)*z
        r = sqrt(x0(i)**2+y0(i)**2+z0(i)**2)
        x0(i) = x0(i)/r
        y0(i) = y0(i)/r
        z0(i) = z0(i)/r

36      format(2i5,x,4e11.4)
 
c       ...rotate to find c0,d0

        c1 =  c*cos(-psi) + sin(-psi)
        c2 = -c*sin(-psi) + cos(-psi)
        c0 = c1/c2
        d0 = d/c2
 
c       ... calculate hte

        delta = sin(amu)*cos(alpha)
        r1tmp = sqrt(1. - delta*delta)
        p = delta*sin(amu)
        q = delta*cos(amu)
        rnorm = sqrt((x-p)**2+(y-q)**2+z**2)
        r1norm = sqrt((x1-p)**2+(y1-q)**2+z1**2)
        atmp = (x-p)*(x1-p)
        btmp = (y-q)*(y1-q)
        ctmp = z*z1
        delta = (atmp + btmp + ctmp)/(rnorm*r1norm)
c       if (delta.gt.1.0d0) then
c         print *, i,j,delta,abs(x1-cos(alpha))
        delta = min(delta,1.0d0)
c       endif
        delta = acos(delta)
        htep(i) = r1tmp*delta*radius

c       ... calculate temporaries needed to calculate htn

        atmp = r3tmp*sin(anu)
        btmp = r3tmp*cos(anu)
        xi(i) = x - atmp
        yi(i) = y
        zi(i) = z - btmp

      enddo

c     ... calculate htn

      do i = 1,nlng

        im1 = i-1
        if (i.eq.1) im1 = nlng
        atmp = xi(i)*xi(im1)
        btmp = yi(i)*yi(im1)
        ctmp = zi(i)*zi(im1)
        delta = acos((atmp + btmp + ctmp)/(r2tmp**2))
        htnp(i) = r2tmp*delta*radius

      enddo

      end

c*********************************************************************

 
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c     end file dipole.f
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

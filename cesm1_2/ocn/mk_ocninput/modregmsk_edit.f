      program region

      implicit none
      integer :: i, j, k, recl_factor
      integer, parameter :: imt=NX, jmt=NY
      integer*4, dimension(imt,jmt) :: kmt, mask 
      double precision ulat(imt,jmt)

      integer temp(imt,jmt)

      open (1, file='kmt.da', 
     &      access='direct',
     &      form='unformatted', recl=jmt*imt*4,
     &      status='old')

      open (2, file='grid.pop.da', 
     &      access='direct',
     &      form='unformatted', recl=jmt*imt*8,
     &      status='old')

      open (3, file='region.ieeei4', 
     &      access='direct',
     &      form='unformatted', recl=jmt*imt*4,
     &      status='unknown',
     &      action='write')

c read in topo levels:
      read  (1, rec=1) kmt
c read in grid:
      read  (2, rec=1) ulat
     
      mask = kmt
      do i = 1,imt 
      do j = 1,jmt 
ccthis should work since ulat=0 at eq and has same j as kmt box just north.
       if(ulat(i,j).ge.0.d0.and.kmt(i,j).gt.0) mask(i,j)=1 
       if(ulat(i,j).lt.0.d0.and.kmt(i,j).gt.0) mask(i,j)=2 
c this should work since ulat=0 at eq and has same j as kmt box just south.

c north panthalassa
c        if(ulat(i,j).gt.0.d0.and.kmt(i,j).gt.0) mask(i,j)=1 
c south panthalassa
c        if(ulat(i,j).le.0.d0.and.kmt(i,j).gt.0) mask(i,j)=2 

c tethys
c        if(i.ge.0.and.i.le.57.and.j.ge.50.and.j.le.280.and.
c     &     kmt(i,j).gt.0) mask(i,j)=3
c        if(i.ge.280.and.i.le.320.and.j.ge.75.and.j.le.254.and.
c     &     kmt(i,j).gt.0) mask(i,j)=3

c return some cells to north panthalassa
c        if(i.ge.45.and.i.le.60.and.j.ge.260.and.j.le.284.and.
c     &     kmt(i,j).gt.0) mask(i,j)=1

      enddo
      enddo

      write (unit=3, rec=1)mask
      end

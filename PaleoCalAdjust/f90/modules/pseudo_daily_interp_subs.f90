module pseudo_daily_interp_subs
    
implicit none

    logical             :: debug=.false.
    integer             :: debug_unit = 10
    integer             :: out_unit = 6

contains

subroutine hdaily(nm,nd,xm,monlen,no_negatives,xdh)

    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nd,nm
    real(8), intent(in)     :: xm(nm)
    integer(4), intent(in)  :: monlen(nm)
    logical, intent(in)     :: no_negatives
    real(8), intent(out)    :: xdh(nd)

    real(8)             :: a(0:nh),b(0:nh)
    integer(4)          :: m
    
    if (debug) then
        do m=1,nm
            write (debug_unit,*) m,xm(m),monlen(m)
        end do
    end if

    ! interpolate daily values
    call harmonic_coeffs(nm,xm,a,b)
    call xdhat(nm,nd,monlen,a,b,xdh)
    if (no_negatives) call dzero(nm,nd,monlen,xm,xdh)

end subroutine hdaily

subroutine harmonic_coeffs(nm,y,a,b)
! Calculates a's and b's of an "adjusted" harmonic fit to monthly values of a variable,
! which preserves the monthly (and annual) mean values by interpolated daily values
! adapted from Epstein (1991, On obtaining daily climatological values from monthly means,
! J. Climate 4:365-368).

    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nm
    real(8), intent(in)     :: y(nm)
    real(8), intent(out)    :: a(0:nh),b(0:nh)
    
    real(8)                 :: pi
    real(8)                 :: asum,bsum,c0,c1,c2,c3,c4
    integer(4)              :: j,t
    
    a=0.0d0; b=0.0d0
    pi=4.0d0*datan(1.0d0)
    
    ! a0
    do t=1,nm
        a(0)=a(0)+y(t)
    end do
    a(0)=a(0)/dble(nm) 
    
    ! a's and b's
    do j=1,nh-1
        if (debug) write (debug_unit,'(a)') " "
        asum=0.0d0; bsum=0.0d0
        c1=pi*(dble(j)/dble(nm))
        do t=1,nm
            c0=dble(t)/dble(nm)
            c2=(2.0d0*pi*dble(j)*dble(t))/dble(nm) 
            asum=asum+y(t)*dcos(c2)/(dble(nh))
            bsum=bsum+y(t)*dsin(c2)/(dble(nh))
            if (debug) write (debug_unit,'("j,t,c0,c2,asum,bsum: ",2i3,5f10.4)') j,t,y(t),c0,c2,asum,bsum
        end do
        a(j)=(c1/dsin(c1))*asum
        b(j)=(c1/dsin(c1))*bsum
        if (debug) write (debug_unit,'("j,c1,a,b:",i3,3f10.4)') j,c1,a(j),b(j)
    end do
    
    if (debug) write (debug_unit,'(a)') " "
    asum=0.0d0
    do t=1,nm
        c3=cos(pi*dble(t))/dble(nm) 
        asum=asum+(y(t)*c3)
        if (debug) write (debug_unit,'("t,y,c3,asum: ",i3,5f10.4)') t,y(t),c3,asum,cos(pi*dble(t)),pi*dble(t)
    end do
    c4=((pi/2.0d0)/sin(pi/2.0d0))
    a(nh)=c4*asum !((pi/2.0)/sin(pi/2.0))*asum
    b(nh)=0.0d0
    if (debug) write (debug_unit,'(a)') " "
    if (debug) write (debug_unit,'("c4,asum: ",2f10.4)') c4,asum
    
    if (debug) then
        write (debug_unit,'(a)') " "
        do j=0,nh
                write (debug_unit,'("j,a,b", i3,2f10.4)') j,a(j),b(j)
        end do
        write (debug_unit,'(a)') " "
    end if 
    
end subroutine harmonic_coeffs

subroutine xdhat(nm,nd,monlen,a,b,yhat)
! Calculates/interpolates pseudo-daily values of variable using the a's and b's from harmonic_coeffs()
! adapted from Epstein (1991, On obtaining daily climatological values from monthly means,
! J. Climate 4:365-368).

    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: a(0:nh),b(0:nh)
    real(8), intent(out)    :: yhat(nd)
    
    integer(4)              :: i,j,m,ii
    real(8)                 :: t,pi
    real(8)                 :: c2
    
    pi=4.0d0*datan(1.0d0)
    yhat=0.0
    ii=0
    do i=1,nm
        do m=1,monlen(i)
            ii=ii+1
            t=(dble(i)-0.5d0)+(dble(m)-0.5d0)/dble(monlen(i))
            do j=0,nh
                c2=((2.0d0*pi*dble(j)*t)/dble(nm))
                yhat(ii)=yhat(ii)+a(j)*dcos(c2)+b(j)*dsin(c2)
            end do
            if (debug) write (debug_unit,'("i,m,ii,t,c2,yhat: ",3i4,3f10.4)') i,m,ii,t,c2,yhat(ii)
        end do
    end do
    
end subroutine xdhat

subroutine dayinterp(nm,nd,monlen,zm,zd)
! Interpolate pseudo-daily values of monthly data.  Not mean-preserving.

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: zm(nm)
    real(8), intent(out)    :: zd(nd)
    
    integer(4)              :: nm1
    integer(4)              :: midmon(nm),midmon2(0:nm+1)
    real(8)                 :: zm2(0:nm+1)
    integer                 :: i,m
    
    call midmonth_int(nm,monlen,midmon)
    
    ! pad data at beginning (m=0) and end (m=13)
    nm1=nm+1
    zm2(1:nm)=zm
    zm2(0)=zm(nm)
    zm2(nm1)=zm(1)
    midmon2(1:nm)=midmon
    midmon2(0)=1-(nd-midmon(nm))-2
    midmon2(nm1)=nd+midmon(1)   
    
    ! linear pseudo-daily interpolation
    do i=1,nd
        ! find month day i lies in
        do m=1,nm+1
            if (i.gt.midmon2(m-1) .and. i.le.midmon2(m)) exit
        end do   
        zd(i)=(dble(i-midmon2(m-1))/dble(midmon2(m)-midmon2(m-1)))*(zm2(m)-zm2(m-1))+zm2(m-1)       
    end do
    
end subroutine dayinterp

subroutine dayspread(nm,nd,monlen,zm,zd)
! block fill daily values from monthly means

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    real(4), intent(in)     :: monlen(nm)
    real(8), intent(in)     :: zm(nm)
    real(8), intent(out)    :: zd(nd)
    
    integer                 :: i,j,m
    integer(4)              :: imonlen(nm)

    imonlen = int(monlen)
    
    i=0
    do m=1,nm
        do j=1,imonlen(m)
            i=i+1
            zd(i)=zm(m)
        end do
    end do
    
end subroutine dayspread

subroutine midmonth_int(nm,monlen,midmon)
! gets mid-month day number 

    implicit none
    
    integer(4), intent(in)  :: nm
    integer(4), intent(in)  :: monlen(nm)
    integer(4), intent(out) :: midmon(nm)
    
    integer(4)              :: m,endday(nm)

    ! midmonth day numbers
    m=1
    midmon(m)=ceiling(dble(monlen(m))/2.0d0)
    endday(m)=monlen(m)
    do m=2,nm
        midmon(m)=ceiling(dble(monlen(m))/2.0d0)+endday(m-1)
        endday(m)=endday(m-1)+monlen(m)
    end do
    
end subroutine midmonth_int

subroutine midmonth_real(nm,veq_mon,veq_midmon_day,rmonlen,rmidmon,rmonbeg,rmonend)
! gets mid-month day number 

    implicit none
    
    integer(4), intent(in)  :: nm
    integer(4), intent(in)  :: veq_mon
    real(8), intent(in)     :: veq_midmon_day
    real(8), intent(in)     :: rmonlen(nm)
    real(8), intent(out)    :: rmidmon(nm),rmonbeg(nm),rmonend(nm)
    
    integer(4)              :: m, debug_unit = 10
    
    logical                 :: debug_write = .false.
    
    ! midmonth target values 
    ! if first month is Jan: 0ka equinox:  31.0+28.0+21.5 = 80.5; 0ka March midmonth day:  31.0+28.0+15.5 = 74.5
    ! if first month is Dec: 0ka equinox:  31.0+31.0+28.0+21.5 = 111.5; 0ka March midmonth day:  31.0+31.0+28.0+15.5 = 105.5
    
    rmidmon(veq_mon)=veq_midmon_day ! fixed March midmonth relative to noleap VE
    rmonbeg(veq_mon)=rmidmon(veq_mon)-(rmonlen(veq_mon)/2.0d0)
    rmonend(veq_mon)=rmidmon(veq_mon)+(rmonlen(veq_mon)/2.0d0)
    do m=veq_mon-1,1,-1
        rmidmon(m)=rmidmon(m+1)-(rmonlen(m+1)/2.0d0)-(rmonlen(m)/2.0d0) 
        rmonbeg(m)=rmidmon(m)-(rmonlen(m)/2.0d0)
        rmonend(m)=rmidmon(m)+(rmonlen(m)/2.0d0)
    end do
    do m=veq_mon+1,nm
        rmidmon(m)=rmidmon(m-1)+(rmonlen(m-1)/2.0d0)+(rmonlen(m)/2.0d0) 
        rmonbeg(m)=rmidmon(m)-(rmonlen(m)/2.0d0)
        rmonend(m)=rmidmon(m)+(rmonlen(m)/2.0d0)
    end do
    
    if (debug_write) then
        write (debug_unit,'("veq_mon, veq_midmon_day: ",i3,f11.6)') veq_mon,veq_midmon_day
        write (debug_unit,'("rmonlen   ",12f11.6)') rmonlen
        write (debug_unit,'("rmidmon   ",12f11.6)') rmidmon
        write (debug_unit,'("rmonbeg   ",12f11.6)') rmonbeg
        write (debug_unit,'("rmonend   ",12f11.6)') rmonend
    end if
                
end subroutine midmonth_real


subroutine dzero(nm,nd,monlen,xm,xd0)
! enforces 0.0 values of interpolated daily data when the monthly mean is 0.0

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: xm(nm)
    real(8), intent(inout)  :: xd0(nd)
    
    real(8)                 :: xdm(nm),diff,totaldiff
    integer(4)              :: i,m,j,nonzero(nm),l,maxiter=30
    integer(4)              :: imonlen(nm)
    
    imonlen = int(monlen)
    
    ! zero all negative daily values
    do i=1,nd
        if (xd0(i).le.0.0) xd0(i)=0.0
    end do
    
    do l=1,maxiter  
      ! zero daily values in months where xm=0.0
      i=0
      do m=1,nm
          do j=1,imonlen(m)
              i=i+1
              if (xm(m).eq.0.0) xd0(i)=0.0          
          end do
      end do
  
      i=0
      xdm=0.0d0; nonzero=0; totaldiff=0.0d0
      do m=1,nm
          do j=1,imonlen(m)
              i=i+1
              xdm(m)=xdm(m)+xd0(i)
              if (xdm(m).gt.0.0d0) nonzero(m)=nonzero(m)+1
          end do
          xdm(m)=xdm(m)/dble(monlen(m))
          diff=dabs((xm(m)-xdm(m)))
          totaldiff=totaldiff+diff
      end do
      if (totaldiff.le.0.0001) exit
    
      i=0
      do m=1,nm
          if (nonzero(m).ne.0) then
              do j=1,imonlen(m)
                  i=i+1
                  if (xd0(i).gt.0d0) then
                      xd0(i)=xd0(i)+(xm(m)-xdm(m)) !/nonzero(m)
                  end if
              end do
          end if
      end do    
    end do    
  
  ! zero daily values in months where xm=0.0
    i=0
    do m=1,nm
        do j=1,imonlen(m)
            i=i+1
            if (xm(m).eq.0.0) xd0(i)=0.0          
        end do
    end do
    
    ! zero all negative daily values
    do i=1,nd
        if (xd0(i).le.0.0) xd0(i)=0.0
    end do
    
end subroutine dzero

subroutine monmean(nm,nd,monlen,xd,xm)
! gets monthly means of interpolated daily data

    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: xd(nd)
    real(8), intent(out)    :: xm(nm)
    integer                 :: i,j,m
    
    xm=0.0d0
    i=0
    do m=1,nm
        do j=1,monlen(m)
            i=i+1
            xm(m)=xm(m)+xd(i)
        end do
        xm(m)=xm(m)/dble(monlen(m))
    end do

end subroutine monmean

subroutine ann_wmean(n,x,w,xm)
! gets a weighted (by month length) annual mean value from monthly data

    implicit none
    
    integer(4), intent(in)  :: n
    real(8), intent(in)     :: x(n),w(n)
    real(8), intent(out)    :: xm
    real(8)                 :: wsum
    integer                 :: i
    
    xm=0.0d0; wsum=0.0d0
    do i=1,n
        xm=xm+x(i)*w(i)
        wsum=wsum+w(i)
    end do

    xm=xm/wsum
    
end subroutine ann_wmean
    
subroutine ann_mean(n,x,xm)
! gets an unweighted annual mean value from monthly or daily data

    implicit none
    
    integer(4), intent(in)  :: n
    real(8), intent(in)     :: x(n)
    real(8), intent(out)    :: xm
    integer                 :: i
    
    xm=0.0d0
    do i=1,n
        xm=xm+x(i)
    end do

    xm=xm/dble(n)
    
end subroutine ann_mean   

end module pseudo_daily_interp_subs
module calendar_effects_subs
! This module contains the two main subroutines used to adjust monthly data for the paleo calendar effect:
! mon_to_day_ts() interpolates input monthly data to pseudo-daily values using the Epstein (1991) mean-
!   preserving "harmonic" interplation method, implemented in the module pseudo_daily_interp()
! day_to_mon() aggregates daily or pseudo daily values to monthly totals using real-number definitions of 
!   month lengths.
    
implicit none

    integer             :: debug_unit=10
    integer             :: console_unit=6
    
    logical             :: debug_write_cal_effects=.false.

contains

subroutine mon_to_day_ts(nt,imonlen,xm_in,xfill,no_negatives,smooth,restore,ndtot,nw,nsw,xd_out)
! Daily interpolation of a monthly time series
! Interpolation is done one year at a time, and so there can be small discontinuities between years.
! This version makes one pass over the input time series, and optionally smooths and restores the long-term mean.

    use pseudo_daily_interp_subs

    implicit none
    
    integer(4), parameter   :: nm=12, ndl=366
    integer(4), intent(in)  :: nt, ndtot, nw, nsw
    real(8), intent(in)     :: xm_in(nt), xfill
    integer(4), intent(in)  :: imonlen(nt)
    logical, intent(in)     :: no_negatives,smooth,restore
    real(8), intent(out)    :: xd_out(ndtot)
    
    real(8)                 :: xm(nm),xdh(ndl),xd_temp(ndtot)
    real(8)                 :: xm_ltm,xd_ltm,monlentot,ltmdiff
    real(8)                 :: pi,x,wgt(nw),wsum
    integer(4)              :: iml(nm),nd
    
    integer(4)              :: nyrs
    integer(4)              :: i,j,jj,jjj,js,m,mm,n,nn,nzero,nnonfill,nfill
    
    nyrs=nt/nm    
    
    if (debug_write_cal_effects) write (debug_unit,'(a)') "mon_to_day"
    if (debug_write_cal_effects) write (debug_unit,*) nt,nyrs    
    
    pi=4.0d0*datan(1.0d0)
    
    ! generate smoothing weights
    do j=1,nw
        jj=j-(nw/2)-1
        x=(dble(jj)/((nw-1)/2.0d0))*4.0d0
        wgt(j)=(1.0d0/2.0d0*pi)*(exp(-0.5d0*x**2))
        if (debug_write_cal_effects) write (debug_unit,*) j,jj,x,wgt(j)
    end do
       
    ! intialize output variables
    xd_out=xfill; xd_temp=xfill
    
    ! main loop 
    n=0; mm=0
    do nn=1,nyrs         
        ! collect nm monthly values to interpolate this year, along with month lengths
        nd=0; nfill=0
        do m=1,nm
            mm=mm+1
            xm(m)=xm_in(mm)
            if (xm(m) .eq. xfill) nfill=nfill+1
            iml(m)=imonlen(mm)
            nd=nd+iml(m)
        end do
        if (debug_write_cal_effects) write (debug_unit,'(i12,i12,i8,12g14.6)') n,nn,mm,(xm(m),m=1,nm)
        if (debug_write_cal_effects) write (debug_unit,'(24x,13i14)') (iml(m),m=1,nm),nd
        
        ! check for fill value in any month, skip whole year if found
        if (nfill .eq. 0) then      
            
            ! mean-preserving daily interpolation
            call hdaily(nm,nd,xm,iml,no_negatives,xdh)           

            ! save daily values
            do i=1,nd
                n=n+1
                xd_out(n)=xdh(i)
            end do
        else        
            do i=1,nd
                n=n+1
                xd_out(n)=xfill
            end do      
        end if
    end do
    
    if (smooth) then
    
        ! save xd_out
        xd_temp = xd_out
        
        ! smooth across years
        n=imonlen(1)+imonlen(2)+imonlen(3)+imonlen(4)+imonlen(5)+imonlen(6) &
            +imonlen(7)+imonlen(8)+imonlen(9)+imonlen(10)+imonlen(11)+imonlen(12)
        mm=12
        if (debug_write_cal_effects) write (debug_unit,'("smooth:  n,mm: ",2i6)') n,mm
        do nn=1,nyrs-1
            jjj=n-(nsw/2)-1       
            do js=1,nsw
                jjj=jjj+1
                wsum=0.0d0; xd_out(jjj)=0.0d0
                jj=jjj-((nw-1)/2)-1
                do j=1,nw
                    jj=jj+1
                    if (xd_temp(jj) .ne. xfill) then
                        xd_out(jjj)=xd_out(jjj)+xd_temp(jj)*wgt(j)
                        wsum=wsum+wgt(j)
                    end if
                end do
                if (wsum.ne.0.0d0) then
                    xd_out(jjj)=xd_out(jjj)/wsum
                else
                    xd_out(jjj)=xfill
                end if

            end do

            n=n+imonlen(mm+1)+imonlen(mm+2)+imonlen(mm+3)+imonlen(mm+4)+imonlen(mm+5)+imonlen(mm+6) &
                +imonlen(mm+7)+imonlen(mm+8)+imonlen(mm+9)+imonlen(mm+10)+imonlen(mm+11)+imonlen(mm+12)
            mm=mm+12
        end do
        
    end if
    
    if (restore) then
    
        ! restore long-term mean
        if (debug_write_cal_effects) write (debug_unit,'("restore")') 
        xm_ltm=0.0d0; xd_ltm=0.0d0; monlentot=0.0d0
        if (no_negatives) then
            do mm=1,nt
                if (xm_in(mm).gt.0.0d0) then
                    if (xm_in(mm).ne.xfill) then
                        xm_ltm=xm_ltm+xm_in(mm)*dble(imonlen(mm))
                        monlentot=monlentot+dble(imonlen(mm))
                    end if
                end if
            end do
            nzero=0
            do n=1,ndtot
                if (xd_out(n).gt.0.0d0) then
                    if (xd_out(n).ne. xfill) then
                        xd_ltm=xd_ltm+xd_out(n)
                        nzero=nzero+1
                    end if
                end if
            end do

            if (monlentot.ne.0.0d0) then
                xm_ltm=xm_ltm/monlentot
            else 
                xm_ltm=xfill
            end if
            if (nzero.ne.0) then
                xd_ltm=xd_ltm/dble(nzero)
            else
                xd_ltm=xfill
            end if
            ltmdiff=xm_ltm-xd_ltm
        else
            do mm=1,nt      
                if (xm_in(mm).ne.xfill) then
                    xm_ltm=xm_ltm+xm_in(mm)*dble(imonlen(mm))
                    monlentot=monlentot+dble(imonlen(mm))
                end if
            end do
            nnonfill=0
            do n=1,ndtot
                if (xd_out(n).ne.xfill) then
                    xd_ltm=xd_ltm+xd_out(n)
                    nnonfill = nnonfill +1
                end if
            end do

            if (monlentot.ne.0.0d0) then
                xm_ltm=xm_ltm/monlentot
            else 
                xm_ltm=xfill
            end if
            if (nnonfill.ne.0) then
                xd_ltm=xd_ltm/dble(nnonfill)
            else
                xd_ltm=xfill
            end if
            ltmdiff=xm_ltm-xd_ltm
        end if
    
        if (no_negatives) then
            do n=1,ndtot
                if (xd_out(n).ne.xfill) then
                    if (xd_out(n).gt.0.0d0) then
                        if (xd_out(n)+ltmdiff.gt.0.0d0) then
                            xd_out(n)=xd_out(n)+ltmdiff
                        end if
                    end if
                end if
            end do
        else
            do n=1,ndtot
                if (xd_out(n).ne.xfill) then
                    xd_out(n)=xd_out(n)+ltmdiff
                end if
            end do
        end if
        
    end if
    
    if (debug_write_cal_effects) write (10,'(10g14.6)') xdh
    
end subroutine mon_to_day_ts

subroutine day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xd,xfill,xm_adj)
! aggregation of pseudo- or actual daily data to months using a paleo calendar

    implicit none
    
    integer(4), parameter       :: nm=12, nd=366
    integer(4), intent(in)      :: ny, ndtot        ! number of years, total number of days
    integer(4), intent(in)      :: ndays(ny)        ! number of days in each year
    real(8), intent(in)         :: rmonbeg(ny,nm), rmonend(ny,nm)   ! beginning and ending days of each month
    real(8), intent(in)         :: xd(ndtot)        ! daily values
    real(8), intent(in)         :: xfill            ! _FillValue
    real(8), intent(out)        :: xm_adj(ny*nm)    ! (aggregated) average monthly values
    
    ! variables used to calculate monthly means
    integer(4)              :: ibegday, iendday             ! beginning day and ending day of each year
    integer(4)              :: ibeg(nm), iend(nm)           ! beginning and ending (integer) day of each month
    integer(4)              :: ndays_in_month(nm)           ! integer number of days in month
    real(8)                 :: xdx(-29:nd+30)               ! daily data for current year, padded by data from adjacent years
    real(8)                 :: wgt(-29:nd+30), wsum         ! weights (for interpolating over fractional days)
    integer(4)              :: nfill                        ! number of days with fill values
    
    integer(4)              :: n, m, i, nn
    
    logical                 :: debug_write_cal_effects = .false.
    
    debug_write_cal_effects = .false.
    if (debug_write_cal_effects) write (10,'(a)') "day_to_mon"
    if (debug_write_cal_effects) write (10,*) ny, ndtot
    
    ! loop over years, collecting daily data for each year, and getting monthly means
    iendday = 0; nn = 0
    xm_adj=0.0d0
    do n=1,ny
        if (debug_write_cal_effects) write (10,'("n, ndays:", 2i5)') n,ndays(n)
        ibegday = iendday + 1
        iendday = ibegday + ndays(n) - 1
        if (debug_write_cal_effects) write (10,*) n,ibegday,iendday
        
        if (ny .eq. 1) then       ! single-year Aclim data  
            ! wrap the input daily data
            xdx(-29:0)=xd(ndays(n)-30+1:ndays(n))
            xdx(1:ndays(n))=xd(1:ndays(n))
            xdx(ndays(n)+1:ndays(n)+30)=xd(1:30)
        else 
            ! copy current year into xdx
            xdx(1:ndays(n)) = xd(ibegday:iendday)
            ! pad beginning and end of xdx
            if (n .eq. 1) then
                xdx(-29:0) = xd(ndays(n)-30+1:ndays(n))
                xdx(ndays(n)+1:ndays(n)+30) = xd(iendday+1:iendday+30)
            elseif (n .eq. ny) then
                xdx(-29:0) = xd(ibegday-30:ibegday-1)
                xdx(ndays(n)+1:ndays(n)+30) = xd(ibegday+1:ibegday+30)
            else
                xdx(-29:0) = xd(ibegday-30:ibegday-1)
                xdx(ndays(n)+1:ndays(n)+30) = xd(iendday+1:iendday+30)
            end if
        end if
    
        ! integer beginning and end of each month, and number of days in each month
        ! ndays_in_month should be equal to the integer month length + 1
        ibeg=ceiling(rmonbeg(n,:)); iend=ceiling(rmonend(n,:)); ndays_in_month=(iend-ibeg+1)
        if (debug_write_cal_effects) write (10,'("rmonbeg:  ",12f14.8)') rmonbeg(n,:)
        if (debug_write_cal_effects) write (10,'("rmonend:  ",12f14.8)') rmonend(n,:)
        if (debug_write_cal_effects) write (10,'("ibeg:  ",12i4)') ibeg
        if (debug_write_cal_effects) write (10,'("iend:  ",12i4)') iend
        if (debug_write_cal_effects) write (10,'("ndays: ",13i4)') ndays_in_month, sum(ndays_in_month)
 
        ! monthly means
        do m=1,nm
            nn = nn + 1
            nfill = 0; wgt=1.0d0; wsum=0.0d0
            wgt(ibeg(m))=abs(rmonbeg(n,m)-dble(ibeg(m)))
            wgt(iend(m))=abs(rmonend(n,m)-dble(iend(m)-1))
            do i=ibeg(m),iend(m)
                if (xdx(i) .ne. xfill) then
                    xm_adj(nn)=xm_adj(nn)+xdx(i)*wgt(i)
                    wsum=wsum+wgt(i)
                    if (debug_write_cal_effects) &
                    write (10,'("m, i, xd(i), xm_adj(nn), wgt(i), wsum: ",2i4,2f12.6,2f12.6)') &
                        m,i,xdx(i),xm_adj(nn),wgt(i),wsum
                else
                    nfill = nfill + 1
                end if
            end do
            if (wsum .ne. 0.0d0 .and. nfill .eq. 0) then
                xm_adj(nn)=xm_adj(nn)/wsum
            else
                xm_adj(nn)=xfill
            end if
            if (debug_write_cal_effects) write(10,'("m, xm_adj(nn): ",i4,2f12.6)') &
                m,xm_adj(nn),sngl(xm_adj(nn))
        end do
    end do 
    if (debug_write_cal_effects) write (10,'(12g14.6)') xm_adj

end subroutine day_to_mon_ts


end module calendar_effects_subs

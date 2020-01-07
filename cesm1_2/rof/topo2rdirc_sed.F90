program topo2rdirc

  implicit none
  include 'netcdf.inc'

! -------------------------------------------------------------------
! Notes
! -------------------------------------------------------------------
! Start with topography and end with river directions map.
! Simpler procedure than Graham et al. (1999), because
! we have no idea where the rivers were in the geologic past.
!
! The output consists of two ascii files: fort.10 and fort.11.
! The former is the river directions map in the format required by CLM.
! The latter lists all grid cells involved in infinite river loops.
! The user must redirect rivers in the vicinity of infinite loops.
! Alternatively, the user may start with topography which always slopes to
! the ocean, doesn't have internal basins, and doesn't contain completely
! flat plateaus.
!
! Currently written for Cretaceous topography.
!
! slevis, feb 2003
! -------------------------------------------------------------------
! Variable declarations
! -------------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

! IF YOU CHOOSE A DIFFERENT RESOLUTION (nlon, nlat),
! YOU MUST ALSO READ A DIFFERENT INPUT FILE (filei).
! integer, parameter :: nlon = 180        !input grid : longitude points
! integer, parameter :: nlat =  90        !input grid : latitude  points
!  integer, parameter :: nlon = 720        !input grid : longitude points
!  integer, parameter :: nlat = 360        !input grid : latitude  points
 integer, parameter :: nlon = 360        !input grid : longitude points
 integer, parameter :: nlat = 180        !input grid : latitude  points
  integer, parameter :: maxsecpass = 5000 !max # of infinite loops

! input file
  character(len=200) :: filei              !input filename
  integer :: ncid                         !netCDF file id for filei
  integer :: lat_id                       !netCDF latitude id
  integer :: lon_id                       !netCDF longitude id
  integer :: topo_id                      !netCDF topography id

! input variables from filei
  real(r8) :: lon(nlon)                   !longitude dimension array
  real(r8) :: lat(nlat)                   !latitude dimension array
  real(r8) :: topo(nlon,nlat)             !topography array

! output related variables
  integer  :: nowat(maxsecpass,2)         !'now at' index
  integer  :: rdirc(nlon,nlat)            !river directions
  integer  :: intbasin(nlon,nlat)         !internal basins
  integer  :: source(nlon,nlat)           !river sources
  real(r8) :: temp(nlon*nlat,3)           !temporary
  real(r8) :: tempo(maxsecpass,4)         !temporary
  real(r8) :: tmpo(4)                     !temporary
  real(r8) :: tmp(nlon*nlat,5)            !temporary

  logical :: dontcount                    !index related to next variable
  integer :: csecpass                     !count 2nd passes => inf. loops
  integer :: secondpass, ocean, ret       !indices
  integer :: ii, jj, i, j, count, line    !indices
  integer :: sort                         !index and
  integer :: locsmallest                  !variables
  integer, dimension(1) :: minlocarray    !used in the
  real(r8) :: smallest                    !sorting loops
  real(r8) :: topomin(8)     !surrounding topography in ascending order

  integer :: iii
  integer :: jjj
  integer :: isinf(nlon,nlat)
  integer :: mindist,mindir,x,found


! -------------------------------------------------------------------
! Read in netcdf topography file
! -------------------------------------------------------------------

! VARIABLE NAME IS topo_depth.
! This must be changed in script. 
  filei = 'input_topo_data'

! data:   lat (1D -89.75 to 89.75), lon (1D 0.25 to 359.75), topo (2D 360x720)
! OR 2x2: lat (1D -89 to 89),       lon (1D 1 to 359),       topo (2D  90x180)

  ret = nf_open (filei, nf_nowrite, ncid)
  if (ret == nf_noerr) then
    write(6,*)'Successfully opened ',trim(filei)
    call wrap_inq_varid (ncid, 'lat', lat_id)
    call wrap_inq_varid (ncid, 'lon', lon_id)
    call wrap_inq_varid (ncid, 'topo', topo_id)        !filei
    call wrap_get_var8 (ncid, lat_id, lat)
    call wrap_get_var8 (ncid, lon_id, lon)
    call wrap_get_var8 (ncid, topo_id, topo)
    write(*,*)'successfully read netcdf topography file'
  else
    write(6,*)'cannot open netcdf topography file'
    call endrun
  end if
  ret = nf_close (ncid)
  write(*,*)'netcdf topography file read successfully'
  write(*,*)'========================================'

  rdirc    = 0  !initialize
  intbasin = 0  !these
  temp     = 0. !four
  tmp      = 0. !variables
  isinf = 0.0

! take care of the edges

do i = 1,nlon                 ! south and north edges
  if (topo(i,1) > 0.) then    ! land points (Graham et al. 1999)
    rdirc(i,1) = 1            ! point north only
  else                        ! ocean points
    rdirc(i,1) = 0
  end if                      ! land mask

  if (topo(i,nlat) > 0.) then ! land points (Graham et al. 1999)
    rdirc(i,nlat) = 5         ! point south only
  else                        ! ocean points
    rdirc(i,nlat) = 0
  end if                      ! land mask
end do                        ! south and north edges

do j = 2,nlat-1               ! east and west edges
  if (topo(1,j) > 0.) then    ! land points (Graham et al. 1999)
    topomin(1) = topo(1,j+1)
    topomin(2) = topo(2,j+1)
    topomin(3) = topo(2,j)
    topomin(4) = topo(2,j-1)
    topomin(5) = topo(1,j-1)
    topomin(6) = topo(nlon,j-1)
    topomin(7) = topo(nlon,j)
    topomin(8) = topo(nlon,j+1)
    do sort = 1,8
      smallest = minval(topomin(sort:8))
      minlocarray = minloc(topomin(sort:8))
      locsmallest = (sort-1) + minlocarray(1)
      topomin(locsmallest) =topomin(sort)
      topomin(sort) = smallest
    end do
    count = 0
    do while (rdirc(1,j) == 0 .and. count < 9)
      count = count + 1
      if     (topo(1,j+1) == topomin(count)) then
        rdirc(1,j) = 1   ! point N
      elseif (topo(2,j+1) == topomin(count)) then
        rdirc(1,j) = 2   ! point NE
      elseif (topo(2,j) == topomin(count)) then
        rdirc(1,j) = 3   ! point E
      elseif (topo(2,j-1) == topomin(count)) then
        rdirc(1,j) = 4   ! point SE
      elseif (topo(1,j-1) == topomin(count) .and. rdirc(1,j-1) /= 1) then
        rdirc(1,j) = 5   ! point S
      elseif (topo(nlon,j-1) == topomin(count) .and. rdirc(nlon,j-1) /= 2) then
        rdirc(1,j) = 6   ! point SW
      elseif (topo(nlon,j) == topomin(count)) then
        rdirc(1,j) = 7   ! point W
      elseif (topo(nlon,j+1) == topomin(count)) then
        rdirc(1,j) = 8   ! point NW
      end if
    end do
    if (rdirc(1,j) == 0) then
      intbasin(1,j) = 1
      pause 'internal basin2'
    end if
  else                        ! ocean points
    rdirc(1,j) = 0
  end if                      ! land mask for western edge of domain

  if (topo(nlon,j) > 0.) then ! land points (Graham et al. 1999)
    topomin(1) = topo(nlon,j+1)
    topomin(2) = topo(1,j+1)
    topomin(3) = topo(1,j)
    topomin(4) = topo(1,j-1)
    topomin(5) = topo(nlon,j-1)
    topomin(6) = topo(nlon-1,j-1)
    topomin(7) = topo(nlon-1,j)
    topomin(8) = topo(nlon-1,j+1)
    do sort = 1,8
      smallest = minval(topomin(sort:8))
      minlocarray = minloc(topomin(sort:8))
      locsmallest = (sort-1) + minlocarray(1)
      topomin(locsmallest) =topomin(sort)
      topomin(sort) = smallest
    end do
    count = 0;
    do while (rdirc(nlon,j) == 0 .and. count < 9)
      count = count + 1
      if (topo(nlon,j+1) == topomin(count)) then
        rdirc(nlon,j) = 1  ! point N
      elseif (topo(1,j+1) == topomin(count)) then
        rdirc(nlon,j) = 2  ! point NE
      elseif (topo(1,j) == topomin(count) .and. rdirc(1,j) /= 7) then
        rdirc(nlon,j) = 3  ! point E
      elseif (topo(1,j-1) == topomin(count) .and. rdirc(1,j-1) /= 8) then
        rdirc(nlon,j) = 4  ! point SE
      elseif (topo(nlon,j-1) == topomin(count) .and. rdirc(nlon,j-1) /= 1) then
        rdirc(nlon,j) = 5  ! point S
      elseif (topo(nlon-1,j-1) == topomin(count)) then
        rdirc(nlon,j) = 6  ! point SW
      elseif (topo(nlon-1,j) == topomin(count)) then
        rdirc(nlon,j) = 7  ! point W
      elseif (topo(nlon-1,j+1) == topomin(count)) then
        rdirc(nlon,j) = 8  ! point NW
      end if
    end do
    if (rdirc(nlon,j) == 0) then
      intbasin(nlon,j) = 1
      pause 'internal basin3'
    end if
  else                        ! ocean points
    rdirc(nlon,j) = 0
  end if                      ! land mask for eastern edge of domain
end do                        ! east and west edges

! loop through all data except the domain's edges

do i = 2,nlon-1
  do j = 2,nlat-1
    if (topo(i,j) > 0.) then  ! land points (Graham et al. 1999)
      topomin(1) = topo(i,j+1)
      topomin(2) = topo(i+1,j+1)
      topomin(3) = topo(i+1,j)
      topomin(4) = topo(i+1,j-1)
      topomin(5) = topo(i,j-1)
      topomin(6) = topo(i-1,j-1)
      topomin(7) = topo(i-1,j)
      topomin(8) = topo(i-1,j+1)
      do sort = 1,8
        smallest = minval(topomin(sort:8))
        minlocarray = minloc(topomin(sort:8))
        locsmallest = (sort-1) + minlocarray(1)
        topomin(locsmallest) =topomin(sort)
        topomin(sort) = smallest
      end do
      count = 0;
      do while (rdirc(i,j) == 0 .and. count < 9)
        count = count + 1
        if     (topo(i,j+1) == topomin(count)) then
          rdirc(i,j) = 1  ! point N
        elseif (topo(i+1,j+1) == topomin(count)) then
          rdirc(i,j) = 2  ! point NE
        elseif (topo(i+1,j) == topomin(count)) then
          rdirc(i,j) = 3  ! point E
        elseif (topo(i+1,j-1) == topomin(count)) then
          rdirc(i,j) = 4  ! point SE
        elseif (topo(i,j-1) == topomin(count) .and. rdirc(i,j-1) /= 1) then
          rdirc(i,j) = 5  ! point S
        elseif (topo(i-1,j-1) == topomin(count) .and. rdirc(i-1,j-1) /= 2) then
          rdirc(i,j) = 6  ! point SW
        elseif (topo(i-1,j) == topomin(count) .and. rdirc(i-1,j) /= 3) then
          rdirc(i,j) = 7  ! point W
        elseif (topo(i-1,j+1) == topomin(count) .and. rdirc(i-1,j+1) /= 4) then
          rdirc(i,j) = 8  ! point NW
        end if
      end do
      if (rdirc(i,j) == 0) then
        intbasin(i,j) = 1
        pause 'internal basin1'
      end if
    else                ! ocean points
      rdirc(i,j) = 0
    end if              ! land mask
  end do                ! nlat
end do                  ! nlon

write(*,*) 'begin error checks'

! Error checks
! Look for arrows pointing at each other.
! Don't worry about the first and last rows (j=1 and j=360 at half degree),
! because they have been hardwired to point north and south respectively.

do i = 2,nlon-1
  do j = 2,nlat-1
    if     (rdirc(i,j) == 1 .and. rdirc(i,j+1) == 5) then
      pause 'error1'
    elseif (rdirc(i,j) == 2 .and. rdirc(i+1,j+1) == 6) then
      pause 'error2'
    elseif (rdirc(i,j) == 3 .and. rdirc(i+1,j) == 7) then
      print*,i,j
      pause 'error3'
    elseif (rdirc(i,j) == 4 .and. rdirc(i+1,j-1) == 8) then
      pause 'error4'
    elseif (rdirc(i,j) == 5 .and. rdirc(i,j-1) == 1) then
      pause 'error5'
    elseif (rdirc(i,j) == 6 .and. rdirc(i-1,j-1) == 2) then
      pause 'error6'
    elseif (rdirc(i,j) == 7 .and. rdirc(i-1,j) == 3) then
      pause 'error7'
    elseif (rdirc(i,j) == 8 .and. rdirc(i-1,j+1) == 4) then
      pause 'error8'
    end if
  end do
end do

do j = 2,nlat-1
  if    (rdirc(1,j) == 1 .and. rdirc(1,j+1) == 5) then
      pause 'error9'
  elseif (rdirc(1,j) == 2 .and. rdirc(2,j+1) == 6) then
      pause 'error10'
  elseif (rdirc(1,j) == 3 .and. rdirc(2,j) == 7) then
      pause 'error11'
  elseif (rdirc(1,j) == 4 .and. rdirc(2,j-1) == 8) then
      pause 'error12'
  elseif (rdirc(1,j) == 5 .and. rdirc(1,j-1) == 1) then
      pause 'error13'
  elseif (rdirc(1,j) == 6 .and. rdirc(nlon,j-1) == 2) then
      pause 'error14'
  elseif (rdirc(1,j) == 7 .and. rdirc(nlon,j  ) == 3) then
      pause 'error15'
  elseif (rdirc(1,j) == 8 .and. rdirc(nlon,j+1) == 4) then
      pause 'error16'
  end if
  if     (rdirc(nlon,j) == 1 .and. rdirc(nlon,j+1) == 5) then
      pause 'error17'
  elseif (rdirc(nlon,j) == 2 .and. rdirc(1,j+1) == 6) then
      pause 'error18'
  elseif (rdirc(nlon,j) == 3 .and. rdirc(1,j  ) == 7) then
      pause 'error19'
  elseif (rdirc(nlon,j) == 4 .and. rdirc(1,j-1) == 8) then
      pause 'error20'
  elseif (rdirc(nlon,j) == 5 .and. rdirc(nlon,j-1) == 1) then
      pause 'error21'
  elseif (rdirc(nlon,j) == 6 .and. rdirc(nlon-1,j-1) == 2) then
      pause 'error22'
  elseif (rdirc(nlon,j) == 7 .and. rdirc(nlon-1,j  ) == 3) then
      print*,nlon,j
      pause 'error23'
  elseif (rdirc(nlon,j) == 8 .and. rdirc(nlon-1,j+1) == 4) then
      pause 'error24'
  end if
end do

write(*,*) 'end error checks'

! reorder lon from (0 to 360) to (-180 to 180) as required
! for this file by CLM

line = 0
do j = 1,nlat
  do i = nlon/2+1,nlon      ! ie, 361 to 720 at half degree
    line = line + 1
    temp(line,1) = lat(j)
    temp(line,2) = lon(i)-360.
    temp(line,3) = rdirc(i,j)
    write(10,*) temp(line,:)
    if (intbasin(i,j) == 1) then
      tmp(line,1) = float(j)
      tmp(line,2) = float(i)-float(nlon)/2.
      tmp(line,3) = lat(j)
      tmp(line,4) = lon(i)-360.
      tmp(line,5) = rdirc(i,j)
      write(12,*) tmp(line,:)
    end if
  end do
  do i = 1,nlon/2           ! ie, 1 to 360 at half degree
    line = line + 1
    temp(line,1) = lat(j)
    temp(line,2) = lon(i)
    temp(line,3) = rdirc(i,j)
    write(10,*) temp(line,:)
    if (intbasin(i,j) == 1) then
      tmp(line,1) = float(j)
      tmp(line,2) = float(i)+float(nlon)/2.
      tmp(line,3) = lat(j)
      tmp(line,4) = lon(i)
      tmp(line,5) = rdirc(i,j)
      write(12,*) tmp(line,:)
    end if
  end do
end do

! Identify infinite loops in 2D river direction map.
! Then fix them by hand. Haven't thought of a method
! to fix them automatically...
!
! Follow all streams to the ocean.
! If pass through same grid cell twice => infinite loop.
! Start from river sources which are land points with no
! neighbors pointing to them.

source = 0              ! initialize
csecpass = 0            ! variables
                        ! ***           TO ALL USERS                ***
do i = 2,nlon-1         ! *** Pls check manually for infinite loops ***
  do j = 2,nlat-1       ! ***      at the edges of the domain       ***
    if (rdirc(i,j) > 0) then
      if (rdirc(i,j+1) /= 5 .and. rdirc(i+1,j+1) /= 6 .and. &
          rdirc(i+1,j) /= 7 .and. rdirc(i+1,j-1) /= 8 .and. &
          rdirc(i,j-1) /= 1 .and. rdirc(i-1,j-1) /= 2 .and. &
          rdirc(i-1,j) /= 3 .and. rdirc(i-1,j+1) /= 4) then
        source(i,j) = 1 ! river source
      end if
    end if
  end do
end do

write(*,*) 'Found river sources. Now follow rivers to the ocean.'

do i = 2,nlon-1
  do j = 2,nlat-1
!    if (source(i,j) == 1) then
      nowat = 0; jj = j; ii = i; secondpass = 0; ocean = 0; count = 1
      do while (secondpass == 0 .and. ocean == 0)
        if (ii < 1) then
          ii = nlon
        elseif (ii > nlon) then
          ii = 1
        end if
        count = count + 1
        nowat(count,1) = jj
        nowat(count,2) = ii
        if     (rdirc(ii,jj) == 1) then
          jj = jj + 1;
        elseif (rdirc(ii,jj) == 2) then
          jj = jj + 1;
          ii = ii + 1;
        elseif (rdirc(ii,jj) == 3) then
          ii = ii + 1;
        elseif (rdirc(ii,jj) == 4) then
          ii = ii + 1;
          jj = jj - 1;
        elseif (rdirc(ii,jj) == 5) then
          jj = jj - 1;
        elseif (rdirc(ii,jj) == 6) then
          jj = jj - 1;
          ii = ii - 1;
        elseif (rdirc(ii,jj) == 7) then
          ii = ii - 1;
        elseif (rdirc(ii,jj) == 8) then
          jj = jj + 1;
          ii = ii - 1;
        elseif (rdirc(ii,jj) == 0) then
          ocean = 1                       ! reached the ocean
        end if
        do sort = 1,count-1
          if (nowat(sort,1) == nowat(count,1) .and. &
              nowat(sort,2) == nowat(count,2)) then
            secondpass = 1;               ! found an infinite loop
            tmpo(1) = jj
            tmpo(2) = ii
            tmpo(3) = lat(jj)
            tmpo(4) = lon(ii)
            if (csecpass == 0) then       ! first infinite loop to be counted
              tempo(1,1) = 0
              tempo(1,2) = 0
              tempo(1,3) = 0
              tempo(1,4) = 0
            end if
            dontcount = .false.

            do line = 1,maxsecpass
              if (tempo(line,1)==tmpo(1) .and. tempo(line,2)==tmpo(2) .and. &
                  tempo(line,3)==tmpo(3) .and. tempo(line,4)==tmpo(4)) then
                dontcount = .true.        ! avoid double counting
              end if
            enddo

              if (isinf(ii,jj) == 1) then
                dontcount = .true.        ! avoid double counting of loops
              end if
        isinf(ii,jj)=1
        jjj=0
        iii=0
        if     (rdirc(ii,jj) == 1) then
          jjj = 1;
        elseif (rdirc(ii,jj) == 2) then
          jjj = 1;
          iii = 1;
        elseif (rdirc(ii,jj) == 3) then
          iii = 1;
        elseif (rdirc(ii,jj) == 4) then
          iii = 1;
          jjj = -1;
        elseif (rdirc(ii,jj) == 5) then
          jjj = -1;
        elseif (rdirc(ii,jj) == 6) then
          jjj = -1;
          iii = -1;
        elseif (rdirc(ii,jj) == 7) then
          iii = -1;
        elseif (rdirc(ii,jj) == 8) then
          jjj = 1;
          iii = -1;
        endif
        if (isinf(ii+iii,jj+jjj) == 1) then
          dontcount = .true.        ! avoid double counting of loops
        end if
        isinf(ii+iii,jj+jjj)=1

        if (.not. dontcount) then
              csecpass = csecpass + 1;    !count infinite loops
              tempo(csecpass,1) = jj      !NB: I have shifted ii & lon(ii)
              tempo(csecpass,2) = ii      !in fort.10 (not fort.11) to go from
              tempo(csecpass,3) = lat(jj) !-180 to 180 degrees (see comment
              tempo(csecpass,4) = lon(ii) !after 'end error checks')
              write(11,*) tempo(csecpass,:)
            end if
          end if
        end do
      end do                              ! while secondpass = 0 AND ocean = 0
!    end if                                ! source
  end do                                  ! j loop
end do                                    ! i loop


do x=1,csecpass

mindir=1
mindist=10000

!OK, first travel north:
count=1
found=0
do while(found == 0 .and. tempo(x,1)+count.le.nlat )
  if ( rdirc(tempo(x,2),tempo(x,1)+count).eq.0) then    
    found=1
    if (count.lt.mindist) then
      mindist=count
      mindir=1  
    endif
  else
    count=count+1
  endif
enddo

!OK, south:
count=1
found=0
do while(found == 0 .and. tempo(x,1)-count.gt.0 )
  if ( rdirc(tempo(x,2),tempo(x,1)-count).eq.0) then    
    found=1
    if (count.lt.mindist) then
      mindist=count
      mindir=5  
    endif
  else
    count=count+1
  endif
enddo

!OK, west:
count=1
found=0
do while(found == 0 .and. tempo(x,2)-count.gt.0 )
  if ( rdirc(tempo(x,2)-count,tempo(x,1)).eq.0) then    
    found=1
    if (count.lt.mindist) then
      mindist=count
      mindir=7  
    endif
  else
    count=count+1
  endif
enddo

!OK, east:
count=1
found=0
do while(found == 0 .and. tempo(x,2)+count.le.nlon )
  if ( rdirc(tempo(x,2)+count,tempo(x,1)).eq.0) then    
    found=1
    if (count.lt.mindist) then
      mindist=count
      mindir=3  
    endif
  else
    count=count+1
  endif
enddo

if (mindir.eq.1) then
rdirc(tempo(x,2),tempo(x,1):tempo(x,1)+(mindist-1))=mindir
endif

if (mindir.eq.5) then
rdirc(tempo(x,2),tempo(x,1)-(mindist-1):tempo(x,1))=mindir
endif

if (mindir.eq.7) then
rdirc(tempo(x,2)-(mindist-1):tempo(x,2),tempo(x,1))=mindir
endif

if (mindir.eq.3) then
rdirc(tempo(x,2):tempo(x,2)+(mindist-1),tempo(x,1))=mindir
endif

enddo


line = 0
do j = 1,nlat
  do i = nlon/2+1,nlon      ! ie, 361 to 720 at half degree
    line = line + 1
    temp(line,1) = lat(j)
    temp(line,2) = lon(i)-360.
    temp(line,3) = rdirc(i,j)
    write(13,*) temp(line,:)
    if (intbasin(i,j) == 1) then
      tmp(line,1) = float(j)
      tmp(line,2) = float(i)-float(nlon)/2.
      tmp(line,3) = lat(j)
      tmp(line,4) = lon(i)-360.
      tmp(line,5) = rdirc(i,j)
      write(15,*) tmp(line,:)
    end if
  end do
  do i = 1,nlon/2           ! ie, 1 to 360 at half degree
    line = line + 1
    temp(line,1) = lat(j)
    temp(line,2) = lon(i)
    temp(line,3) = rdirc(i,j)
    write(13,*) temp(line,:)
    if (intbasin(i,j) == 1) then
      tmp(line,1) = float(j)
      tmp(line,2) = float(i)+float(nlon)/2.
      tmp(line,3) = lat(j)
      tmp(line,4) = lon(i)
      tmp(line,5) = rdirc(i,j)
      write(15,*) tmp(line,:)
    end if
  end do
end do




end program topo2rdirc

!===============================================================================

subroutine handle_error(ret)
  implicit none
  include 'netcdf.inc'
  integer :: ret
  write(6,*) 'BAILING'
  if (ret .ne. nf_noerr) then
     write(6,*) 'NCDERR: ERROR: ',nf_strerror(ret)
     call abort
  endif
  write(6,*) 'BAILING'
end subroutine handle_error

!==============================================================================

subroutine wrap_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret

  ret = nf_inq_varid (nfid, varname, varid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_inq_varid

!=============================================================================

subroutine wrap_get_var8 (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real*8 arr(*)

  integer ret

  ret = nf_get_var_double (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_get_var8

!==============================================================================

subroutine endrun
  implicit none
  include 'netcdf.inc'

  call abort
  stop 999
end subroutine endrun


!===============================================================================




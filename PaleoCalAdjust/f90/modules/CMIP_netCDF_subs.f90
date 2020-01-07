module CMIP_netCDF_subs
! This module contains several subroutines to facilitate reading an existing netCDF file of daily or monthly data
!   and creating a new file (of adjusted data), amending and copying the attributes of the original file.
! get_var_diminfo() gets information on dimension variables in the netCDF file (including the number of dimensions).
! copy_dims_and_glatts() copies dimension, dimension variables and attributes from the old file to the new.
! new_time_day() and new_time_month() redefine the time-dimension variable to reflect the adjusted month lengths.
! define_outvar() defines the output (adjusted) variable.
! check() and current_time() are utility subroutines.

    use netcdf
    use typesizes

    implicit none

    integer, parameter          :: maxdims = 15, maxvars = 150, maxatts = 100, maxdimsize = 1440

    ! dimensions
    integer(4)                  :: ncid_in, ncid_out, ndim, nvar, nglatt, unlimid, ncformat
    integer(4)                  :: dimid(maxdims), dimlen(maxdims)
    character(nf90_max_name)    :: dimname(maxdims)

    ! variables
    integer(4)                  :: varid_in, varid_out
    integer(4)                  :: xtype(maxvars), nvardims(maxvars), vardimids(maxvars,nf90_max_var_dims), nvaratts(maxvars)
    character(nf90_max_name)    :: varname(maxvars)
    character(2048)             :: varatt(maxvars,maxatts)
    integer(4)                  :: dataid(maxvars)
    
    ! input variable dimensions
    integer(4)                  :: invar_ndim
    integer(4)                  :: invar_dimid(maxdims), invar_dimlen(maxdims)
    character(nf90_max_name)    :: invar_dimname(maxdims)

    ! global attributes
    character(nf90_max_name)    :: glattname, addglattname
    character(2048)             :: addglatt
    character(19)               :: current      ! current time

    ! input and output variable names
    character(nf90_max_name)    :: varinname, varoutname
    character(nf90_max_name)    :: addvarattname
    character(2048)             :: addvaratt

    logical                     :: nc_print = .false.

    contains   
    
subroutine get_var_diminfo(ncid_in, invar_name, invar_ndim, invar_dimid, invar_dimlen, invar_dimname)
! gets the dimension characteristics of a variable

    implicit none
    
    ! netCDF file and variable name
    integer(4), intent(in)      :: ncid_in
    character(*), intent(in)    :: invar_name
    
    ! dimensions
    integer(4), intent(out)     :: invar_ndim
    integer(4), intent(out)     :: invar_dimid(maxdims), invar_dimlen(maxdims)
    character(nf90_max_name), intent(out)    :: invar_dimname(maxdims)
    
    integer(4)                  :: invarid, i
    
    if (nc_print) print '("Variable dimensions")'
    call check( nf90_inq_varid(ncid_in, invar_name, invarid) )
    call check( nf90_inquire_variable(ncid_in, invarid, ndims=invar_ndim) )
    if (nc_print) print '("invar_name, invarid, invar_ndim: ",a,1x,i5,i3)', trim(invar_name), invarid, invar_ndim
    
    call check( nf90_inquire_variable(ncid_in, invarid, dimids = invar_dimid(:invar_ndim)) )
    
    do i = 1, invar_ndim
        call check( nf90_inquire_dimension(ncid_in, invar_dimid(i), invar_dimname(i), invar_dimlen(i)) )
        if (nc_print) print '("i, invar_dimid, invar_dimlen, invar_dimname: ",i3,i5,i6,1x,a)', &
        i, invar_dimid(i), invar_dimlen(i), trim(invar_dimname(i))
    end do
    
end subroutine get_var_diminfo

subroutine copy_dims_and_glatts(ncid_in, ncid_out, addglattname, addglatt, nt, time, time_bnds, comment, varid_out)
! copies dimensions and global attributes from an existing to a new netCDF file
! while replacing the old time values

    implicit none

    integer(4), intent(in)      :: ncid_in, ncid_out
    character(*), intent(in)    :: addglattname, addglatt
    integer(4), intent(in)      :: nt
    real(8), intent(in)         :: time(nt), time_bnds(2,nt)
    character(*), intent(in)    :: comment
    integer(4), intent(out)     :: varid_out

    real(8), allocatable        :: var1d(:), var2d(:,:), var3d(:,:,:) ! see select case() statements below for dimensions
    real(8)                     :: var0d
    integer(4)                  :: i, ii
    character(256)              :: attname

    ! structure of input file
    if (nc_print) print '("Define and copy dimensions, dimension variables, and global attributes")'
    call check( nf90_inquire(ncid_in, ndim, nvar, nglatt, unlimid, ncformat) )
    if (nc_print) print '(" ndim, nvar, nglatt, unlimid: ",4i6)', ndim, nvar, nglatt, unlimid

    ! get dimensions
    if (nc_print) print '("Define dimensions:")'
    write (*,*) nt
    do i = 1,ndim
       dimid(i) = i
        call check( nf90_inquire_dimension(ncid_in, dimid(i), dimname(i), dimlen(i)) )
        if (nc_print) print '("   in:  dimid, dimlen, dimname = " ,2i7,1x,a)', i, dimlen(i), trim(dimname(i))
        
        
        if (i .eq. unlimid) then
           call check( nf90_def_dim(ncid_out, dimname(i), nf90_unlimited, dimid(i)) )
        else
            call check( nf90_def_dim(ncid_out, dimname(i), dimlen(i), dimid(i)) )
         end if
        ! if dimension is time, replace length with new length (nt)
        if (dimname(i) .eq. 'time') dimlen(i) = nt
        
    end do
    ! define dimension variables
    if (nc_print) print '("Define dimension variables:")'
    varid_out = 0
    do i = 1,nvar

        varid_in = i
        call check( nf90_inquire_variable(ncid_in, varid_in, varname(i), xtype(i), ndims=nvardims(i), natts=nvaratts(i)) )
        if (nc_print) print '(" i, xtype, nvardims, nvaratts = ", 4i6, 1x, a)', i,xtype(i),nvardims(i),nvaratts(i),trim(varname(i))
        call check( nf90_inquire_variable(ncid_in, varid_in, dimids=vardimids(i,:nvardims(i))) )
        if (nc_print) print '("    vardimids = ", 6i6)', vardimids(i,:nvardims(i))

        ! define dimension variables only
        select case (varname(i))
        case ('lon', 'lon_bnds', 'lat', 'lat_bnds', 'time', 'time_bnds', 'climatology_bnds', 'height', &
             'plev', 'j', 'i', 'vertices', 'lon_vertices', 'lat_vertices', 'lev', 'lev_bnds', 'gw', 'hyam', 'hybm', 'P0', 'ilev', &
             'hyai', 'hybi', 'ndbase', 'nsbase', 'nbdate', 'nbsec', 'mdt', 'date', 'datesec', &
             'ndcur', 'nscur', 'co2vmr', 'ch4vmr', 'n2ovmr', 'f11vmr', 'f12vmr', 'sol_tsi', 'nsteph', &
             'TLON', 'TLAT', 'ULON', 'ULAT', 'NCAT', 'VGRDi', 'VGRDs', 'VGRDb', 'tmask', 'tarea', 'uarea', &
             'dxt', 'dyt', 'dxu', 'dyu', 'HTN', 'HTE', 'ANGLE', 'ANGLET', 'lont_bounds', 'latt_bounds', &
             'lonu_bounds', 'latu_bounds', 'time_bounds', &
!             'moc_components', 'transport_components', 'transport_regions', 'z_t', 'z_t_150m', 'z_w', &
             'z_t', 'z_t_150m', 'z_w', &
             'z_w_top', 'z_w_bot', 'lat_auz_grid', 'moc_z', 'dz', 'dzw', 'ULONG', 'TLONG', 'KMT', 'KMU', &
             'REGION_MASK', 'UAREA', 'TAREA', 'HU', 'HT', 'DXU', 'DYU', 'DXT', 'DYT', 'HUS', 'HUW', &
             'days_in_norm_year', 'grav', 'omega', 'radius', 'cp_sw', 'sound', 'vonkar', 'cp_air', 'rho_air', &
             'rho_sw', 'rho_fw', 'stefan_boltzmann', 'latent_heat_vapor', 'latent_heat_fusion', 'latent_heat_fusion_mks', &
             'ocn_ref_salinity', 'sea_ice_salinity', 'T0_Kelvin', 'salt_to_ppt', 'ppt_to_salt', 'mass_to_Sv', &
             'heat_to_PW', 'salt_to_Svppt', 'salt_to_mmday', 'momentum_factor', 'hflux_factor', 'fwflux_factor', &
             'salinity_factor', 'sflux_factor', 'nsurface_t', 'nsurface_u', 'time_bound')
        
            ! define variable
            varid_out=varid_out+1

            call check( nf90_def_var(ncid_out, varname(i), xtype(i), vardimids(i,:nvardims(i)), dataid(i)) )
            if (nc_print) print '("     defining varid_out : ",i4,1x,a)', varid_out, trim(varname(i))

            ! copy attributes
            do ii=1,nvaratts(i)
                call check( nf90_inq_attname(ncid_in, varid_in, ii, attname) )
                call check( nf90_copy_att(ncid_in, varid_in, attname, ncid_out, varid_out) )
                if (nc_print) print '("  ii, attname: ",i4,1x,a)', ii, trim(attname)
            end do

            ! dimension variable is time, add comment
            if (varname(i) .eq. 'time') then
                call check ( nf90_put_att(ncid_out, varid_out, 'comment', comment) )
                if (nc_print) print '("  ii, attname: ",i4,1x,a)', ii, "comment"
            end if
            

        case default
            continue
        end select
    end do

    ! copy global attributes
    if (nc_print) print '("Copy global attributes")' 
    do i = 1, nglatt
        call check ( nf90_inq_attname(ncid_in, nf90_global, i, glattname) )
        if (nc_print) print '(" global attribute: ", a)', trim(glattname)
        call check ( nf90_copy_att(ncid_in, nf90_global, trim(glattname), ncid_out, nf90_global) )
    end do

    ! add new global attribute
    call check ( nf90_put_att(ncid_out, nf90_global, addglattname, addglatt) )

    ! end definition of new netCDF file (temporarily)
    call check( nf90_enddef(ncid_out) )

    ! copy dimension variable values, replacing time and time_bnds
    if (nc_print) print '("Copy dimension variables")'
    varid_out = 0
    do i = 1,nvar

        varid_in = i
        select case (varname(i))
        case ('lon', 'lon_bnds', 'lat', 'lat_bnds', 'time', 'time_bnds', 'climatology_bnds', 'height', &
             'plev', 'j', 'i', 'vertices', 'lon_vertices', 'lat_vertices', 'lev', 'lev_bnds', 'gw', 'hyam', 'hybm', 'P0', 'ilev', &
             'hyai', 'hybi', 'ndbase', 'nsbase', 'nbdate', 'nbsec', 'mdt', 'date', 'datesec', &
             'ndcur', 'nscur', 'co2vmr', 'ch4vmr', 'n2ovmr', 'f11vmr', 'f12vmr', 'sol_tsi', 'nsteph', &
             'TLON', 'TLAT', 'ULON', 'ULAT', 'NCAT', 'VGRDi', 'VGRDs', 'VGRDb', 'tmask', 'tarea', 'uarea', &
             'dxt', 'dyt', 'dxu', 'dyu', 'HTN', 'HTE', 'ANGLE', 'ANGLET', 'lont_bounds', 'latt_bounds', &
             'lonu_bounds', 'latu_bounds', 'time_bounds', &
!             'moc_components', 'transport_components', 'transport_regions', 'z_t', 'z_t_150m', 'z_w', &
             'z_t', 'z_t_150m', 'z_w', &
             'z_w_top', 'z_w_bot', 'lat_auz_grid', 'moc_z', 'dz', 'dzw', 'ULONG', 'TLONG', 'KMT', 'KMU', &
             'REGION_MASK', 'UAREA', 'TAREA', 'HU', 'HT', 'DXU', 'DYU', 'DXT', 'DYT', 'HUS', 'HUW', &
             'days_in_norm_year', 'grav', 'omega', 'radius', 'cp_sw', 'sound', 'vonkar', 'cp_air', 'rho_air', &
             'rho_sw', 'rho_fw', 'stefan_boltzmann', 'latent_heat_vapor', 'latent_heat_fusion', 'latent_heat_fusion_mks', &
             'ocn_ref_salinity', 'sea_ice_salinity', 'T0_Kelvin', 'salt_to_ppt', 'ppt_to_salt', 'mass_to_Sv', &
             'heat_to_PW', 'salt_to_Svppt', 'salt_to_mmday', 'momentum_factor', 'hflux_factor', 'fwflux_factor', &
             'salinity_factor', 'sflux_factor', 'nsurface_t', 'nsurface_u', 'time_bound')

        
            varid_out=varid_out + 1
            if (nc_print) print '(" varname: ",a)', trim(varname(i))
            ! how many dimensions?
            select case(nvardims(i))
            case (0)
                if (nc_print) print '("  i: nvardims(i): ", 2i6, 1x, a)', &
                    i,nvardims(i),trim(varname(i))
                call check( nf90_get_var(ncid_in, varid_in, var0d) )
                
                ! if variable is time, replace the existing values with new ones
                if (varname(i) .eq. 'time') var1d=time
                call check( nf90_put_var(ncid_out, varid_out, var0d) )
            case (1)
                if (nc_print) print '("  i, nvardims(i), vardimids, dimlen: ", 4i6, 1x, a)', &
                    i, nvardims(i), vardimids(i,1), dimlen(vardimids(i,1)),trim(varname(i))
                allocate(var1d(dimlen(vardimids(i,1))))
                call check( nf90_get_var(ncid_in, varid_in, var1d) )
                
                ! if variable is time, replace the existing values with new ones
                if (varname(i) .eq. 'time') var1d=time
                call check( nf90_put_var(ncid_out, varid_out, var1d) )
                deallocate(var1d)
            case (2)
                if (nc_print) print '("  i, nvardims(i), vardimids, dimlen: ", 6i6, 1x, a)', &
                    i, nvardims(i), vardimids(i,1), vardimids(i,2), dimlen(vardimids(i,1)), dimlen(vardimids(i,2)), trim(varname(i))
                allocate(var2d(dimlen(vardimids(i,1)), dimlen(vardimids(i,2))))
                call check( nf90_get_var(ncid_in, varid_in, var2d) )
                
                ! if variable is time_bnds or climatology_bnds, replace the existing values with new ones
                if (varname(i) .eq. 'time_bnds') var2d=time_bnds
                if (varname(i) .eq. 'time_bounds') var2d=time_bnds
                if (varname(i) .eq. 'time_bound') var2d=time_bnds
                if (varname(i) .eq. 'climatology_bnds') var2d=time_bnds
                call check( nf90_put_var(ncid_out, varid_out, var2d) )
                deallocate(var2d)
            case (3)
                if (nc_print) print '("  i, nvardims(i), vardimids, dimlen: ", 8i6, 1x, a)', &
                    i, nvardims(i), vardimids(i,1), vardimids(i,2), vardimids(i,3), &
                    dimlen(vardimids(i,1)), dimlen(vardimids(i,2)), dimlen(vardimids(i,3)), trim(varname(i))
                allocate(var3d(dimlen(vardimids(i,1)), dimlen(vardimids(i,2)), dimlen(vardimids(i,3))))
                call check( nf90_get_var(ncid_in, varid_in, var3d) )
                
                ! if variable is time_bnds or climatology_bnds, replace the existing values with new ones
                if (varname(i) .eq. 'time_bnds') var2d=time_bnds
                if (varname(i) .eq. 'time_bounds') var2d=time_bnds
                if (varname(i) .eq. 'time_bound') var2d=time_bnds
                if (varname(i) .eq. 'climatology_bnds') var2d=time_bnds
                call check( nf90_put_var(ncid_out, varid_out, var3d) )
                deallocate(var3d)               
                
            case default
                continue
            end select

        case default
            continue
        end select
    end do

    varid_out = varid_out + 1   ! id of next variable to be written

end subroutine copy_dims_and_glatts

subroutine new_time_day(ncid_in, ny, nm, nt, ndtot, &
    imonmid_ts, imonbeg_ts, imonend_ts, ndays_ts, time, time_bnds)

    implicit none

    integer(4), intent(in)      :: ncid_in
    integer(4), intent(in)      :: ny, nm, nt           ! number of years, months and total number of months
    integer(4), intent(in)      :: ndtot                ! total number of days
    integer(4), intent(in)      :: imonmid_ts(nt)       ! month mid-days as time series
    integer(4), intent(in)      :: imonbeg_ts(nt)       ! month beginning days as time series
    integer(4), intent(in)      :: imonend_ts(nt)       ! month ending days as time series
    integer(4), intent(in)      :: ndays_ts(nt)         ! number of days in year

    real(8), intent(out)        :: time(nt), time_bnds(2,nt)                ! new time variables

    integer(4)                  :: timeid               ! variable id
    real(8)                     :: day_time(ndtot)      ! (old) time value for each day
    integer(4)                  :: ndyr                 ! number of days in previous years

    integer(4)                  :: imid, ibeg, iend     ! indices
    integer(4)                  :: n

    write (*,'("new time values: ny, nm, nt, ndtot: ",4i8)') ny,nm,nt,ndtot

    ! get the existing daily time values
    call check ( nf90_inq_varid(ncid_in, 'time', timeid) )
    call check ( nf90_get_var(ncid_in, timeid, day_time) )

    ! new time variables -- copy appropriate existing daily time values
    imid = imonmid_ts(1); ibeg = imonbeg_ts(1); iend = imonend_ts(1)
    if (ibeg .lt. 1) ibeg = 1
    time(1) = day_time(imid)
    time_bnds(1,1) = day_time(ibeg)
    time_bnds(2,1) = day_time(iend)
    do n = 2,nt
        ndyr = ndays_ts(n) - ndays_ts(1)
        imid = imonmid_ts(n); ibeg = imonbeg_ts(n); iend = imonend_ts(n)
        time(n) = day_time(imid + ndyr)
        time_bnds(1,n) = time_bnds(2,n-1)
        time_bnds(2,n) = day_time(iend + ndyr)
    end do

end subroutine new_time_day

subroutine new_time_month(calendar_type, ncid_in, ny, nm, nt, &
     rmonmid_ts, rmonbeg_ts, rmonend_ts, ndays_ts, time, time_bnds)

    ! redefines monthly time variables

    implicit none

    character(*), intent(in)    :: calendar_type
    integer(4), intent(in)      :: ncid_in
    integer(4), intent(in)      :: ny, nm, nt           ! number of years, months and total number of months
    real(8), intent(in)         :: rmonmid_ts(nt)       ! month mid-days as time series
    real(8), intent(in)         :: rmonbeg_ts(nt)       ! month beginning days as time series
    real(8), intent(in)         :: rmonend_ts(nt)       ! month ending days as time series
    integer(4), intent(in)      :: ndays_ts(nt)         ! number of days in year

    real(8), intent(out)        :: time(nt), time_bnds(2,nt)                ! new time variables

    integer(4)                  :: timeid               ! variable id
    real(8)                     :: mon_time(nt)         ! (old) time value for each month
    real(8)                     :: ref_time             ! reference time (e.g. Jan "0" of the first year)
    integer(4)                  :: ndyr                 ! number of days in previous years

    integer(4)                  :: n

    write (*,'("new_time: ny, nm, nt: ",3i8)') ny,nm,nt

    ! get the existing monthly time values
    call check ( nf90_inq_varid(ncid_in, 'time', timeid) )
    call check ( nf90_get_var(ncid_in, timeid, mon_time) )
    
    ! new time variables -- calculate appropriate monthly values
    if (trim(calendar_type) .eq. '360_day') then
        ref_time = mon_time(1) - 15.0
    else
        ref_time = mon_time(1) - 31.0  !mvr  15.5
    end if

!    time(1) = rmonmid_ts(1) + ref_time
    time_bnds(1,1) = rmonbeg_ts(1) + ref_time 
    time_bnds(2,1) = rmonend_ts(1) + ref_time 
    time(1) = time_bnds(2,1)

    do n = 2,nt
        ndyr = ndays_ts(n) - ndays_ts(1)
!        time(n) = rmonmid_ts(n) + ref_time + dble(ndyr)
        time_bnds(1,n) = time_bnds(2,n-1)
        time_bnds(2,n) = rmonend_ts(n) + ref_time + dble(ndyr)
        time(n) = time_bnds(2,n)
    end do

end subroutine new_time_month

subroutine define_outvar(ncid_in, ncid_out, varinname, varid_out, varoutname, addvarattname, addvaratt, varid_in) 

    implicit none

    integer(4), intent(in)      :: ncid_in, ncid_out
    character(*), intent(in)    :: varinname, varoutname
    integer(4), intent(in)      :: varid_out
    character(*), intent(in)    :: addvarattname, addvaratt
    integer(4), intent(out)     :: varid_in 

    integer(4)                  :: i, ii

    character(256)              :: attname

    ! find input variable number
    do i=1,nvar
        if (trim(varname(i)) .eq. trim(varinname)) exit
    end do
    varid_in = i
    if (nc_print) print '(" i, varid, nvardims(i), varname(i), varinname: ",3i4,1x,a,1x,a)', &
        i, varid_in, nvardims(i), trim(varname(i)), trim(varinname)

    ! enter define mode again
    call check( nf90_redef(ncid_out) )

    ! define variable
    call check( nf90_def_var(ncid_out, varoutname, xtype(i), vardimids(i,:nvardims(i)), dataid(i)) )
    if (nc_print) print '(" defining output variable : ",i4,1x,a)', varid_out, trim(varname(i))

    ! copy attributes
    do ii=1,nvaratts(i)
        call check( nf90_inq_attname(ncid_in, varid_in, ii, attname) )
        call check( nf90_copy_att(ncid_in, varid_in, attname, ncid_out, varid_out) )
        if (nc_print) print '("  ii, attname: ",i4,1x,a)', ii, trim(attname)
    end do

    ! add new variable attribute
    call check ( nf90_put_att(ncid_out, varid_out, addvarattname, addvaratt) )

    ! end definition of new netCDF file
    call check( nf90_enddef(ncid_out) )

end subroutine define_outvar

subroutine check(status)
! netCDF error message handler

    use netcdf
    use typesizes

    implicit none

    integer, intent(in) :: status
    if (status.ne.nf90_noerr) then
        print *, status,trim(nf90_strerror(status))
        stop 'stopped on error'
    end if

end subroutine check

subroutine current_time(current)
! gets the current time

    character(19), intent(out) :: current
    character(8)        :: cdate
    character(10)       :: ctime

    call date_and_time(cdate,ctime)

    current = cdate(1:4)//"-"//cdate(5:6)//"-"//cdate(7:8)//" "// &
        ctime(1:2)//":"//ctime(3:4)//":"//ctime(5:6)

end subroutine current_time

end module CMIP_netCDF_subs

from netCDF4 import Dataset
import sys

def nccopy(filein,fileout,unpackshort=True,
    zlib=True,complevel=6,shuffle=True,fletcher32=False,
    clobber=False,lsd_dict=None,nchunk=10,quiet=False,classic=0,
    vars=None,istart=0,istop=-1):
    """convert a netcdf 3 file (filein) to a netcdf 4 file
    The default format is 'NETCDF4', but can be set
    to NETCDF4_CLASSIC if classic=1.
    If unpackshort=True, variables stored as short
    integers with a scale and offset are unpacked to floats.
    in the netcdf 4 file.  If the lsd_dict is not None, variable names
    corresponding to the keys of the dict will be truncated to the decimal place
    specified by the values of the dict.  This improves compression by
    making it 'lossy'..
    If vars is not None, only variable names in the list 
    will be copied (plus all the dimension variables).
    The zlib, complevel and shuffle keywords control
    how the compression is done.

    This code is basically the nc3tonc4 script from python netCDF4 libraray
    with slight modifications.
    """

    ncfilein = Dataset(filein, 'r')
    if classic:
        ncfileout = Dataset(fileout,'w',clobber=clobber,format='NETCDF4_CLASSIC')
    else:
        ncfileout = Dataset(fileout,'w',clobber=clobber,format='NETCDF4')
    
    mval = 1.e30 # missing value if unpackshort=True
    # create dimensions. Check for unlimited dim.
    unlimdimname = False
    unlimdim = None
    # create global attributes.
    if not quiet: sys.stdout.write('copying global attributes ..\n')
    #for attname in ncfilein.ncattrs():
    #    setattr(ncfileout,attname,getattr(ncfilein,attname))
    ncfileout.setncatts(ncfilein.__dict__) 


    if not quiet: sys.stdout.write('copying dimensions ..\n')
    for dimname,dim in ncfilein.dimensions.items():
        if dim.isunlimited():
            unlimdimname = dimname
            unlimdim     = dim
            ncfileout.createDimension(dimname, None)
            if istop == -1: istop = len(unlimdim)
        else:
            ncfileout.createDimension(dimname, len(dim))
   
    # create variables.
    if vars is None:
       varnames = ncfilein.variables.keys()
    else:
       # variables to copy specified
       varnames = vars
       # add dimension variables
       for dimname in ncfilein.dimensions.keys():
           if dimname in ncfilein.variables.keys() and dimname not in varnames:
               varnames.append(dimname)
    
    # Copy variables
    for varname in varnames:
        ncvar = ncfilein.variables[varname]
        if not quiet: sys.stdout.write('copying variable %s\n' % varname)
        # quantize data?
        if lsd_dict is not None and lsd_dict.has_key(varname):
            lsd = lsd_dict[varname]
            if not quiet: sys.stdout.write('truncating to least_significant_digit = %d\n'%lsd)
        else:
            lsd = None # no quantization.
        
        # unpack short integers to floats?
        if unpackshort and hasattr(ncvar,'scale_factor') and hasattr(ncvar,'add_offset'):
            dounpackshort = True
            datatype = 'f4'
        else:
            dounpackshort = False
            datatype = ncvar.dtype
        
        # is there an unlimited dimension?
        if unlimdimname and unlimdimname in ncvar.dimensions:
            hasunlimdim = True
        else:
            hasunlimdim = False
        
        if dounpackshort:
            if not quiet: sys.stdout.write('unpacking short integers to floats ...\n')
            sys.stdout.write('')
        
        if hasattr(ncvar, '_FillValue'):
            FillValue = ncvar._FillValue
        else:
            FillValue = None 
        
        var = ncfileout.createVariable(varname,datatype,ncvar.dimensions, fill_value=FillValue, 
                                                                          least_significant_digit=lsd,
                                                                          zlib=zlib,
                                                                          complevel=complevel,
                                                                          shuffle=shuffle,
                                                                          fletcher32=fletcher32)
        
        # fill variable attributes.
        attdict = ncvar.__dict__
        if '_FillValue' in attdict: del attdict['_FillValue']
        if dounpackshort and 'add_offset' in attdict: del attdict['add_offset']
        if dounpackshort and 'scale_factor' in attdict: del attdict['scale_factor']
        if dounpackshort and 'missing_value' in attdict: attdict['missing_value']=mval
        var.setncatts(attdict)
        #for attname in ncvar.ncattrs():
        #    if attname == '_FillValue': continue
        #    if dounpackshort and attname in ['add_offset','scale_factor']: continue
        #    if dounpackshort and attname == 'missing_value':
        #        setattr(var,attname,mval)
        #    else:
        #        setattr(var,attname,getattr(ncvar,attname))
        # fill variables with data.
        if hasunlimdim: # has an unlim dim, loop over unlim dim index.
            # range to copy
            if nchunk:
                start = istart; stop = istop; step = nchunk
                if step < 1: step = 1
                for n in range(start, stop, step):
                    nmax = n+nchunk
                    if nmax > istop: nmax=istop
                    idata = ncvar[n:nmax]
                    if dounpackshort:
                        tmpdata = (ncvar.scale_factor*idata.astype('f')+ncvar.add_offset).astype('f')
                        if hasattr(ncvar,'missing_value'):
                            tmpdata = NP.where(idata == ncvar.missing_value, mval, tmpdata)
                    else:
                        tmpdata = idata
                    var[n-istart:nmax-istart] = tmpdata
            else:
                idata = ncvar[:]
                if dounpackshort:
                    tmpdata = (ncvar.scale_factor*idata.astype('f')+ncvar.add_offset).astype('f')
                    if hasattr(ncvar,'missing_value'):
                        tmpdata = NP.where(idata == ncvar.missing_value, mval, tmpdata)
                else:
                    tmpdata = idata
                var[0:len(unlimdim)] = tmpdata
        else: # no unlim dim or 1-d variable, just copy all data at once.
            idata = ncvar[:]
            if dounpackshort:
                tmpdata = (ncvar.scale_factor*idata.astype('f')+ncvar.add_offset).astype('f')
                if hasattr(ncvar,'missing_value'):
                    tmpdata = NP.where(idata == ncvar.missing_value, mval, tmpdata)
            else:
                tmpdata = idata
            var[:] = tmpdata
        ncfileout.sync() # flush data to disk
    
    # close files.
    ncfilein.close()
    ncfileout.close()

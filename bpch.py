__all__ = ['bpch']

import os

# part of the default Python distribution
from collections import defaultdict
from warnings import warn

# numpy is a very common installed library
import numpy as np
from numpy import ndarray, fromfile, memmap, dtype, arange, zeros, ceil, diff, concatenate, append, pi, sin

from matplotlib import use
# PseudoNetCDF is my own home grown
# https://dawes.sph.unc.edu/trac/PseudoNetCDF
#from PseudoNetCDF import PseudoNetCDFVariable, PseudoNetCDFFile
class PseudoNetCDFDimension(object):
    """
    Dimension object responds like that of netcdf4-python
    """
    def __init__(self, group, name, size):
        self._len = int(size)
    def isunlimitied(self):
        return False
    def __len__(self):
        return self._len

class PseudoNetCDFFile(object):
    """
    PseudoNetCDFFile provides an interface and standard set of
    methods that a file should present to act like a netCDF file
    using the Scientific.IO.NetCDF.NetCDFFile interface.
    """
    def __new__(cls, *args, **kwds):
        new = object.__new__(cls, *args, **kwds)
        new.variables={}
        new.dimensions={}
        new._ncattrs = ()
        return new
    
    def __init__(self, *args, **properties):
        for k, v in properties.iteritems():
            setattr(self, k, v)

    def __setattr__(self, k, v):
        if not (k[:1] == '_' or k in ('dimensions', 'variables', 'groups')):
            self._ncattrs += (k,)
        object.__setattr__(self, k, v)
    def createDimension(self,name,length):
        """
        name - string name for dimension
        length - maximum length of dimension
        """
        self.dimensions[name]=PseudoNetCDFDimension(self, name, length)

    def createVariable(self, name, type, dimensions, **properties):
        """
        name - string
        type - numpy dtype code (e.g., 'f', 'i', 'd')
        dimensions - tuple of dimension keys that can be
                     found in objects' dimensions dictionary
        """
        var = self.variables[name] = PseudoNetCDFVariable(self,name,type,dimensions, **properties)
        return var

    def close(self):
        """
        Does nothing.  Implemented for continuity with Scientific.IO.NetCDF
        """
        pass

    def ncattrs(self):
        return self._ncattrs

    sync=close
    flush=close

class PseudoNetCDFVariable(ndarray):
    """
    PseudoNetCDFVariable presents the Scientific.IO.NetCDF.NetCDFVariable interface,
    but unlike that type, provides a contructor for variables that could be used
    without adding it to the parent file
    """
    def __setattr__(self, k, v):
        if not hasattr(self, k) and k[:1] != '_':
            self._ncattrs += (k,)
        ndarray.__setattr__(self, k, v)
    def ncattrs(self):
        return self._ncattrs
    def __new__(subtype,parent,name,typecode,dimensions,**kwds):
        """
        Creates a variable using the dimensions as defined in
        the parent object

        parent: an object with a dimensions variable
        name: name for variable
        typecode: numpy style typecode
        dimensions: a typle of dimension names to be used from
                    parrent
        kwds: Dictionary of keywords to be added as properties
              to the variable.  **The keyword 'values' is a special
              case that will be used as the starting values of
              the array

        """
        if 'values' in kwds.keys():
            result=kwds.pop('values')
        else:
            shape=[]
            for d in dimensions:
                dim = parent.dimensions[d]

                # Adding support for netCDF3 dimension objects
                if not isinstance(dim, int):
                    dim = len(dim)
                shape.append(dim)

            result=zeros(shape,typecode)

        result=result[...].view(subtype)

        if hasattr(result, '__dict__'):
            result.__dict__['typecode'] = lambda: typecode
            result.__dict__['dimensions'] = tuple(dimensions)
        else:
            result.__dict__={
                'typecode': lambda: typecode,
                'dimensions': tuple(dimensions)
            }

#        object.__setattr__(result, '_ncattrs', ())

        for k,v in kwds.iteritems():
            setattr(result,k,v)
        return result

    def __array_finalize__(self, obj):
        assert(hasattr(self, '_ncattrs') == False)
        self._ncattrs = ()
        if obj is None: return
        if hasattr(obj, '_ncattrs'):
            for k in obj._ncattrs:
                setattr(self, k, getattr(obj, k))
        
    
    def getValue(self):
        """
        Return scalar value
        """
        return self.item()

    def assignValue(self,value):
        """
        assign value to scalar variable
        """
        self.itemset(value)


# These variables define the binary format of the header blocks
# and are only for internal
_general_header_type = dtype('>i4, S40, >i4, >i4, S80, >i4')
_datablock_header_type = dtype('>i4, S20, 2>f4, >i4, >i4, >i4, >i4, S40, >i4, S40, >f8, >f8, S40, 3>i4, 3>i4, >i4, >i4')
_first_header_size = _general_header_type.itemsize + _datablock_header_type.itemsize

class defaultdictfromkey(defaultdict):
    """
    defaultdictfromkey dynamically produces dictionary items
    for keys, using the default_factor function called
    with the key as an argument
    """

    def __missing__(self, key):
        """
        __missing__(key) # Called by __getitem__ for missing key; pseudo-code:
        if self.default_factory is None: raise KeyError((key,))
        self[key] = value = self.default_factory()
        return value
        """
        return self.default_factory(key)

class defaultdictfromthesekeys(defaultdict):
    """
    defaultdictfromthesekeys dynamically produces dictionary items
    for known keys, using the default_factor function called
    with the key as an argument
    """
    def __init__(self, keys, default_factory = None):
        """
        keys - iterable of keys that default_factory can produce
        default_factory - function that takes a key as an argument
                          to create a dictionary item
        """
        self._keys = set([k for k in keys])
        defaultdict.__init__(self, default_factory)
    
    def __iter__(self):
        for i in self._keys:
            yield i
            
    def iterkeys(self):
        for i in self._keys:
            yield i
    
    def itervalues(self):
        for k in self.iterkeys():
            yield self[k]
    
    def iteritems(self):
        for k in self.iterkeys():
            yield (k, self[k])
    
    def keys(self):
        return [k for k in self]
    
    def __setitem__(self, key, value):
        val = defaultdict.__setitem__(self, key, value)
        self._keys.add(key)
        return val
    
    def __delitem__(self, key):
        self._keys.discard(key)
        return defaultdict.__delitem__(self, key)
    
    def pop(self, key):
        val = defaultdict.pop(self, key)
        self._keys.discard(key)
        return val
        
    def __missing__(self, key):
        """
        __missing__(key) # Called by __getitem__ for missing key; pseudo-code:
        if self.default_factory is None: raise KeyError((key,))
        self[key] = value = self.default_factory()
        return value
        """
        if key in self._keys:
            return self.default_factory(key)
        else:
            raise KeyError("%s not found" % (key, ))

    for k in '__setitem__ __delitem__ pop __iter__ iterkeys itervalues iteritems keys'.split():
        exec('%s.__doc__ = defaultdict.%s.__doc__' % (k, k))

class defaultpseudonetcdfvariable(defaultdictfromthesekeys):
    """
    Overwrites __repr__ function to show variables
    """
    def __repr__(self):
        out = "{"
        for k in self._keys:
            out += "\n  '%s': PseudoNetCDFVariable(...)" % k
        out += "\n}"
        return out

    def __str__(self):
        return self.__repr__()

class _diag_group(PseudoNetCDFFile):
    """
    This object acts as a PseudoNetCDF file that gets data from the parent object.
    """
    def __init__(self, parent, groupname, groupvariables):
        """
        parent - a PseudoNetCDFFile
        groupname - a string describing the group
        """
        template = '%s_%%s' % groupname
        def getvar(key):
            try:
                return parent.variables[template % key]
            except (KeyError, ValueError):
                return parent.variables[key]
        self._parent = parent
        self.variables = defaultpseudonetcdfvariable(list(groupvariables) + ['time', 'lev', 'tau0', 'tau1', 'crs', 'AREA', 'lat', 'lon', 'lat_bnds', 'lon_bnds'], getvar)
    
    def __getattr__(self, key):
        try:
            return object.__getattr__(self, key)
        except AttributeError:
            return getattr(self._parent, key)
    
# This class is designed to operate like a dictionary, but
# dynamically create variables to return to the user
class _tracer_lookup(defaultpseudonetcdfvariable):
    """
    _tracer_lookup: finds geos_chem tracer indices from names and returns
                    netcdf like variable
    """
    def __init__(self, parent, datamap, tracerinfo, diaginfo, keys, noscale = False):
        """
        parent: NetCDF-like object to serve dimensions
        datamap: array of pre-dimensioned and datatyped values dim(tstep, i, j, k)
        tracerinfo: dictionary of tracer data keyed by ordinal
        keys: list of keys to serve
        """
        self._noscale = noscale
        self._tracer_data = tracerinfo
        self._diag_data = diaginfo
        self._memmap = datamap
        self._parent = parent
        self._special_keys = ['time', 'lev', 'tau0', 'tau1', 'crs', 'lat', 'lon', 'lat_bnds', 'lon_bnds']
        self._keys = keys + self._special_keys
        self._example_key = keys[0]
        
    def __missing__(self, key):
        if key in ('lat', 'lat_bnds'):
            yres = self._parent.modelres[1]
            if self._parent.halfpolar == 1:
                data = concatenate([[-90.], arange(-90. + yres / 2., 90., yres), [90.]])
            else:
                data = arange(-90, 90 + yres, yres)
            
            dims = ('lat',)
            dtype = 'i'
            kwds = dict(units = 'degrees north', long_name = key, var_desc = key)
            if key == 'lat':
                data = data[:-1] + diff(data) / 2.
                kwds['bounds'] = 'lat_bnds'
            else:
                dims += ('nv',)
                data = data.repeat(2,0)[1:-1].reshape(-1, 2)
            example = self[self._example_key]
            sj = getattr(example, 'STARTJ', 0)
            data = data[sj:sj + example.shape[2]]
            kwds = dict(standard_name = "latitude", long_name = "latitude", units = "degrees_north", axis = "Y", bounds = "lat_bnds")
        elif key in ('lon', 'lon_bnds'):
            xres = self._parent.modelres[0]
            i = arange(0, 360 + xres, xres)
            data = i - (180 + xres / 2. * self._parent.center180)
            dims = ('lon',)
            dtype = 'i'
            kwds = dict(units = 'degrees east', long_name = key, var_desc = key)
            if key == 'lon':
                data = data[:-1] + diff(data) / 2.
                kwds['bounds'] = 'lon_bnds'
            else:
                dims += ('nv',)
                data = data.repeat(2,0)[1:-1].reshape(-1, 2)
            example = self[self._example_key]
            si = getattr(example, 'STARTI', 0)
            data = data[si:si + example.shape[3]]
            kwds = dict(standard_name = "longitude", long_name = "longitude", units = "degrees_east", axis = "X", bounds = "lon_bnds")
        elif key == 'AREA':
           lon = self['lon']
           xres = self._parent.modelres[0]
           nlon = 360. / xres
           latb = self['lat_bnds']
           Re = self['crs'].semi_major_axis
           latb = append(latb[:, 0], latb[-1, 1])
           latr = pi / 180. * latb
           data = 2. * pi * Re * Re / (nlon) * ( sin( latr[1:] ) - sin( latr[:-1] ) )
           data = data[:, None].repeat(lon.size, 1)
           kwds = dict(units = 'm**2')
           dtype = 'i'
           dims = ('J', 'I')
        elif key == 'crs':
          dims = ()
          kwds = dict(grid_mapping_name = "latitude_longitude",
                      semi_major_axis = 6375000.0,
                      inverse_flattening = 0)
          dtype = 'i'
          data = zeros(1, dtype = dtype)
        elif key == 'time':
            tmp_key = self._keys[0]
            data = self._memmap[tmp_key]['header']['f10'] + .5
            dims = ('time',)
            dtype = 'i'
            kwds = dict(units = 'hours since 0 GMT 1/1/1985', standard_name = key, long_name = key, var_desc = key)
        elif key == 'lev':
            tmp_key = self._keys[0]
            data = arange(self._parent.dimensions['lev'], dtype = 'i')
            dims = ('lev',)
            dtype = 'i'
            kwds = dict(units = 'model layer', standard_name = key, long_name = key, var_desc = key)
        elif key == 'tau0':
            tmp_key = self._keys[0]
            data = self._memmap[tmp_key]['header']['f10']
            dims = ('time',)
            dtype = 'i'
            kwds = dict(units = 'hours since 0 GMT 1/1/1985', standard_name = key, long_name = key, var_desc = key)
        elif key == 'tau1':
            tmp_key = self._keys[0]
            data = self._memmap[tmp_key]['header']['f11']
            dims = ('time',)
            dtype = 'i'
            kwds = dict(units = 'hours since 0 GMT 1/1/1985', standard_name = key, long_name = key, var_desc = key)
        else:
            dtype = 'f'
            header = self._memmap[key]['header'][0]
            sl, sj, si = header['f14'][::-1] - 1
            group = header['f7'].strip()
            offset = self._diag_data.get(group, {}).get('offset', 0)
            ord = header['f8'] + offset
            base_units = header['f9']
            scale = self._tracer_data[ord]['SCALE']
            carbon = self._tracer_data[ord]['C']
            units = self._tracer_data[ord]['UNIT']
            kwds = dict(scale = scale, carbon = carbon, units = units, base_units = base_units, standard_name = key, long_name = key, var_desc = key, coordinates = "time lev lat lon")
            tmp_data = self._memmap[key]['data']
            dims = ('time', 'lev', 'lat', 'lon')
            if len(['lev' in dk_ for dk_ in self._parent.dimensions]) > 1:
                dims = ('time', 'lev%d' % tmp_data.dtype['f1'].shape[0], 'lat', 'lon')
                
            assert((tmp_data['f0'] == tmp_data['f2']).all())
            if self._noscale:
                if scale != 1.:
                    warn("Not scaling variables; good for writing")
                data = tmp_data['f1']
            else:
                data = tmp_data['f1'] * scale
            if any([sl != 0, sj != 0, si != 0]):
                nl, nj, ni = header['f13'][::-1]
                #import pdb; pdb.set_trace()
                #tmp_data = zeros((data.shape[0], self._parent.dimensions['lev'], self._parent.dimensions['lat'], self._parent.dimensions['lon']), dtype = data.dtype)
                #el, ej, ei = data.shape[1:]
                #el += sl
                #ej += sj
                #ei += si
                #tmp_data[:, sl:el, sj:ej, si:ei] = data[:]
                #data = tmp_data
                kwds['STARTI'] = si 
                kwds['STARTJ'] = sj
                kwds['STARTK'] = sl
        return PseudoNetCDFVariable(self._parent, key, dtype, dims, values = data, **kwds)
            
class bpch(PseudoNetCDFFile):
    """
    NetCDF-like class to interface with GEOS-Chem binary punch files
    
    f = bpch(path_to_binary_file)
    dim = f.dimensions[dkey] # e.g., dkey = 'lon'
    
    # There are two ways to get variables.  Directly from
    # the file using the long name
    var = f.variables[vkey] # e.g., vkey = 'IJ-AVG-$_NOx'
    
    # Or through a group using the short name
    g = f.groups[gkey] # e.g., gkey = 'IJ-AVG-$'
    var = g.variables[vkey] # e.g., vkey = 'NOx'

    # The variable returned is the same either way    
    print f.dimensions
    print var.unit
    print var.dimensions
    print var.shape
    
    """

    def __init__(self, bpch_path, tracerinfo = None, diaginfo = None, mode = 'r', timeslice = slice(None), noscale = False):
        """
        bpch_path: path to binary punch file
        tracerinfo: path to ascii file with tracer definitions
        diaginfo: path to ascii file with diagnostic group definitions
        mode : {'r+', 'r', 'w+', 'c'}, optional
         |      The file is opened in this mode:
         |  
         |      +------+-------------------------------------------------------------+
         |      | 'r'  | Open existing file for reading only.                        |
         |      +------+-------------------------------------------------------------+
         |      | 'r+' | Open existing file for reading and writing.                 |
         |      +------+-------------------------------------------------------------+
         |      | 'w+' | Create or overwrite existing file for reading and writing.  |
         |      +------+-------------------------------------------------------------+
         |      | 'c'  | Copy-on-write: assignments affect data in memory, but       |
         |      |      | changes are not saved to disk.  The file on disk is         |
         |      |      | read-only.                                                  |
         |      +------+-------------------------------------------------------------+        
         timeslice: If the file is larger than 2GB, timeslice provides a way to subset results.
                    The subset requested depends on the data type of timeslice:
                        - int: return the a part of the file if it was broken into 2GB chunks (0..N-1)
                        - slice: return the times that correspond to that slice (i.e., range(ntimes)[timeslice])
                        - list/tuple/set: return specified times where each time is in the set (0..N-1)
        """
        self._ncattrs = () 
        self._noscale = noscale
        # Read binary data for general header and first datablock header
        header_block = fromfile(bpch_path, dtype = 'bool', count = _first_header_size)
        
        # combine data for convenience
        header = tuple(header_block[:_general_header_type.itemsize].view(_general_header_type)[0]) + \
                 tuple(header_block[_general_header_type.itemsize:].view(_datablock_header_type)[0])
        
        # Verify that all Fortran unformatted buffers match 
        try:
            assert(header[0] == header[2])
            assert(header[3] == header[5])
        except AssertionError:
            raise ValueError("BPCH Files fails header check")
        
        # Assign data from header to global attributes
        self.ftype = header[1]
        self.toptitle = header[4]
        self.modelname, self.modelres, self.halfpolar, self.center180 = header[7:11]
        dummy, dummy, dummy, self.start_tau0, self.start_tau1, dummy, dim, dummy, dummy = header[13:-1]
        for dk, dv in zip('lon lat lev'.split(), dim):
            self.createDimension(dk, dv)
        self.createDimension('nv', 2)
        tracerinfo = tracerinfo or os.path.join(os.path.dirname(bpch_path), 'tracerinfo.dat')
        if not os.path.exists(tracerinfo):
            tracerinfo = 'tracerinfo.dat'
        if os.path.exists(tracerinfo):
            if os.path.isdir(tracerinfo): tracerinfo = os.path.join(tracerinfo, 'tracerinfo.dat')
            tracer_data = dict([(int(l[52:61].strip()), dict(NAME = l[:8].strip(), FULLNAME = l[9:39].strip(), MOLWT = float(l[39:49]), C = int(l[49:52]), TRACER = int(l[52:61]), SCALE = float(l[61:71]), UNIT = l[72:].strip())) for l in file(tracerinfo).readlines() if l[0] not in ('#', ' ')])
            tracer_names = dict([(k, v['NAME']) for k, v in tracer_data.iteritems()])
        else:
            warn('Reading file without tracerinfo.dat means that names and scaling are unknown')
            tracer_data = defaultdict(lambda: dict(SCALE = 1., C = 1.))
            tracer_names = defaultdictfromkey(lambda key: key)
        
        diaginfo = diaginfo or os.path.join(os.path.dirname(bpch_path), 'diaginfo.dat')
        if not os.path.exists(diaginfo):
            diaginfo = 'diaginfo.dat'
        if os.path.exists(diaginfo):
            if os.path.isdir(diaginfo): diaginfo = os.path.join(diaginfo, 'diaginfo.dat')
            diag_data = dict([(l[9:49].strip(), dict(offset = int(l[:8]), desc = l[50:].strip())) for l in file(diaginfo).read().strip().split('\n') if l[0] != '#'])
        else:
            warn('Reading file without diaginfo.dat loses descriptive information')
            diag_data = defaultdictfromkey(lambda key: dict(offset = 0, desc = key))
            
        if len(tracer_names) == 0 and not isinstance(tracer_names, defaultdictfromkey):
            raise IOError("Error parsing %s for Tracer data")
        file_size = os.stat(bpch_path).st_size
        offset = _general_header_type.itemsize
        data_types = []
        first_header = None
        keys = []
        self._groups = defaultdict(set)
        while first_header is None or \
              offset < file_size:
            header = memmap(bpch_path, offset = offset, shape = (1,), dtype = _datablock_header_type, mode = mode)[0]
            
            group = header[7].strip()
            tracer_number = header[8]
            unit = header[9].strip()
            if not isinstance(diag_data, defaultdictfromkey):
                goffset = diag_data.get(group, {}).get('offset', 0)
                try:
                    tracername = tracer_names[tracer_number + goffset]
                except:
                    # There are some cases like with the adjoint where the tracer does not have
                    # a tracerinfo.dat line.  In this case, the name matches the tracer with that
                    # number (no offset).  The scaling, however, is not intended to be used.
                    # The unit for adjoint, for instance, is unitless.
                    if tracer_number not in tracer_names:
                        tracername = str(tracer_number)
                        tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = 0, UNIT = unit)
                    else:
                        tracername = tracer_names[tracer_number]
                        tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = tracer_data[tracer_number]['C'], UNIT = unit)

            else:
                warn('%s is not in diaginfo.dat; names and scaling cannot be resolved' % group)
                goffset = 0
                tracername = str(tracer_number)
                diag_data[group] = dict(offset = 0, desc = group)
                tracer_data[tracer_number + goffset] = dict(SCALE = 1., C = 1., UNIT = unit)

            self._groups[group].add(tracername)
            if first_header is None:
                first_header = header
            elif (header[7], header[8]) == (first_header[7], first_header[8]):
                break
            dim = header[13][::-1]
            start = header[14][::-1]
            data_type = dtype('>i4, %s>f4, >i4' % str(tuple(dim[:])))
            assert(data_type.itemsize == header[-2])
            data_types.append(data_type)
            keys.append('%s_%s' % (group, tracername))
            offset += _datablock_header_type.itemsize + header[-2]

        time_type = dtype([(k, dtype([('header', _datablock_header_type), ('data', d)])) for k, d in zip(keys, data_types)])
        field_shapes = set([v[0].fields['data'][0].fields['f1'][0].shape for k, v in time_type.fields.iteritems()])
        field_levs = set([s_[0] for s_ in field_shapes])
        field_rows = set([s_[1] for s_ in field_shapes])
        field_cols = set([s_[2] for s_ in field_shapes])
        if len(field_levs) == 1:
            pass
        elif len(field_levs) == 2:
            self.dimensions['lev'] = max(field_levs)
            self.dimensions['srf_lev'] = min(field_levs)
        else:
            field_levs = list(field_levs)
            field_levs.sort()
            for fl in field_levs:
                self.createDimension('lev%d' % fl, fl)

        assert((float(os.path.getsize(bpch_path)) - _general_header_type.itemsize) % time_type.itemsize == 0.)
        # load all data blocks  
        try:
            datamap = memmap(bpch_path, dtype = time_type, offset = _general_header_type.itemsize, mode = mode)
        except OverflowError:
            hdrsize = _general_header_type.itemsize
            items = (2*1024**3-hdrsize) // time_type.itemsize
            if timeslice != slice(None):
                filesize = os.stat(bpch_path).st_size
                datasize = (filesize - hdrsize)
                all_times = range(datasize / time_type.itemsize)
                if isinstance(timeslice, int):
                    timeslice = slice(items * timeslice, items * (timeslice + 1))
                if isinstance(timeslice, (list, tuple, set)):
                    times = timeslice
                else:
                    times = all_times[timeslice]

                outpath = bpch_path + '.tmp.part'
                mint = times[0]
                maxt = times[-1]
                nt = maxt - mint + 1
                
                if nt > items:
                    warn('Requested %d items; only returning %d items due to 2GB limitation' % (nt, items))
                    times = times[:items]

                outfile = file(outpath, 'w')
                infile = file(bpch_path, 'r')
                hdr = infile.read(hdrsize)
                outfile.write(hdr)
                for t in all_times:
                    if t in times:
                        outfile.write(infile.read(time_type.itemsize))
                    else:
                        infile.seek(time_type.itemsize, 1)
                outfile.close()
                #print mint, maxt, nt, nt * time_type.itemsize
                #cmd = 'dd if=%s ibs=%d skip=1 obs=%d | dd of=%s bs=%d skip=%d count=%d' % (bpch_path, hdrsize, time_type.itemsize, outpath, time_type.itemsize, mint, nt)
                #print cmd
                #os.system(cmd)
                datamap = memmap(outpath, dtype = time_type, mode = mode, offset = hdrsize)
            else:
                datamap = memmap(bpch_path, dtype = time_type, shape = (items,), offset = _general_header_type.itemsize, mode = mode)
                warn('Returning only the first 2GB of data')

        
        # Create variables and dimensions
        self.variables = _tracer_lookup(parent = self, datamap = datamap, tracerinfo = tracer_data, diaginfo = diag_data, keys = keys, noscale = self._noscale)
        del datamap
        self.createDimension('time', self.variables['tau0'].shape[0])
        self.groups = dict([(k, _diag_group(self, k, v)) for k, v in self._groups.iteritems()])
        self.Conventions = "CF-1.6"


    def __repr__(self):
        return PseudoNetCDFFile.__repr__(self) + str(self.variables)

def tileplot(f, toplot, vmin = None, vmax = None, xmin = None, xmax = None, ymin = None, ymax = None, title = '', unit = '', maptype = 0, log = False, cmap = None):
    # Example: spatial plotting
    from pylab import figure, colorbar, axis
    from matplotlib.colors import Normalize, LogNorm
    has_map = False
    lat = np.append(f.variables['lat_bnds'][0, 0], f.variables['lat_bnds'][:, 1])
    lon = np.append(f.variables['lon_bnds'][0, 0], f.variables['lon_bnds'][:, 1])
    parallels = arange(lat.min(),lat.max() + 15,15)
    meridians = arange(lon.min(), lon.max() + 30,30)
    fig = figure(figsize = (9,4))
    ax = fig.add_axes([.1, .1, .9, .8])
    m = ax
    if maptype == 0:
        try:
            from mpl_toolkits.basemap import Basemap
    
            # I'm actually not sure what the right projection for this 
            # data is.  So the next few lines might need to change.
            m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                    llcrnrlon=-180 - f.modelres[1] * f.center180 / 2.,\
                    urcrnrlon=180 - f.modelres[1] * f.center180 / 2.,\
                    resolution='c')
            has_map = True
            x, y = np.meshgrid(*m(lon, lat))
        except Exception, e:
            warn(str(e))
            x, y = np.meshgrid(lon, lat)
    
        if has_map:
            m.drawcountries()
            m.drawstates()
            m.drawcoastlines()
            m.drawparallels(parallels)
            m.drawmeridians(meridians)
        x = lon
        y = lat
        ax.set_xticks(meridians[::2])
        ax.set_yticks(parallels[::2])
        ax.set_xlim(meridians.min(), meridians.max())
        ax.set_ylim(parallels.min(), parallels.max())
    elif maptype == 1:
        x = lat
        ax.set_xticks(parallels)
        y = np.arange(toplot.shape[0])
        if 'BXHGHT-$_BXHEIGHT' in f.variables.keys():
            y = np.cumsum(np.append(0, f.variables['BXHGHT-$_BXHEIGHT'].mean(0).mean(1).mean(1))) / 1000.
            y = y[:toplot.shape[0]]
            ax.set_ylabel('km agl')
        else:
            ax.set_ylabel('layer')
        ax.set_xlabel('degrees north')
    elif maptype == 2:
        x = lon
        y = np.arange(toplot.shape[0])
        ax.set_xticks(meridians)
        if 'BXHGHT-$_BXHEIGHT' in f.variables.keys():
            y = np.cumsum(np.append(0, f.variables['BXHGHT-$_BXHEIGHT'].mean(0).mean(1).mean(1))) / 1000.
            y = y[:toplot.shape[0]]
            ax.set_ylabel('km')
        else:
            ax.set_ylabel('layer')
        ax.set_xlabel('degrees east')
    elif maptype == 3:
        from datetime import datetime, timedelta
        x = np.append(f.variables['tau0'][0], f.variables['tau1'][:]) / 24.
        start_date = datetime(1985, 1, 1) + timedelta(hours = f.variables['tau0'][0])
        x = x - x.min()
        y = np.arange(toplot.shape[0])
        if 'BXHGHT-$_BXHEIGHT' in f.variables.keys():
            y = np.cumsum(np.append(0, f.variables['BXHGHT-$_BXHEIGHT'].mean(0).mean(1).mean(1))) / 1000.
            y = y[:toplot.shape[0]]
            ax.set_ylabel('km')
        else:
            ax.set_ylabel('layer')
        ax.set_xticks(np.linspace(x[0], x[-1], min(len(x), len(ax.get_xticks()))))
        ax.set_xlabel('days since %s' % start_date.strftime('%Y%m%dT%H%M%S'))
        toplot = toplot.T
    elif maptype == 4:
        from datetime import datetime, timedelta
        x = lat
        ax.set_xlabel('degrees north')
        ax.set_xticks(parallels)
        y = np.append(f.variables['tau0'][0], f.variables['tau1'][:]) / 24.
        start_date = datetime(1985, 1, 1) + timedelta(hours = f.variables['tau0'][0])
        y = y - y.min()
        ax.set_yticks(np.linspace(y[0], y[-1], min(len(y), len(ax.get_yticks()))))
        ax.set_ylabel('days since %s' % start_date.strftime('%Y%m%dT%H%M%S'))
    elif maptype == 5:
        from datetime import datetime, timedelta
        x = lon
        ax.set_xlabel('degrees east')
        ax.set_xticks(meridians)
        y = np.append(f.variables['tau0'][0], f.variables['tau1'][:]) / 24.
        start_date = datetime(1985, 1, 1) + timedelta(hours = f.variables['tau0'][0])
        y = y - y.min()
        ax.set_yticks(np.linspace(y[0], y[-1], min(len(y), len(ax.get_yticks()))))
        ax.set_ylabel('days since %s' % start_date.strftime('%Y%m%dT%H%M%S'))
    
    poly = m.pcolor(x, y, toplot, cmap = cmap)
    ax.collections[-1].set_norm((LogNorm if log else Normalize)(vmin = vmin, vmax = vmax))
    cb = colorbar(poly, ax = ax)
    cb.ax.set_xlabel(toplot.units.strip())
    if vmax is None:
        vmax = toplot.max()
    if vmin is None:
        vmin = np.ma.masked_values(toplot, 0).min()
    if log:
        if (np.log10(vmax) - np.log10(vmin)) < 4:
            ticks = np.logspace((np.log10(vmin)), (np.log10(vmax)), 10.)
            cb.set_ticks(ticks)
            cb.set_ticklabels(['%.1f' % x for x in ticks])
    axis('tight')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(title)
    return fig

def getvar(path_to_test_file = '', group_key = None, var_key = None):
    # Example: file open and variable aquisition
    f = None
    while f is None:
        try:
            f = bpch(path_to_test_file)
        except:
            path_to_test_file = raw_input('Enter path to a valid GEOS-Chem file\n:')
    if group_key is None:
        group_names = f.groups.keys()
        if len(group_names) == 1:
            group_key, = group_names
        else:
            group_names.sort()
            while group_key not in group_names:
                group_key = raw_input('Enter a group name: %s\n:' % ', '.join(group_names))
    g = f.groups[group_key]
    
    if var_key is None:
        var_names = [k for k in g.variables.keys() if k not in ['crs', 'AREA', 'lat_bnds', 'lon', 'lon_bnds', 'tau0', 'tau1', 'lev', 'time', 'lat']]
        if len(var_names) == 1:
            var_key, = var_names
        else:
            var_names.sort()
            var_names = map(str, var_names)
            while var_key not in var_names:
                var_key = raw_input('Enter a variable name: %s\n:' % ', '.join(map(str, var_names)))

    var = g.variables[int(var_key) if var_key.isdigit() else var_key]
    return f, group_key, var_key, var

def pad(nplots, option, option_key, default):
    nopts = len(option)
    if nplots != nopts:
        if nopts == 1:
            option = option * nplots
        else:
            if nopts != 0:
                warn("Padding %s with default layer = %s" % (option_key, str(default)))
            
            option = option + [default] * (nplots - nopts)
    return option

def reduce_dim(toplot, eval_str, axis):
    if eval_str == 'animate':
        return toplot
    if eval_str in dir(np):
        toplot = getattr(np, eval_str)(toplot[:], axis = axis)[(slice(None),) * axis + (None,)]
    elif eval_str.isdigit():
        eval_idx = int(eval_str)
        toplot = np.take(toplot, [eval_idx], axis = axis)
    else:
        eval_idx = eval(eval_str)
        if isinstance(eval_idx, tuple) and len(eval_idx) == 3:
            toplot = getattr(np, eval_idx[0])(toplot[(slice(None),) * axis + (slice(eval_idx[1], eval_idx[2]),)], axis = axis)[(slice(None),) * axis + (None,)]
        else:
            toplot = toplot[(slice(None),) * axis + (eval_idx,)]
    return toplot

def run():
    from optparse import OptionParser
    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    parser = MyParser(description = """
bpch.py provides scriptable plotting and acts as a NetCDF-like library for Binary Punch files. 

For plotting, the a binary punch file must be provided. If 

For use as a library, use "from bpch import bpch" in a python script. For more information, on the bpch reader execute help(bpch) after importing as described above.

""",
    epilog="""
If -d and -p are not provided, a variable (specified by -v)
from a group (specified by -g) from each bpchpath (1...N) 
is plotted in 2 dimensions after using -t, -l, -r, and -c to 
reduce other dimensions.

Options -n -x -t -l -r -c and --title can be specified once 
for each for each bpchpath. If only one value is provided, it 
is repeated for each bpchpath

-t -l -r -c (or --time, --layer, --row, --col) are all dimension 
options. If an integer is provided, a single slice of that 
dimensions is taken using 0-based index (i.e., the first time 
is time 0). If a function is provided, it will be applied to 
that dimension (e.g., mean, median, max, min, std...). If a non-function 
value is provided, it will be converted to a python object 
and be used to acquire a "slice." A special case is "animate" which will
produce mp4 files instead of png files. "animate" requires that you have
a working installation of ffmpeg.


**Time is processed before layer, layer is processed before row, 
row is processed before col. If -t is std and -t is std:
    plotvals = timestd(layerstd(var))
  
Examples:  
    
    1. This example produce a Lat-Lon, time average, mean 
       layer with a log color-scale.

       $ python bpch.py -g IJ-AVG-$ -v O3  -t mean -l mean --log ctm.bpch
       Successfully created ctm.bpch_IJ-AVG_O3_timemean_layermean_rowall_colall.png


    2. This example produces a zonal mean, time average plot 
       stopping at 20km (if BOXHEIGHT available) or the 21st layer.
 
       $ python bpch.py -g IJ-AVG-$ -v O3 -t mean -c mean --ymax 20 ctm.bpch
       Successfully created ctm.bpch_IJ-AVG_O3_timemean_layerall_rowall_colmean.png

    
    3. This example produces a 1st layer Latitude-Time Hovmoller Diagram. 
    
       $ python bpch.py -g IJ-AVG-$ -v O3 -l 0 -c mean ctm.bpch
       Successfully created ctm.bpch_IJ-AVG_O3_timeall_layer0_rowall_colmean.png


    4. This example produces a 1st layer Longitude-Time Hovmoller Diagram. 
    
       $ python bpch.py -g IJ-AVG-$ -v O3 -l 0 -r mean ctm.bpch
       Successfully created ctm.bpch_IJ-AVG_O3_timeall_layer0_rowmean_colall.png

    
    5. This example would produce two Ox figures. Both from 
       time 1 (default), but the first file from layer 1 and
       the second file from layer 2. The first figure has a 
       minimum (max) value of 20 (60) and the second has a 
       minimum (maximum) of 25 (65).

       $ python bpch.py -g IJ-AVG-$ -v O3 -n 20 -x 60 -t 0 -l 0 -n 25 -x 65 -t 0 -l 1 ctm.bpch ctm.bpch2
       Successfully created ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png
       Successfully created ctm.bpch2_O3_time0_layer1_rowall_colall.png

       Successfully created ctm.bpch1_Ox_time0_layer0.png
       Successfully created ctm.bpch2_Ox_time0_layer1.png    

    6. This example would produce one Ox difference figure with 
       a minimum of -2 and a maximum of 2.

       $ python bpch.py -d -g IJ-AVG-$ -v O3 -n -2 -x 2 -t 0 -l 0 ctm.bpch ctm.bpch2
       Successfully created ctm.bpch-ctm.bpch2-diff_O3_time0_layer0_rowall_colall.png

    7. This example would produce one O3 animation across time

       $ python bpch.py -d -g IJ-AVG-$ -v O3 -t animate -l 0 ctm.bpch
       Successfully created ctm.bpch_O3_timeanimation_layer0_rowall_colall.mp4
""")
    parser.set_usage("Usage: python bpch.py [-dpnx] [-t TIME] [-l LAYER] [-r ROW] [-c COL] [-g GROUP] [-v VARIABLE] [--title=TITLE] [bpchpath1 [bpchpath2 ... bpchpathN]]")
    parser.add_option("-d", "--difference", dest="difference",action="store_true",default=False,
                        help="Plot (bpchpath1 - bpchpath2)")

    parser.add_option("-p", "--percent-difference", dest="percent",action="store_true",default=False,
                        help="Plot bpchpath1 / bpchpath2 or, if used with -d, (bpchpath1 - bpchpath2) / bpchpath2 * 100.")
    
    parser.add_option("", "--title", dest = "title", type = "string", help = "Title for plot.", action = "append", default = [])
    parser.add_option("-n", "--min", dest = "vmin", type = "float", action = 'append', default = [], help = "Minimum for the color-scale")
    parser.add_option("-x", "--max", dest = "vmax", type = "float", action = 'append', default = [], help = "Maximum for the color-scale")
    parser.add_option("", "--xmin", dest = "xmin", type = "float", action = 'append', default = [], help = "Minimum for the x-axis")
    parser.add_option("", "--xmax", dest = "xmax", type = "float", action = 'append', default = [], help = "Maximum for the x-axis")
    parser.add_option("", "--ymin", dest = "ymin", type = "float", action = 'append', default = [], help = "Minimum for the y-axis")
    parser.add_option("", "--ymax", dest = "ymax", type = "float", action = 'append', default = [], help = "Maximum for the y-axis")
    parser.add_option("", "--backend", dest = "backend", type = "string", default = 'Agg', help = "Set the backend for use. See matplotlib for more details.")
    parser.add_option("", "--log", dest = "log", action = "store_true", default = False, help = "Put data on a log scale.")
    parser.add_option("-g", "--group", dest = "group", default = None,
                        help = "bpch variables are organized into groups whose names are defined in diaginfo.dat; for a list simply do not provide the group and you will be prompted.")

    parser.add_option("-v", "--variable", dest = "variable", default = None,
                        help = "bpch variables have names defined in tracerinfo.dat; for a list simply do not provide the variable and you will be prompted.")
    
    parser.add_option("-t", "--time", dest = "time", type = "string", action = "append", default = [],
                        help = "bpch variables have time intervals and plotting currently only supports single time slices, or aggregate operations (mean, sum, median, min, max, etc).")

    parser.add_option("-l", "--layer", dest = "layer", type = "string", action = "append", default = [],
                        help = "bpch variables have eta-sigma pressure layers and plotting currently only supports single layer slices, or aggregate operations (mean, sum, median, min, max, etc).")

    parser.add_option("-r", "--row", dest = "row", type = "string", action = "append", default = [],
                        help = "bpch variables has rows for latitudes and plotting currently only supports single row slices, or aggregate operations (mean, sum, median, min, max, etc).")

    parser.add_option("-c", "--column", dest = "column", type = "string", action = "append", default = [],
                        help = "bpch variables has rows for latitudes and plotting currently only supports single row slices, or aggregate operations (mean, sum, median, min, max, etc).")

    parser.add_option("", "--format", dest="format", type="string", default="png", help = "Sets the output format (png, pdf, etc.)")

    parser.add_option("", "--mhchem", dest="mhchem", action = "store_true", default = False, help = "Use the mhchem latex options for typesetting.")

    parser.add_option("", "--latex", dest="latex", action = "store_true", default = False, help = "Use the latex options for typesetting.")

    parser.add_option("", "--dpi", dest="dpi", type="int", default=300, help = "Sets the dots per inch for output figures.")

    parser.add_option("", "--cmap", dest="cmap", type="string", default='jet', help = "Sets cmap property of tile plots. See matplotlib for more details. Good options (jet, bwr)")

    parser.add_option("", "--frames-per-second", dest="fps", type="float", default=0.5, help = "Used in combination with animate on time, layer, row, or column to regulate the number of frames per second.")


    (options, args) = parser.parse_args()
    use(options.backend)
    from pylab import rcParams
    rcParams['text.usetex'] = options.latex | options.mhchem
    nfiles = len(args)
    if nfiles == 0:
        parser.print_help()
        exit()
    if options.difference or options.percent:
        fpath1 = args[0]
        fpath2 = args[1]        
        f, group_key, var_key, var1 = getvar(fpath1, group_key = options.group, var_key = options.variable)
        f2, group_key2, var_key2, var2 = getvar(fpath2, group_key = options.group, var_key = options.variable)
        var_key = var1.long_name.strip() + ' - ' + var2.long_name.strip()
        if options.difference:
            var = var1[:] - var2[:]
        else:
            var = var1[:]
        if options.percent:
            var = var / var2[:] * 100
            if options.difference:
                var.units = 'pctdiff'
            else:
                var.units = 'pct'
        
        if options.difference and options.percent:
            keyappend = 'pctdiff'
        elif options.difference:
            keyappend = 'diff'
        elif options.percent:
            keyappend = 'pct'
        
        plotfvars = [('%s-%s-%s' % (os.path.basename(fpath1), os.path.basename(fpath2), keyappend), f, '', var_key, var)]
    else:
        plotfvars = [(fpath,) + getvar(fpath, group_key = options.group, var_key = options.variable) for fpath in args]
    nplots = len(plotfvars)
    if all(map(lambda x: x == [], (options.time, options.layer, options.row, options.column))):
        warn("Using defaults for reducing dimensions: time=mean, layer = 0 (first)")
        options.time = ['mean']
        options.layer = ['0']
    times = pad(nplots, options.time, "time", "slice(None)")
    layers = pad(nplots, options.layer, "layer", "slice(None)")
    cols = pad(nplots, options.column, "column", "slice(None)")
    rows = pad(nplots, options.row, "row", "slice(None)")
    vmins = pad(nplots, options.vmin, "vmin", None)
    vmaxs = pad(nplots, options.vmax, "vmax", None)
    xmins = pad(nplots, options.xmin, "xmin", None)
    xmaxs = pad(nplots, options.xmax, "xmax", None)
    ymins = pad(nplots, options.ymin, "ymin", None)
    ymaxs = pad(nplots, options.ymax, "ymax", None)
    titles = pad(nplots, options.title, "titles", None)
    
    for time_str, layer_str, row_str, col_str, vmin, vmax, xmin, xmax, ymin, ymax, title_str, (fpath, f, group_key, var_key, var) in zip(times, layers, rows, cols, vmins, vmaxs, xmins, xmaxs, ymins, ymaxs, titles, plotfvars):
        try:
            fig_path = ('%s_%s_%s_time%s_layer%s_row%s_col%s.%s' % (os.path.basename(fpath), group_key, var_key, time_str, layer_str, row_str, col_str, options.format)).replace('-$', '').replace('$', '').replace(' ', '').replace('slice(None)', 'all')
            
            toplot = reduce_dim(var[:], time_str, axis = 0)
            toplot = reduce_dim(toplot, layer_str, axis = 1)
            toplot = reduce_dim(toplot, row_str, axis = 2)
            toplot = reduce_dim(toplot, col_str, axis = 3)
            if toplot.shape[2] != 1 and toplot.shape[3] != 1:
                maptype = 0
            elif toplot.shape[1] != 1 and toplot.shape[2] != 1:
                maptype = 1
            elif toplot.shape[1] != 1 and toplot.shape[3] != 1:
                maptype = 2
            elif toplot.shape[0] != 1 and toplot.shape[1] != 1:
                maptype = 3
            elif toplot.shape[0] != 1 and toplot.shape[2] != 1:
                maptype = 4
            elif toplot.shape[0] != 1 and toplot.shape[3] != 1:
                maptype = 5
            else:
                print toplot.shape
            
            anim_idx = None
            if 'animate' in (time_str, layer_str, row_str, col_str):
                anim_idx = [time_str, layer_str, row_str, col_str].index('animate')
            elif toplot.squeeze().ndim > 2:
                raise UserWarning("Cannot create plot with more than 2 dimensions;\n\tUse the -t, -l, -r, -c options to reduce dimensions.")
            else:
                toplot = toplot.squeeze()
            var_label = var_key
            if options.mhchem:
                var_label = r'\ce{%s}' % var_label
                
            if title_str is None:
                title_str = ('%s (Time: %s, Layer: %s, Row: %s, Col: %s)' % (var_label, time_str, layer_str, row_str, col_str)).replace('$', r'\$').replace('_', ' ').replace('slice(None)', 'all')
            if anim_idx is not None:
                from pylab import figure, draw
                from matplotlib.animation import FuncAnimation
                fig_path = fig_path.replace('.png', '.mp4')
                fig = tileplot(f, toplot[(slice(None),) * anim_idx + (0,)].squeeze(), maptype = maptype, vmin = vmin, vmax = vmax, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, title = title_str, log = options.log, cmap = options.cmap)
                def getframe(i):
                    fig.axes[0].collections[-1].set_array(toplot[(slice(None),) * anim_idx + (i,)].ravel())
                    draw()
                try:
                    animo = FuncAnimation(fig, getframe, frames = np.arange(toplot.shape[anim_idx]), interval = 1, blit = True)
                    animo.save(fig_path, codec = 'libmpeg', fps = .25)
                except IOError, e:
                    print str(e)
                    print 
                    print "-----------------------------"
                    print "Do you have ffmpeg installed?"
                    if raw_input("Should I continue creating %s files for each frame (y/n)?\n:" % options.format).lower() == 'y':
                        for i in np.arange(toplot.shape[anim_idx]):
                            getframe(i)
                            fig.savefig(fig_path.replace('.mp4', '.%d.%s' % (i, options.format)), dpi = options.dpi)
            else:
                fig = tileplot(f, toplot, maptype = maptype, vmin = vmin, vmax = vmax, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, title = title_str, log = options.log, cmap = options.cmap)
                fig.savefig(fig_path, dpi = options.dpi)
                
            f.close()
            print("Successfully created %s" % fig_path)
        except IOError, e:
            print("Unable to produce test figure (maybe you don't have matplotlib or basemap); " + str(e))
    

if __name__ == '__main__':
    run()

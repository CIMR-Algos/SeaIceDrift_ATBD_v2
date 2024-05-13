'''
Purpose: to remap the OSI SAF ice edge product into the ice drift grid.
Two ice masks (one for each day) are needed for the ice drift processing.

The polar observation hole is filled with ice.

Based on processIceMask.c by Thomas Lavergne.
'''

import os
import re
import argparse
from argparse import RawDescriptionHelpFormatter
from datetime import datetime, timedelta
import numpy as np
from numpy import ma
import cartopy.crs as ccrs
from pyresample import parse_area_file, utils, AreaDefinition, kd_tree
from netCDF4 import Dataset, date2num

fill_values = {np.float32:-1e10, np.int16:-32767, np.int32:2147483647,
               int:-32767}


def parse_args():

    p = argparse.ArgumentParser("process_ice_mask",
        formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-i', '--infile', required=True,
                   help="Name for the ice edge or conc file to be "
                   "processed")
    p.add_argument('-o', '--outname', required=False, default=None,
                   help="Full filepath to the output ice edge file, or "
                   "directory path to the output dir")
    p.add_argument('-f', '--gridfile', required=False,
                   default='grids_py.def',
                   help="Path to the gridfile")
    p.add_argument('-g', '--gridname', required=False, default=None,
                   help="Name of the grid to output. If not set will "
                   "use 125 EASE2")

    args = p.parse_args()

    return args


def nc_read(ncfile, var):

    ncdata = {}

    if isinstance(var, list):
        v = var[0]
    else:
        v = var

    # Reading from the NetCDF file
    with Dataset(ncfile, 'r') as dataset:
        ncdata['grid_mapping'] = dataset.variables[v].__dict__['grid_mapping']
        try:
            ncdata['proj4_string'] = dataset.variables[
                ncdata['grid_mapping']].__dict__['proj4_string']
        except:
            ncdata['proj4_string'] = dataset.variables[
                ncdata['grid_mapping']].__dict__['proj4']
        try:
            ncdata['proj_dict'] = utils._proj4.proj4_str_to_dict(
                ncdata['proj4_string'])
        except:
            ncdata['proj_dict'] = utils.proj4.proj4_str_to_dict(
                ncdata['proj4_string'])
        ncdata['xc'] = dataset['xc'][:]
        ncdata['yc'] = dataset['yc'][:]

        try:
            ncdata['fv'] = dataset.variables[v].__dict__['_FillValue']
        except:
            pass

        if isinstance(var, list):
            varlist = var
        else:
            varlist = [var]
        for item in varlist:
            vardata = dataset[item][:]
            if len(vardata.shape) == 3:
                vardata = vardata[0, :, :]
            else:
                vardata = vardata[:, :]

            # NOTE: Be very careful with the fill value here. Trying to
            # use ncdata['fv'] as the fill_value for a status array such
            # as 'flag' means that flag values of 0 (i.e. nominal) are
            # masked out
            if 'fv' in ncdata.keys() and item not in ['statusflag',
                                                      'status_flag', 'flag']:
                ncdata[item] = ma.array(vardata, fill_value=ncdata['fv'])
            else:
                ncdata[item] = ma.array(vardata)

            if var in ['statusflag', 'status_flag', 'flag']:
                ncdata[item].mask = None
            else:
                ncdata[item].mask = ncdata[item].data == ncdata[item].fill_value

        ncdata['lon'] = dataset.variables['lon'][:]
        ncdata['lat'] = dataset.variables['lat'][:]

        # Try fetching time info
        try:
            try:
                d0 = datetime.strptime(dataset.start_date,'%Y-%m-%d %H:%M:%S')
                d1 = datetime.strptime(dataset.stop_date,'%Y-%m-%d %H:%M:%S')
            except:
                try:
                    d0 = datetime.strptime(dataset.start_date_and_time,
                                           '%Y-%m-%dT%H:%M:%SZ')
                    d1 = datetime.strptime(dataset.end_date_and_time,
                                           '%Y-%m-%dT%H:%M:%SZ')
                except:
                    d0 = datetime.strptime(dataset.time_coverage_start,
                                           '%Y-%m-%dT%H:%M:%SZ')
                    d1 = datetime.strptime(dataset.time_coverage_end,
                                           '%Y-%m-%dT%H:%M:%SZ')
            ncdata['sdate'] = d0
            ncdata['edate'] = d1
            ncdata['tspan_hours'] = (d1_00 - d0_00).total_seconds() / (60.*60.)
        except:
            pass

        ncdata['time'] = dataset.variables['time'][0]
        ncdata['time_bnds0'] = dataset.variables['time_bnds'][0][0]
        ncdata['time_bnds1'] = dataset.variables['time_bnds'][0][1]

    # Grid spacing
    sorted_xc = np.sort(ncdata['xc'])
    sorted_yc = np.sort(ncdata['yc'])
    smallest_xc = sorted_xc[0]
    second_smallest_xc = sorted_xc[1]
    smallest_yc = sorted_yc[0]
    second_smallest_yc = sorted_yc[1]
    ncdata['ax'] = second_smallest_xc - smallest_xc
    ncdata['ay'] = second_smallest_yc - smallest_yc

    # Area definitions and extents
    if abs(float(ncdata['xc'][0])) < 10000.:
        sf = 1000.
    else:
        sf = 1.
    ncdata['area_extent'] = (float(ncdata['xc'][0] * sf),
                             float(ncdata['yc'][-1] * sf),
                             float(ncdata['xc'][-1] * sf),
                             float(ncdata['yc'][1] * sf))
    ncdata['area_def'] = AreaDefinition('data', 'data', 'data',
                                        ncdata['proj_dict'],
                                        ncdata['xc'].shape[0],
                                        ncdata['yc'].shape[0],
                                        ncdata['area_extent'])
    ncdata['data_crs'] = ncdata['area_def'].to_cartopy_crs()

    if ncdata['grid_mapping'] in ['Polar_Stereographic_Grid',
                                  'projection_stere']:
        data_globe = ccrs.Globe(semimajor_axis=ncdata['proj_dict']['a'],
                                semiminor_axis=ncdata['proj_dict']['b'])
        if ncdata['lat'][0, 0] > 0:
            ncdata['data_ccrs'] = ccrs.NorthPolarStereo(central_longitude=-45.0,
                                                        globe=data_globe)
            ncdata['hemi'] = 'nh'
        else:
            ncdata['data_ccrs'] = ccrs.SouthPolarStereo(central_longitude=0.0,
                                                        globe=data_globe)
            ncdata['hemi'] = 'sh'
    elif ncdata['grid_mapping'] in ['LambertAzimuthalEqualArea',
                                    'Lambert_Azimuthal_Equal_Area',
                                    'Lambert_Azimuthal_Grid',
                                    'projection_laea']:
        if ncdata['lat'][0, 0] > 0:
            ncdata['data_ccrs'] = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=90,
                false_easting=0, false_northing=0)
            ncdata['hemi'] = 'nh'
        else:
            ncdata['data_ccrs'] = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=-90,
                false_easting=0, false_northing=0)
            ncdata['hemi'] = 'sh'
    else:
        raise ValueError("Unrecognised grid mapping {}".format(
            ncdata['grid_mapping']))

    return ncdata


def write_icemask(indata, outfile, gridname, mode='edge'):

    with Dataset(outfile, 'w') as dataset:

        # dimension and auxiliary datasets
        dimx = dataset.createDimension('xc', len(indata['new_xc']))
        dimy = dataset.createDimension('yc', len(indata['new_yc']))
        dimt = dataset.createDimension('time', 1)

        # The "if" clause deals with "+no_defs" etc
        pdict = dict([el.split('=')
                      for el in indata['new_area_def'].proj4_string.split()
                      if (len(el.split('=')) == 2)])
        if re.search('ease', indata['new_area_def'].area_id):
            crsname = 'Lambert_Azimuthal_Grid'
            gridmapname = 'lambert_azimuthal_grid'
        else:
            crsname = 'Polar_Stereographic_Grid'
            gridmapname = 'polar_stereographic'

        var_crs = dataset.createVariable('crs', np.int32)
#        var_crs.crs_wkt = indata['data_ccrs'].crs_wkt
        if '+lat_ts' in pdict.keys():
            var_crs.standard_parallel = np.float32(pdict['+lat_ts'])
        if '+a' in pdict.keys():
            var_crs.semi_major_axis = np.float32(pdict['+a'])
        if '+b' in pdict.keys():
            var_crs.semi_minor_axis = np.float32(pdict['+b'])
            crs_type = 'b'
        elif '+rf' in pdict.keys():
            var_crs.inverse_flattening = np.float32(pdict['+rf'])
            crs_type = 'rf'
        if '+datum' in pdict.keys():
            var_crs.reference_ellipsoid_name = pdict['+datum']
        var_crs.reference_ellipsoid_name = "WGS 84"
        var_crs.longitude_of_prime_meridian = 0.
        var_crs.prime_meridian_name = "Greenwich"
        var_crs.geographic_crs_name = "unknown"
        var_crs.horizontal_datum_name = "World Geodetic System 1984"
        var_crs.projected_crs_name = "unknown"
        var_crs.grid_mapping_name = gridmapname
        var_crs.latitude_of_projection_origin = np.float32(pdict['+lat_0'])
        var_crs.longitude_of_projection_origin = np.float32(pdict['+lon_0'])
        # Note - original had
        #var_crs.straight_vertical_longitude_from_pole = np.float32(pdict['+lon_0'])
        var_crs.false_easting = 0.0
        var_crs.false_northing = 0.0
        var_crs.proj4_string = str(indata['new_area_def'].proj4_string)
        var_crs.area_id = gridname

        tunits = "seconds since 1978-01-01 00:00:00"
        var_time = dataset.createVariable('time', np.float64, ('time',))
        var_time.axis = "T"
        var_time.units = tunits
        var_time.calendar = "standard"
        var_time.standard_name = "time"
        var_time.long_name = "reference time of map"
        var_time[:] = date2num(indata['edate'], tunits)

        var_xc = dataset.createVariable('xc', np.float64, ('xc'))
        var_xc.axis = "X"
        var_xc.units = "km"
        var_xc.long_name = "x coordinate of projection (eastings)"
        var_xc.standard_name = "projection_x_coordinate"
        var_xc[:] = indata['new_xc'] #* unitconv

        var_yc = dataset.createVariable('yc', np.float64, ('yc'))
        var_yc.axis = "Y"
        var_yc.units = "km"
        var_yc.long_name = "y coordinate of projection (northings)"
        var_yc.standard_name = "projection_y_coordinate"
        var_yc[:] = indata['new_yc'] #* unitconv

        # Data

        var_ice_edge = dataset.createVariable('ice_edge', int,
                                              ('time', 'yc', 'xc'),
                                              fill_value=int(fill_values[int]))
        var_ice_edge.grid_mapping = crsname
        var_ice_edge.coordinates = "time yc xc"
        var_ice_edge[:] = indata['new_ice_edge']

        # Global Attributes
        dataset.Conventions = "CF-1.8"
        dataset.geospatial_lat_min = np.nanmin(indata['new_lat'])
        dataset.geospatial_lat_max = np.nanmax(indata['new_lat'])
        dataset.geospatial_lon_min = np.nanmin(indata['new_lon'])
        dataset.geospatial_lon_max = np.nanmax(indata['new_lon'])
        dataset.time_coverage_start = "{:%Y-%m-%dT%H:%M:%SZ}".format(indata['sdate'])
        dataset.time_coverage_end = "{:%Y-%m-%dT%H:%M:%SZ}".format(indata['edate'])
        dataset.instruments = "ice_{}".format(mode)



def process_ice_mask(infile, outname, gridfile, gridname):

    # The original code processed three types of input file:
    # 1. ice_edge, new (CF) netcdf format
    # 2. ice_conc, new (CF) netcdf format
    # 3. ice_edge, old (MERSEA) netcdf format (currently not implemented here)

    # Reading the input file
    if 'edge' in infile:
        mode = 'edge'
        indata = nc_read(infile, ['ice_edge', 'status_flag'])
    elif 'conc' in infile:
        mode = 'conc'
        indata = nc_read(infile, ['ice_conc', 'status_flag'])
    else:
        raise ValueError("Unknown file type {}, should contain 'ede' or "
                         "'conc'".format(infile))
    # If the mode is conc, convert this to ice edge
    if mode == 'conc':
        datashape = indata['ice_conc'].shape
        indata['ice_edge'] = ma.array(np.zeros(datashape, dtype=int))
        indata['ice_edge'][indata['ice_conc'] < 0.] = -1
        open_water = np.logical_and(indata['ice_conc'] >= 0,
                                    indata['ice_conc'] < 30)
        indata['ice_edge'][open_water] = 0 #1
        open_ice = np.logical_and(indata['ice_conc'] >= 30,
                                  indata['ice_conc'] < 70)
        indata['ice_edge'][open_ice] = 3 #2
        close_ice = np.logical_and(indata['ice_conc'] >= 70,
                                   indata['ice_conc'] < 120)
        indata['ice_edge'][close_ice] = 3
        msk = indata['ice_conc'] < 0.
        indata['ice_edge'][msk] = 0
        indata['ice_edge'].mask = msk

        # This doesn't work because the field has a multiplication factor
        #indata['ice_edge'][indata['ice_conc'] == indata['fv']] = 0
#        indata['ice_edge'][indata['ice_conc'].mask] = int(indata['fv'])
#        indata['ice_edge'][indata['ice_edge'] == 2] = 3

#    elif mode == 'edge':
#        icey = np.logical_or(indata['ice_edge'] == 1,
#                             indata['ice_edge'] == 2)
#        indata['ice_edge'][icey] = 3
#        indata['ice_edge'][indata['ice_edge'] == 0] = 1

    # Processing the status flags
    if mode == 'conc':
        flag_land = 100
        flag_coast = 102
        flag_miss = 101
    elif mode == 'edge':
        flag_land = 1
        flag_coast = 100000 # Not used here
        flag_miss = 8192
    indata['ice_edge'][indata['status_flag'] == flag_land] = 9
    indata['ice_edge'][indata['status_flag'] == flag_coast] = 10
    indata['ice_edge'][indata['status_flag'] == flag_miss] = 0

    # Reading the grid info from the file
    if gridname is None:
        gridname = '{}-ease2-125'.format(indata['hemi'])
    if not (os.path.exists(gridfile) and os.path.isfile(gridfile)):
        raise FileNotFoundError("Grid definition file {} is not found."
                                "".format(gridfile))
    indata['new_area_def'] = parse_area_file(gridfile, gridname)[0]

    # Remap the data
    new_ice_edge = kd_tree.resample_nearest(indata['area_def'],
                                            indata['ice_edge'],
                                            indata['new_area_def'],
                                            radius_of_influence=30000,
                                            fill_value=int(indata['fv']),
                                            reduce_data=False)
    indata['new_ice_edge'] = new_ice_edge.copy()

    indata['new_lon'], indata['new_lat'] = indata['new_area_def'].get_lonlats()
    indata['new_xc'] = indata['new_area_def'].projection_x_coords
    indata['new_yc'] = indata['new_area_def'].projection_y_coords

    # Fill the pole hole
    notice = np.logical_or(indata['new_ice_edge'] == 0,
                           indata['new_ice_edge'] >= 9)
    polefill = np.logical_and(indata['new_lat'] > 85, notice)
    indata['new_ice_edge'][polefill] = 3

    # Apply a median masking to undefined output data points
    edshape = indata['new_ice_edge'].shape
    edgefill = indata['new_ice_edge'].ravel()
    for i, ed in enumerate(edgefill):
        if ed == indata['fv']:
            if i == 0:
                edgefill[i] = 1
            else:
                edgefill[i] = edgefill[i - 1]
    indata['new_ice_edge'] = edgefill.reshape(edshape)

    # Figuring out the output filename
    if os.path.isfile(outname):
        outfile = outname
    elif os.path.isdir(outname):
        outf = 'icemask-multi-{}-{:%Y%m%d}12.nc'.format(gridname,
                        indata['edate'] - timedelta(seconds=60))
        outfile = os.path.join(outname, outf)
    else:
        raise ValueError("Output name is neither a file or a directory")

    # Writing the icemask file out
    write_icemask(indata, outfile, gridname, mode=mode)


def main():

    args = parse_args()

    infile = args.infile
    outname = args.outname
    gridfile = args.gridfile
    gridname = args.gridname

    process_ice_mask(infile, outname, gridfile, gridname)


if __name__ == '__main__':

    main()

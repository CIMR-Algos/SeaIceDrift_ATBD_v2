'''Short code to copy a concentration file from thredds and change the dates in the file'''

import os
from netCDF4 import Dataset, date2num, num2date
import numpy as np
from datetime import datetime, timedelta
import argparse


def parse_args():

    parser = argparse.ArgumentParser(description='cp_and_date_change_iceconc')
    parser.add_argument('-i', '--infile',
                        help="Input ice concentration file")
    parser.add_argument('-o', '--outdir',
                        help="Directory to copy file to")
    parser.add_argument('-d', '--newdate',
                        help="New date to give data in file")

    args = parser.parse_args()

    return args


def cp_and_date_change_iceconc(infile, outdir, newdate):

    # Open the input dataset for reading
    with Dataset(infile, 'r') as dataset:

        # Find the date of the original file in python format
        ots = dataset['time'][:][0]
        tunits = dataset['time'].units
        if tunits.startswith('1978'):
            tdel = (datetime(1978, 1, 1) - datetime(1970, 1, 1)).seconds()
            ots = ots - tdel
        ots_cf = num2date(ots, tunits)
        odatetime = datetime.strptime(str(ots_cf),'%Y-%m-%d %H:%M:%S')
        odate0 = odatetime.replace(hour=0)

        # Calculate the new datetimes
        ndatetime = newdate + timedelta(hours=12)
        ndate0 = newdate
        ndate24 = newdate + timedelta(days=1)

        # Create the new filename
        origfile = os.path.basename(infile)
        newfile = origfile.replace(datetime.strftime(odate0, '%Y%m%d'),
                                   datetime.strftime(ndate0, '%Y%m%d'))
        newfile = os.path.join(outdir, newfile)

        # Open the output dataset for writing
        with Dataset(newfile, 'w') as out:

            # Copy dimensions
            for dim in dataset.dimensions:
                #print("Copying dimension {}".format(dim))
                ret = out.createDimension(dim, dataset.dimensions[dim].size)

            # Copy variables
            for var in dataset.variables:
                #print("Copying variable {}".format(var))
                # Creating the variable
                ncas = dataset[var].ncattrs()
                if '_FillValue' in ncas:
                    myvar = out.createVariable(var,
                                               dataset[var][:].dtype,
                                               dataset[var].dimensions,
                                               fill_value=dataset[var]._FillValue,
                                               zlib=True)
                else:
                    myvar = out.createVariable(var,
                                               dataset[var][:].dtype,
                                               dataset[var].dimensions,
                                               zlib=True)
                # Copying the variable attributes
                for myatt in ncas:
                    if myatt not in ['_FillValue']:
                        setattr(myvar, myatt, dataset[var].__dict__[myatt])
                # Copying the data itself
                # Changing the time variables
                if var == 'time':
                    myvar[:] = date2num(ndatetime, tunits)
                elif var == 'time_bnds':
                    myvar[:] = [date2num(ndate0, tunits),
                                date2num(ndate24, tunits)]
                else:
                    myvar[:] = dataset[var][:]

            # Copy metadata
            #print("Copying the global attributes")
            for gatt in dataset.ncattrs():
                # Setting the dates
                if gatt == 'time_coverage_start':
                    setattr(out, gatt,
                            datetime.strftime(ndate0, "%Y-%m-%dT%H:%M:%SZ"))
                elif gatt == 'time_coverage_end':
                    setattr(out, gatt,
                            datetime.strftime(ndate24, "%Y-%m-%dT%H:%M:%SZ"))
                else:
                    setattr(out, gatt, dataset.__dict__[gatt])

    print("Written a version of {} to {}".format(infile, newfile))
    return(newfile)


if __name__ == '__main__':

    args = parse_args()
    infile = args.infile
    outdir = args.outdir
    newdate = datetime.strptime(args.newdate, '%Y%m%d')

    outname = cp_and_date_change_iceconc(infile, outdir, newdate)

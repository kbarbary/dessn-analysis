#!/usr/bin/env python
"""Run nested sampling on a set of SN lightcurves. FILE is the name of a
FITS format file containing the data. START and END define which lightcurves 
in the file to analyze."""

import os
import sys
from optparse import OptionParser

import numpy as np
from astropy.io.misc import fnpickle
import fitsio
import sncosmo
from sncosmo.photdata import standardize_data
from sncosmo.fitting import _nest_lc

from modeldefs import models
from conf import pikdir

parser = OptionParser(usage='%prog FILE START END',
                      description=__doc__)
opts, args = parser.parse_args(sys.argv[1:])
if len(args) != 3:
    parser.print_help()
    exit()
datafilename, startcand, endcand = args
startcand = int(startcand)
endcand = int(endcand)

f = fitsio.FITS(datafilename)

i = 0
ncut = 0
for meta in f[1][startcand:endcand]:
    i += 1

    print i, 'SNID:', meta['snid'], '------------------------------------'

    pikname = '{0}/{1}.pik'.format(pikdir, meta['snid'])
    if os.path.exists(pikname): continue
    allresults = {}

    data = f[2][meta['datastart']:meta['dataend']]

    # Trim data based on status
    mask = data['status'] < 8
    data = data[mask]

    # Standardize data
    data = standardize_data(data)

    # Get min and max data times
    dtmin = np.min(data['time'])
    dtmax = np.max(data['time'])

    # Cut LC based on S/N ratio
    mask = data['flux'] / data['fluxerr'] > 5.
    if (len(np.unique(data['band'][mask])) < 2):
        print "Fails SNR cut"
        ncut += 1
        continue

    # get mwebv
    mwebv = sncosmo.get_ebv_from_map((meta['ra'], meta['dec']),
                                     mapdir='/home/kyle/Data/dust')

    # Evaluate evidence for each model
    for name, m in models.iteritems():

        # Set t0 bounds (its OK that we keep overwriting this dict)
        m['model'].set(z=0.)
        t0off = m['model'].get('t0') - m['model'].mintime() # t0 offset from
                                                            # mintime

        t0min = dtmin - 30. + t0off
        t0max = dtmax - 30. + t0off
        m['bounds']['t0'] = (t0min, t0max)

        # Set mwebv of model.
        m['model'].set(mwebv=mwebv)

        print name, m['type']
        res = _nest_lc(data, m['model'], m['param_names'],
                       bounds=m['bounds'], tied=m['tied'], nobj=50,
                       verbose=True)

        # Set parameters to those of the model (presumably these are 
        # high-likelihood values), but we're mainly doing this to record
        # the non-varied model parameters (such as mwebv)
        res.param_dict = dict(zip(m['model'].param_names,
                                  m['model'].parameters))

        # record prior and type of the model.
        res.mprior = m['mprior']
        res.type = m['type']

        allresults[name] = res

    fnpickle(allresults, pikname)

f.close()

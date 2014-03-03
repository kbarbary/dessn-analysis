#!/usr/bin/env python

import os
import numpy as np
import fitsio
from astropy.io.misc import fnunpickle
import sncosmo

# Import models used in running the typer
from modeldefs import models

resultdir = 'results'
datafname = 'y1lc_1000_20140225.fits'

f = fitsio.FITS(datafname)

for meta in f[1][0:10]:
    pikname = '{0}/{1}.pik'.format(resultdir, meta['snid'])
    if not os.path.exists(pikname): continue

    # Candidate results
    candres = fnunpickle(pikname)

    # Data for this candidate
    data = f[2][meta['datastart']:meta['dataend']]
    mask = data['status'] < 8
    data = data[mask]  # Trim data based on mask

    # Loop over models, calculate chisq / dof
    print "\n",meta['snid']
    for name, res in candres.iteritems():
        res.chisqdof = (-2. * res.loglmax) / res.ndof
        print "    {} : chisq/dof = {:7.1f}".format(name, res.chisqdof)

    # Sort modelnames by chisqdof
    modelnames = sorted(candres.keys(), key=lambda x: candres[x].chisqdof)

    # Make plot for min chisq/dof
    name = modelnames[0]
    model = models[name]['model']
    model.set(**res.param_dict)
    res = candres[name]
    idx = np.argmax(res.likelihoods)
    parameters = res.samples[idx]
    model.set(**dict(zip(res.param_names, parameters)))
    figtext = "SNID = {}\nmodel = {}\n\\chi^2/dof = {}".format(
        meta['snid'], name, res.chisqdof)
    sncosmo.plot_lc(data, model, figtext=figtext,
                    fname="plots/{}.png".format(meta['snid']))

f.close()

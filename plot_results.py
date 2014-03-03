#!/usr/bin/env python

import os
import numpy as np
import fitsio
from astropy.io.misc import fnunpickle
import sncosmo

# Import models used in running the typer
from modeldefs import models
from conf import datafilename, pikdir, plotdir, startcand, endcand

f = fitsio.FITS(datafilename)

for meta in f[1][startcand:endcand]:
    pikname = '{0}/{1}.pik'.format(pikdir, meta['snid'])
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
        res.chisq = -2. * res.loglmax
        res.chisqdof = res.chisq / res.ndof
        print "    {} : chisq/dof = {:7.1f}".format(name, res.chisqdof)

    # Sort modelnames by chisqdof
    modelnames = sorted(candres.keys(), key=lambda x: candres[x].chisqdof)

    # Make plot for min chisq/dof model
    name = modelnames[0]
    model = models[name]['model']
    res = candres[name]
    model.set(**res.param_dict)
    idx = np.argmax(res.logl)
    parameters = res.samples[idx]
    model.set(**dict(zip(res.param_names, parameters)))
    figtext = "SNID = {}\nmodel = {}\ntype={}\n$\\chi^2$/dof = {:.1f}/{:d}".format(
        meta['snid'], name, res.type, res.chisq, res.ndof)
    sncosmo.plot_lc(data, model, figtext=figtext,
                    fname="{0}/{1}.png".format(plotdir, meta['snid']))

f.close()

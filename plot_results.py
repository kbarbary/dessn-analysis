#!/usr/bin/env python
"""Plot results of nested sampling. FILE is the name of a
FITS format file containing the data. START and END define which lightcurves 
in the file to analyze."""

import os
import sys
from optparse import OptionParser

import numpy as np
import fitsio
from astropy.io.misc import fnunpickle
import sncosmo
import matplotlib.pyplot as plt

from modeldefs import models
from conf import pikdir, plotdir

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

# Initialize chisq dist
minchisq = []
minchisqtype = []
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
    name = modelnames[0]

    minchisq.append(candres[name].chisqdof)
    minchisqtype.append(candres[name].type)

    # Make plot for min chisq/dof
    plotname = "{0}/{1}.png".format(plotdir, meta['snid'])
    if not os.path.exists(plotname):
        model = models[name]['model']
        res = candres[name]
        model.set(**res.param_dict)
        idx = np.argmax(res.logl)
        parameters = res.samples[idx]
        model.set(**dict(zip(res.param_names, parameters)))
        figtext = ("SNID = {}\n"
                   "model = {}\n"
                   "type={}\n"
                   "$\\chi^2$/dof = {:.1f}/{:d}"
                   .format(meta['snid'], name, res.type, res.chisq, res.ndof))
        sncosmo.plot_lc(data, model, figtext=figtext, fname=plotname)

# Plot chisq distribution
minchisq = np.array(minchisq)
minchisqtype = np.array(minchisqtype)
minchisqsets = []
minchisqlabels = []
for t in np.unique(minchisqtype):
    mask = minchisqtype == t
    minchisqsets.append(minchisq[mask])
    minchisqlabels.append(t)

plt.hist(minchisqsets, bins=40, range=(0, 20), histtype='barstacked',
         rwidth=1.0, label=minchisqlabels)
plt.legend()
plt.xlabel("lowest $\chi^2$ / dof")
plt.ylabel("# Candidates")
plt.title("{:d} Candidate light curves".format(len(minchisq)))
plt.savefig("chisqdist.png")

f.close()

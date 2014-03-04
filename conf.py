import os

pikdir = 'pik'
plotdir = 'plots'

# Make directories
for d in [pikdir, plotdir]:
    if not os.path.exists(d): os.mkdir(d)

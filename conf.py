import os

datafilename = 'y1lc_1000_20140225.fits'
pikdir = 'pik'
plotdir = 'plots'
startcand = 0
endcand = 40

# Make directories
for d in [pikdir, plotdir]:
    if not os.path.exists(d): os.mkdir(d)

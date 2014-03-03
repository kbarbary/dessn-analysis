dessn-analysis
==============


DES supernova light curve analysis using [sncosmo](http://github.com/sncosmo/sncosmo).

## Get light curve data

This example command retrieves 1000 candidate light curves from Y1
(starting at 2013-08-30) and saves them in FITS format:

```
get-des-lightcurves -n 1000 --format fits --mindate=2013-08-30 \
		    --bandnames=desg,desr,desi,desz -o output.fits
```


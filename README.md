dessn-analysis
==============


DES supernova light curve analysis using
[sncosmo](http://github.com/sncosmo/sncosmo).

### Get light curve data

This example command retrieves 1000 candidate light curves from Y1
(starting at 2013-08-30) and saves them in FITS format:

```
get-des-lightcurves -n 1000 --format fits --mindate=2013-08-30 \
		    --bandnames=desg,desr,desi,desz -o output.fits
```

### Run nested sampling

```
run_nest.py output.fits 0 100
```

This produces a Python pickle file for each of the first 100
candidates in the fits file. Each pickle contains the full nested sampling
chain for each of the models compared.

### Compile results

```
plot_results.py output.fits 0 100
```
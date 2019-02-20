Detrending Periodic Signals (dips)
==================================

_dips_ is an algorithm for detrending timeseries of strictly periodic signals. It does not assume any functional form for the signal or the background or the noise; it disentangles the strictly periodic component from everything else. We use it in astronomy for detrending _Kepler_, _K2_ and _TESS_ timeseries of periodic variable stars, eclipsing binary stars, exoplanets etc. The algorithm is described in detail in Prsa et al. (2019), PASP, in review -- the reference will be updated shortly.

This repository contains a python 3 implementation of _dips_.

Pre-requisites
--------------

To run _dips_, you will need:

* python 3
* numpy
* scipy

Installation
------------

The `dips` program is available from pip. To install, run `pip3 install dips` for a local install, or `sudo pip3 install dips` for a global install.

If you prefer to install `dips` manually, grab the tarball from github, extract it and run `python3 setup.py install` in the top-level `dips` directory (local install), or `sudo python3 setup.py install` in the top-level `dips` directory (global install).

Running _dips_
--------------

The _dips_ program is run from the command line. It takes a filename with the timeseries as input and computes the disentangled synchronous and asynchronous components of the signal as output. The disentangling process is iterative and might take an appreciable amount of time, depending on data length and pdf bin size.

Run _dips_ with:

`dips.py [-h] [-V] [-b BINS] [-t0 ORIGIN] [-P PERIOD] [-l LOGFILE] [-eta TOLERANCE] [-dxk DIFFERENCE] [-xi STEP_SIZE] [-af ATTENUATION] [--allow-upstep] [--cols COLS [COLS ...]] [--disable-mp] [--initial-pdf INITIAL_PDF] [--interim-prefix INTERIM_PREFIX] [--jitter JITTER] [--output-prefix OUTPUT_PREFIX] [--renormalize] [--save-interim SAVE_INTERIM] [--yonly] finput`

The arguments are summarized in the table below.

| Argument | Usage | Type | Default value |
|----------|-------|------|---------------|
| -h, --help | print out the help message and exit | n/a | n/a |
| -V, --version | print dips version and exit | n/a | n/a |
| -b BINS, --bins BINS | assign the number of synchronous pdf bins | int | 200 |
| -t0 ORIGIN, --origin ORIGIN | the zero-point of the timeseries | float | 0.0 |
| -P PERIOD, --period PERIOD | period of the synchronous signal | float | 1.0 |
| -l LOGFILE, --logfile LOGFILE | log file to send output to instead of screen | str | None |
| -eta TOLERANCE, --tolerance TOLERANCE | tolerance for convergence | float | 1e-8 |
| -dxk DIFFERENCE, --difference DIFFERENCE | finite difference size for computing slopes | float | 2e-5 |
| -xi STEP_SIZE, --step-size STEP_SIZE | initial down-step multiplier | float | 1e-3 |
| -af ATTENUATION, --attenuation ATTENUATION | attenuation factor for xi | float | 0.9 |
| --allow-upstep | allow step size to increase during convergence | bool | False |
| --cols COL1 COL2 \[COL3\] | a list of input columns to be parsed, starting from 0 | list of ints | 0 1 |
| --disable-mp | disable multiprocessing (force serial computation) | bool | False |
| --initial-pdf | choice of pdf initialization ('flat', 'mean', 'median', 'random', or external filename) | str | 'median' |
| --interim-prefix | filename prefix for interim results | str | finput |
| --output_prefix PREFIX | filename prefix for saving results (PREFIX.signal, .trend, .ranges) | str | finput |
| --renormalize | force pdf normalization to 1 after every iteration | bool | False |
| --save-interim STEP | save intering solutions every STEP iterations | int | 0 |
| --yonly | use only y-distance instead of full euclidian distance | bool | False |

Distributed with _dips_ (in the tarball's `examples` directory) are three example input files, `synthetic.data`, `kic3953981_sap.data` and `kic3547874_sap.data`.

To run _dips_ on synthetic data (see [http://keplerEBs.villanova.edu/includes/DPS/dps_synthetic.html](here) how the data were created) by using 33 bins, per-bin means as the initial pdf, and with serial calculation (disabling multiprocessing), issue:

```bash
dips synthetic.data -b 33 -P 0.91 --initial-pdf mean --disable-mp
```

To run _dips_ on an eclipsing binary [KIC 3953981](http://keplerEBs.villanova.edu/overview/?k=3953981), using 101 bins, allowing the step size to increase, using per-bin data median as the initial pdf, renormalizing the pdf after each iteration, using only y-direction length and saving every 10th iteration, issue:

```bash
dips kic3953981_sap.data -b 101 -t0 54953.82253243 -P 0.49201716 --allow-upstep --initial-pdf median --save-interim 10 --interim-prefix eb --renormalize --yonly
```

Finally, to run _dips_ on a heartbeat star [KIC 3547874](http://keplerEBs.villanova.edu/overview/?k=3547874), using 200 bins, starting with a flat pdf, computing total length in the y-direction only, renormalizing the synchronous pdf to 1.0 after each iteration, and allowing the step size to increase, issue:

```bash
dips kic3547874_sap.data --cols 0 2 -t0 54989.4209 -P 19.6921722 -b 200 --yonly --initial-pdf flat --renormalize --allow-upstep
```

These examples should provide a basic idea of how to invoke _dips_.
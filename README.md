# DetrendIng Periodic Signals (dips)

_dips_ is a simple algorithm for detrending timeseries of strictly periodic signals. It does not assume any functional form for the signal or the background or the noise; it disentangles the strictly periodic component from everything else. We use it in astronomy for detrending _Kepler_, _K2_ and _TESS_ timeseries of periodic variable stars, eclipsing binary stars, exoplanets etc. The algorithm is described in detail in Prsa et al. (2019), PASP, in review -- the reference will be updated shortly.

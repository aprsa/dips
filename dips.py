#!/usr/bin/python3

import os
import sys
import numpy as np
import argparse
#import matplotlib as mpl
#import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from scipy.stats import binned_statistic as hstats
import multiprocessing as mp

def fold(t, t0, P):
    return ((t-t0) % P) / P

def unfold(t, t0, P, ranges, pdf):
    return interp(0.5*(ranges[1:]+ranges[:-1]), pdf, k=3)(fold(t, t0, P))

def length(t, y, yonly=False):
    if yonly:
        return np.abs(y[1:]-y[:-1]).sum()
    else:
        return ( ( (y[1:]-y[:-1])**2 + (t[1:]-t[:-1])**2 )**0.5 ).sum()

def synclength(pdf, yonly=False):
    if yonly:
        return np.abs(pdf[1:]-pdf[:-1]).sum()
    else:
        return ( ( (pdf[1:]-pdf[:-1])**2 + (1./len(pdf))**2 )**0.5 ).sum()

def slope(k):
    # we assume POSIX compliance here so that main module memory is shared.

    x = pdf.copy()
    x[k] -= args.difference/2
    y = O - unfold(t, t0, P, ranges, x)
    l1 = length(t, y, args.yonly)
    x[k] += args.difference
    y = O - unfold(t, t0, P, ranges, x)
    l2 = length(t, y, args.yonly)
    if args.jitter == 0:
        return (l2-l1)/args.difference
    else:
        return np.random.normal((l2-l1)/args.difference, args.jitter*np.abs((l2-l1)/args.difference))

#mpl.rcParams['font.size'] = 18

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('finput',                  type=str,            help='input file containing time, flux and optional flux error')
    parser.add_argument('-b',   '--bins',          type=int,            help='number of synchronous pdf bins', default=200)
    parser.add_argument('-t0',  '--origin',        type=float,          help='the zero-point of the time-series', default=0.0)
    parser.add_argument('-P',   '--period',        type=float,          help='period of the synchronous signal', default=1.0)
    parser.add_argument('-l',   '--logfile',       type=str,            help='log file to send output to instead of screen', default=None)
    parser.add_argument('-eta', '--tolerance',     type=float,          help='tolerance for convergence', default=1e-8)
    parser.add_argument('-dxk', '--difference',    type=float,          help='finite difference size', default=2e-5)
    parser.add_argument('-xi',  '--step-size',     type=float,          help='initial down-step multiplier', default=1e-3)
    parser.add_argument('-af',  '--attenuation',   type=float,          help='attenuation factor for xi', default=0.9)
    parser.add_argument(        '--allow-upstep',  action='store_true', help='allow step size to increase during convergence', default=False)
    parser.add_argument(        '--cols',          type=int, nargs='+', help='a list of input columns to be parsed, starting from 0', default=[0, 1])
    parser.add_argument(        '--disable-mp',    action='store_true', help='disable multiprocessing (force serial computation)', default=False)
    parser.add_argument(        '--initial-pdf',   type=str,            help='choice of pdf initialization [\'flat\', \'mean\', \'median\', \'random\', or external filename]', default='median')
    parser.add_argument(        '--jitter',        type=float,          help='add jitter to the computed gradients', default=0.0)
    parser.add_argument(        '--output-prefix', type=str,            help='filename prefix for saving results', default=None)
    parser.add_argument(        '--renormalize',   action='store_true', help='force pdf normalization to 1 after every iteration', default=False)
    parser.add_argument(        '--save-interim',  type=int,            help='save interim solutions every N iterations', default=0)
    parser.add_argument(        '--yonly',         action='store_true', help='use only y-distance instead of full euclidian distance', default=False)

    args = parser.parse_args()

    # these are used a lot so assign local variables to them
    bins = args.bins
    xi = args.step_size
    t0 = args.origin
    P = args.period

    if args.logfile is not None:
        log = open(args.logfile, 'w')
    else:
        log = sys.stdout

    log.write('# Issued command:\n')
    log.write('#   %s\n# \n' % ' '.join(sys.argv))

    if os.name != 'posix':
        # as of right now this code won't run on non-posix systems.
        sys.stderr('dips does not run on a non-POSIX system; if you would like\nto run dips on your architecture, please open a ticket on https://github.com/aprsa/dips.\n')
        exit()

    if len(args.cols) == 2:
        t, O = np.loadtxt(args.finput, usecols=args.cols, unpack=True)
    elif len(args.cols) == 3:
        t, O, sO = np.loadtxt(args.finput, usecols=args.cols, unpack=True)
    else:
        raise argparse.ArgumentTypeError('only 2 or 3 columns can be passed to the --cols parameter.')

    log.write('# input data: %d rows, %d columns read in from %s\n# \n' % (len(t), len(args.cols), args.finput))

    ranges = np.linspace(0, 1, bins+1)

    if args.initial_pdf == 'flat':
        pdf = np.ones(bins)
    elif args.initial_pdf == 'mean':
        pdf = hstats(x=fold(t, t0, P), values=O, statistic='mean', bins=bins, range=(0, 1))[0]
    elif args.initial_pdf == 'median':
        pdf = hstats(x=fold(t, t0, P), values=O, statistic='median', bins=bins, range=(0, 1))[0]
    elif args.initial_pdf == 'random':
        pdf = np.random.normal(np.mean(O), np.std(O), bins)
    else:
        log.write('# initial pdf source: %s\n' % args.initial_pdf)
        pdf = np.loadtxt(args.initial_pdf, usecols=(1,))
        if len(pdf) != bins:
            log.write('#   rebinning the input pdf from %d to %d\n' % (len(pdf), bins))
            r = np.linspace(0, 1, len(pdf)+1)
            pdf = np.interp((ranges[1:]+ranges[:-1])/2, (r[1:]+r[:-1])/2, pdf)

    # plt.figure(figsize=(16,6))
    # plt.ylim(0.9, 1.01)
    # plt.xlabel('Phase')
    # plt.ylabel('Normalized flux')
    # plt.plot(fold(t, t0, P), O, 'b.')
    # plt.bar(0.5*(ranges[:-1]+ranges[1:]), pdf, width=1./bins, color='yellow', edgecolor='black', zorder=10, alpha=0.4)
    # plt.show()

    log.write('# number of requested pdf bins: %d\n' % bins)

    nelems_per_bin, _ = np.histogram(fold(t, t0, P), bins=bins)
    log.write('# number of observations per bin:\n')
    log.write('#   min: %d   max: %d   mean: %d\n# \n' % (nelems_per_bin.min(), nelems_per_bin.max(), nelems_per_bin.mean()))

    nprocs = 1 if args.disable_mp else mp.cpu_count()
    log.write('# dips running on %d %s (multiprocessing %s)\n# \n' % (nprocs, 'core' if nprocs == 1 else 'cores', 'off' if args.disable_mp else 'on'))

    Y = O - unfold(t, t0, P, ranges, pdf)
    log.write('# original timeseries length:  %f\n' % length(t, O, args.yonly))
    log.write('# initial asynchronous length: %f\n# \n' % length(t, Y, args.yonly))

    log.write('# computational parameters:\n')
    log.write('#   tolerance (tol):  %6.2e\n' % args.tolerance)
    log.write('#   difference (dxk): %6.2e\n' % args.difference)
    log.write('#   step size (xi):   %6.2e\n' % args.step_size)
    log.write('#   attenuation (af): %6.2e\n' % args.attenuation)
    log.write('#   yonly:            %s\n'    % args.yonly)
    log.write('#   slope jitter:     %2.2f\n' % args.jitter)
    log.write('#   renormalize:      %s\n'    % args.renormalize)
    log.write('# \n')

    i = 0
    l1 = length(t, O - unfold(t, t0, P, ranges, pdf), yonly=args.yonly)
    if args.disable_mp:
        slopes = np.array([slope(k) for k in range(bins)])
    else:
        with mp.Pool() as pool:
            slopes = np.array(pool.map(slope, range(bins)))
            # slopes = np.array(pool.starmap(slope, zip([mp_ranges]*bins, [pdf]*bins, range(bins), [args.difference]*bins, [mp_t]*bins, [mp_O]*bins, [mp_t0]*bins, [mp_P]*bins, [args.yonly]*bins, [args.jitter]*bins)))
    mean_slope = np.abs(slopes).mean()

    log.write('# %3s %14s %12s %14s %14s %14s\n' % ('it', 'async_length', 'sync_length', 'difference', 'step_size', 'mean_slope'))
    while xi*mean_slope > args.tolerance:
        l0 = l1
        while True:
            steps = -xi*slopes
            l1 = length(t, O - unfold(t, t0, P, ranges, pdf+steps), yonly=args.yonly)
            if l1 > l0:
                xi *= args.attenuation
            else:
                pdf += steps
                if args.renormalize:
                    pdf /= pdf.mean()
                if args.allow_upstep:
                    xi /= args.attenuation
                break

        i += 1
        log.write('%5d %14.8f %12.8f %14.8e %14.8e %14.8e\n' % (i, l1, synclength(pdf), l0-l1, xi, mean_slope))

        if args.save_interim > 0 and i % args.save_interim == 0:
            np.savetxt('%s.%04d.ranges' % (args.finput, i), ranges)
            np.savetxt('%s.%04d.signal' % (args.finput, i), np.vstack((0.5*(ranges[:-1]+ranges[1:]), pdf)).T)
            np.savetxt('%s.%04d.trend'  % (args.finput, i), np.vstack((t, O-unfold(t, t0, P, ranges, pdf))).T)

        if args.disable_mp:
            slopes = np.array([slope(k) for k in range(bins)])
        else:
            with mp.Pool() as pool:
                slopes = np.array(pool.map(slope, range(bins)))

        mean_slope = np.abs(slopes).mean()

    # mpl.rcParams['font.size'] = 24

    # plt.figure(figsize=(16,15))
    # plt.axes([0.1, 0.65, 0.8, 0.25])
    # plt.ylabel(r'$O(t)$')
    # plt.plot(t, O, 'b-', label='observed curve, length=%4.4f' % (length(t, O)))
    # plt.legend(loc='lower left')
    # plt.title('KIC 2445134')

    # plt.axes([0.1, 0.4, 0.8, 0.25])
    # plt.xlabel(r'$t$')
    # plt.ylabel(r'$O(t)-X(t)$')
    # plt.plot(t, O-unfold(t, t0, P, ranges, pdf), 'r-', label='asynchronous part, length=%4.4f' % (l1))
    # plt.legend(loc='lower left')

    # plt.axes([0.1, 0.08, 0.8, 0.25])
    # plt.xlabel(r'$\Phi$')
    # plt.ylabel(r'$X(t)$')
    # plt.plot(0.5*(ranges[:-1]+ranges[1:]), pdf, 'r-', lw=2, label='synchronous pdf')
    # plt.legend(loc='upper left')

    # plt.savefig('dps_2445134.pdf')
    # plt.show()

    prefix = args.finput if args.output_prefix is None else args.output_prefix

    np.savetxt('%s.ranges' % (prefix), ranges)
    np.savetxt('%s.signal' % (prefix), np.vstack((0.5*(ranges[:-1]+ranges[1:]), pdf)).T)
    np.savetxt('%s.trend' % (prefix), np.vstack((t, O-unfold(t, t0, P, ranges, pdf))).T)

    if args.logfile is not None:
        log.close()

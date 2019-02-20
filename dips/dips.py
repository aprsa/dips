#!/usr/bin/python3

import os
import sys
import numpy as np
import argparse
# import matplotlib as mpl
# import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from scipy.stats import binned_statistic as hstats
import multiprocessing as mp

class Dips:
    def __init__(self, args):
        self.args = args
        
        # this parameter changes so expose it to the class explicitly:
        self.xi = args['step_size']

        # can't serialize a file stream for multiprocessing, so don't do self.log.
        if args['logfile'] is not None:
            log = open(args['logfile'], 'w')
        else:
            log = sys.stdout

        log.write('# Issued command:\n')
        log.write('#   %s\n# \n' % ' '.join(sys.argv))

        self.data = np.loadtxt(args['finput'], usecols=args['cols'])
        log.write('# input data: %d rows, %d columns read in from %s\n' % (self.data.shape[0], self.data.shape[1], args['finput']))

        self.ranges = np.linspace(0, 1, args['bins']+1)
        self.phases = self.fold(self.data[:,0], args['origin'], args['period'])

        log.write('# initial pdf source: %s\n' % args['initial_pdf'])
        if args['initial_pdf'] == 'flat':
            self.pdf = np.ones(args['bins'])
        elif args['initial_pdf'] == 'mean':
            self.pdf = hstats(x=self.phases, values=self.data[:,1], statistic='mean', bins=args['bins'], range=(0, 1))[0]
        elif args['initial_pdf'] == 'median':
            self.pdf = hstats(x=self.phases, values=self.data[:,1], statistic='median', bins=args['bins'], range=(0, 1))[0]
        elif args['initial_pdf'] == 'random':
            means = hstats(x=self.phases, values=self.data[:,1], statistic='mean', bins=args['bins'], range=(0, 1))[0]
            stds  = hstats(x=self.phases, values=self.data[:,1], statistic='std', bins=args['bins'], range=(0, 1))[0]
            self.pdf = np.random.normal(means, stds)
        else:
            self.pdf = np.loadtxt(args['initial_pdf'], usecols=(1,))
            if len(self.pdf) != args['bins']:
                log.write('#   rebinning the input pdf from %d to %d\n' % (len(self.pdf), args['bins']))
                r = np.linspace(0, 1, len(self.pdf)+1)
                self.pdf = np.interp((self.ranges[1:]+self.ranges[:-1])/2, (r[1:]+r[:-1])/2, self.pdf)

        log.write('# number of requested pdf bins: %d\n' % args['bins'])

        nelems_per_bin, _ = np.histogram(self.phases, bins=args['bins'])
        log.write('# number of observations per bin:\n')
        log.write('#   min: %d   max: %d   mean: %d\n# \n' % (nelems_per_bin.min(), nelems_per_bin.max(), nelems_per_bin.mean()))

        nprocs = 1 if args['disable_mp'] else mp.cpu_count()
        log.write('# dips running on %d %s (multiprocessing %s)\n# \n' % (nprocs, 'core' if nprocs == 1 else 'cores', 'off' if args['disable_mp'] else 'on'))

        # mpl.rcParams['font.size'] = 24
        # plt.figure(figsize=(16,6))
        # plt.ylim(0.9, 1.01)
        # plt.xlabel('Phase')
        # plt.ylabel('Normalized flux')
        # plt.plot(self.phases, self.data[:,1], 'b.')
        # plt.bar(0.5*(self.ranges[:-1]+self.ranges[1:]), self.pdf, width=1./args['bins'], color='yellow', edgecolor='black', zorder=10, alpha=0.4)
        # plt.show()

        Y = self.data[:,1] - self.unfold(self.pdf)
        log.write('# original timeseries length:  %f\n' % self.length(self.data[:,0], self.data[:,1]))
        log.write('# initial asynchronous length: %f\n# \n' % self.length(self.data[:,0], Y))

        log.write('# computational parameters:\n')
        log.write('#   tolerance (tol):  %6.2e\n' % args['tolerance'])
        log.write('#   difference (dxk): %6.2e\n' % args['difference'])
        log.write('#   step size (xi):   %6.2e\n' % args['step_size'])
        log.write('#   attenuation (af): %6.2e\n' % args['attenuation'])
        log.write('#   up-step allowed:  %s\n'    % args['allow_upstep'])
        log.write('#   yonly:            %s\n'    % args['yonly'])
        log.write('#   slope jitter:     %2.2f\n' % args['jitter'])
        log.write('#   renormalize:      %s\n'    % args['renormalize'])
        log.write('# \n')

        if self.args['logfile'] is not None:
            log.close()

    def run(self):
        # can't serialize a file stream for multiprocessing, so don't do self.log.
        if self.args['logfile'] is not None:
            log = open(self.args['logfile'], 'a')
        else:
            log = sys.stdout

        if self.args['save_interim'] > 0:
            interim_prefix = self.args['finput'] if self.args['interim_prefix'] is None else self.args['interim_prefix']

        log.write('# %3s %14s %12s %14s %14s %14s\n' % ('it', 'async_length', 'sync_length', 'difference', 'step_size', 'mean_slope'))

        # starting iterations
        i = 0
        l1 = self.length(self.data[:,0], self.data[:,1] - self.unfold(self.pdf))
        if self.args['disable_mp']:
            slopes = np.array([self.slope(k) for k in range(self.args['bins'])])
        else:
            with mp.Pool() as pool:
                slopes = np.array(pool.map(self.slope, range(self.args['bins'])))
        mean_slope = np.abs(slopes).mean()

        while self.xi*mean_slope > self.args['tolerance']:
            l0 = l1
            safety_counter = 0
            while True:
                steps = -self.xi*slopes
                l1 = self.length(self.data[:,0], self.data[:,1] - self.unfold(self.pdf+steps))
                if l1 > l0 and safety_counter < 1000:
                    self.xi *= self.args['attenuation']
                    safety_counter += 1
                else:
                    self.pdf += steps
                    if self.args['renormalize']:
                        self.pdf /= self.pdf.mean()
                    if self.args['allow_upstep']:
                        self.xi /= self.args['attenuation']
                    break

            if self.args['save_interim'] > 0 and i % self.args['save_interim'] == 0:
                np.savetxt('%s.%05d.ranges' % (interim_prefix, i), self.ranges)
                np.savetxt('%s.%05d.signal' % (interim_prefix, i), np.vstack((0.5*(self.ranges[:-1]+self.ranges[1:]), self.pdf)).T)
                np.savetxt('%s.%05d.trend'  % (interim_prefix, i), np.vstack((self.data[:,0], self.data[:,1]-self.unfold(self.pdf))).T)

            i += 1
            log.write('%5d %14.8f %12.8f %14.8e %14.8e %14.8e\n' % (i, l1, self.synclength(self.pdf), l0-l1, self.xi, mean_slope))

            if self.args['disable_mp']:
                slopes = np.array([self.slope(k) for k in range(self.args['bins'])])
            else:
                with mp.Pool() as pool:
                    slopes = np.array(pool.map(self.slope, range(self.args['bins'])))

            mean_slope = np.abs(slopes).mean()

        prefix = self.args['finput'] if self.args['output_prefix'] is None else self.args['output_prefix']

        np.savetxt('%s.ranges' % (prefix), self.ranges)
        np.savetxt('%s.signal' % (prefix), np.vstack((0.5*(self.ranges[:-1]+self.ranges[1:]), self.pdf)).T)
        np.savetxt('%s.trend' % (prefix), np.vstack((self.data[:,0], self.data[:,1]-self.unfold(self.pdf))).T)

        if self.args['logfile'] is not None:
            log.close()

    def fold(self, t, t0, P):
        return ((t-t0) % P) / P

    def unfold(self, pdf):
        return interp(0.5*(self.ranges[1:]+self.ranges[:-1]), pdf, k=3)(self.phases)

    def length(self, t, y):
        if self.args['yonly']:
            return np.abs(y[1:]-y[:-1]).sum()
        else:
            return ( ( (y[1:]-y[:-1])**2 + (t[1:]-t[:-1])**2 )**0.5 ).sum()

    def synclength(self, pdf):
        if self.args['yonly']:
            return np.abs(pdf[1:]-pdf[:-1]).sum()
        else:
            return ( ( (pdf[1:]-pdf[:-1])**2 + (1./len(pdf))**2 )**0.5 ).sum()

    def slope(self, k):
        x = self.pdf.copy()
        x[k] -= self.args['difference']/2
        l1 = self.length(self.data[:,0], self.data[:,1] - self.unfold(x))
        x[k] += self.args['difference']
        l2 = self.length(self.data[:,0], self.data[:,1] - self.unfold(x))
        if self.args['jitter'] == 0:
            return (l2-l1)/self.args['difference']
        else:
            return np.random.normal((l2-l1)/self.args['difference'], self.args['jitter']*np.abs((l2-l1)/self.args['difference']))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('finput',                   type=str,            help='input file containing time, flux and optional flux error')
    parser.add_argument('-b',   '--bins',           type=int,            help='number of synchronous pdf bins', default=200)
    parser.add_argument('-t0',  '--origin',         type=float,          help='the zero-point of the time-series', default=0.0)
    parser.add_argument('-P',   '--period',         type=float,          help='period of the synchronous signal', default=1.0)
    parser.add_argument('-l',   '--logfile',        type=str,            help='log file to send output to instead of screen', default=None)
    parser.add_argument('-eta', '--tolerance',      type=float,          help='tolerance for convergence', default=1e-8)
    parser.add_argument('-dxk', '--difference',     type=float,          help='finite difference size', default=2e-5)
    parser.add_argument('-xi',  '--step-size',      type=float,          help='initial down-step multiplier', default=1e-3)
    parser.add_argument('-af',  '--attenuation',    type=float,          help='attenuation factor for xi', default=0.9)
    parser.add_argument(        '--allow-upstep',   action='store_true', help='allow step size to increase during convergence', default=False)
    parser.add_argument(        '--cols',           type=int, nargs='+', help='a list of input columns to be parsed, starting from 0', default=[0, 1])
    parser.add_argument(        '--disable-mp',     action='store_true', help='disable multiprocessing (force serial computation)', default=False)
    parser.add_argument(        '--initial-pdf',    type=str,            help='choice of pdf initialization [\'flat\', \'mean\', \'median\', \'random\', or external filename]', default='median')
    parser.add_argument(        '--interim-prefix', type=str,            help='filename prefix for interim results', default=None)
    parser.add_argument(        '--jitter',         type=float,          help='add jitter to the computed gradients', default=0.0)
    parser.add_argument(        '--output-prefix',  type=str,            help='filename prefix for saving results', default=None)
    parser.add_argument(        '--renormalize',    action='store_true', help='force pdf normalization to 1 after every iteration', default=False)
    parser.add_argument(        '--save-interim',   type=int,            help='save interim solutions every N iterations', default=0)
    parser.add_argument(        '--yonly',          action='store_true', help='use only y-distance instead of full euclidian distance', default=False)

    cmdargs = parser.parse_args()
    if len(cmdargs.cols) < 2 or len(cmdargs.cols) > 3:
        raise argparse.ArgumentTypeError('only 2 or 3 columns can be passed to the --cols parameter.')

    dips = Dips(vars(cmdargs))
    dips.run()








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



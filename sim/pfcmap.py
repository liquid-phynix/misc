import numpy as n
from numpy import array, linspace, sin, cos, real, imag, prod, pi, sqrt, cosh
from numpy.fft import *
from pylab import *
import matplotlib.cm as cm
import fftw3
from fftw3.lib import lib
from sys import stdout

class P2Map(object):
    def __init__(self, nsx=10.0, nsy=10.0 ,nx=200):
        nsy /= 2
        self.sigma = 4.0 * pi / sqrt(3.0)
        self.lx = self.sigma * nsx
        self.ly = self.sigma * nsy * sqrt(3.0)
        self.nx = nx
        self.ny = int(nx * self.ly / self.lx)
        vecx = linspace(0,self.lx,num=self.nx,endpoint=False)
        vecy = linspace(0,self.ly,num=self.ny,endpoint=False)
        self.X,self.Y = meshgrid(vecx, vecy)

        nthreads = 1
        # flags = ['measure']
        f = fftw3.fftw_flags
        
        self.y   = fftw3.create_aligned_array((self.ny, self.nx), dtype=typeDict['complex'])
        self.yp  = fftw3.create_aligned_array((self.ny, self.nx), dtype=typeDict['complex'])
        self.y3  = fftw3.create_aligned_array((self.ny, self.nx), dtype=typeDict['complex'])
        self.y3p = fftw3.create_aligned_array((self.ny, self.nx), dtype=typeDict['complex'])

        self.y_to_yp = lib.fftw_plan_dft_2d(self.y.shape[0], self.y.shape[1], self.y, self.yp, -1, f['measure'] | f['preserve input'])
        self.yp_to_y = lib.fftw_plan_dft_2d(self.y.shape[0], self.y.shape[1], self.yp, self.y, 1, f['measure'] | f['preserve input'])
        self.y3_to_y3p = lib.fftw_plan_dft_2d(self.y3.shape[0], self.y3.shape[1], self.y3, self.y3p, -1, f['measure'] | f['preserve input'])
        self.y3p_to_y3 = lib.fftw_plan_dft_2d(self.y3.shape[0], self.y3.shape[1], self.y3p, self.y3, 1, f['measure'] | f['preserve input'])

        xlen, ylen = self.nx, self.ny
        kx = array([x for x in range(xlen//2+1)] + [x for x in range(-(xlen-xlen//2-1), 0)])
        ky = array([x for x in range(ylen//2+1)] + [x for x in range(-(ylen-ylen//2-1), 0)])

        self.lap = zeros_like(self.yp)
        lx, ly = self.lx, self.ly
        for iy in xrange(self.lap.shape[0]):
            for ix in xrange(self.lap.shape[1]):
                self.lap[iy,ix] = (-2*pi*1j*kx[ix]/lx)**2.0+(-2*pi*1j*ky[iy]/ly)**2.0
        self.dop_partial = self.lap * (1 + self.lap)**2.0
        self.dop = None
        
    def init(self, itype, amp, r, psi0):
        self.r = r
        self.psi0 = psi0
        if itype == 'hex':
            sqrt3p2 = sqrt(3.0) / 2.0
            self.y[:] = self.psi0 + amp * (2.0 * cos(self.X[:] * sqrt3p2) * cos(self.Y[:] / 2.0) + cos(self.Y[:])) / 3.0
        elif itype == 'ran':
            self.y[:] = psi0 + 0.001 * (rand(self.ny, self.nx)-0.5)[:]
        else:
            ValueError('not a proper initial condition')
            
        self.dop = self.dop_partial + self.lap * self.r
        lib.fftw_execute(self.y_to_yp)

        # init
        # if itype == 'ran':
        #     self.y[:] = self.psi0 + 0.001 * (rand(self.ny, self.nx)-0.5)[:]
        # elif itype == 'hex':
        #     self.y[:] = self.psi0 + (2.0 * cos(self.X[:] * sqrt3p2) * cos(self.Y[:] / 2.0) + cos(self.Y[:])) / 3.0
        # elif itype == 'hom':
        #     self.y[:] = self.psi0
        # elif itype == 'stripe':
        #     self.y[:] = self.psi0 + 0.5 * cos(self.X[:])
        # elif itype == 'stripe-exp':
        #     self.y[:] = self.psi0 + cos(self.X[:]) * sqrt(2.0 * (- self.r - 3.0 * self.psi0 ** 2.0) / 3.0) * exp(-0.0005 * (self.X[:] - self.lx / 2.0)**2.0)
        # else:
        #     raise ValueError('not a proper initial condition')

    def plot(self):
        # Plot 1
        subplot(211)
        imshow(real(self.y), interpolation='none', cmap=cm.spectral)
        # Plot 2
        subplot(212)
        plot(real(self.y[0,:]))
        pause(0.1)
                

    def plot_save(self):
        for axis in self.axes: axis.cla()
        self.axes[0].imshow(real(self.y), interpolation='none', cmap=cm.spectral)
        self.axes[1].plot(real(self.y.sum(axis = 0) / float(self.y.shape[0])))
        self.axes[2].plot(real(ifft2(self.yp * self.kmask).sum(axis=0)) / float(self.y.shape[0]))
        for fig, prefix in zip(self.figs, ['1', '2', '3']):
            fig.savefig('movie/' + prefix + '.%05d.png' % self.frame_counter)
        pause(0.1)

    def step(self,dt):
        self.y3[:] = self.y[:]**3.0
        lib.fftw_execute(self.y3_to_y3p)
        self.yp[:] = (self.yp[:] + dt * self.lap[:] * self.y3p[:]) / (1.0 - dt * self.dop[:])
        lib.fftw_execute(self.yp_to_y)
        self.y[:] /= prod(self.y.shape)
        
        # self.frame_counter += 1
        # if self.frame_counter % 10 == 0:
        #     self.saved_curves += [real(self.y[0,:]).copy()]
        
    def run(self,dt=0.01,steps=500,freq=100):
        figure()
        self.plot()
        raw_input()
        while steps>0:
            stdout.flush()
            self.step(dt)
            if steps%freq==0:
                print steps,
                self.plot()
            steps-=1

    def run_save(self,dt=0.01,steps=500):
        self.figs = [figure(figsize=(15,4)) for _ in range(3)]
        self.axes = [f.add_subplot(111) for f in self.figs]
        self.plot_save()
        while steps > 0:
            print steps,
            stdout.flush()
            self.step(dt)
            if steps%100==0:
                self.plot_save()
            steps-=1

    def run_aut(self, r, psi0, dt = 0.01):
        steps = 0
        y = real(self.y.copy())
        while True:
            self.step(dt)
            steps += 1
            if steps % 10 == 0:
                diff = (abs(y - real(self.y))).sum() / prod(y.shape)
                if diff < 1e-5: break
                else: y = real(self.y.copy())
        return r, psi0, y.min(), y.max(), y.sum() / prod(y.shape), steps

    def map_space(self, from_r, to_r, from_psi0, to_psi0, divs_r, divs_psi0):
        self.rvec = linspace(from_r, to_r, divs_r)
        self.psi0vec = linspace(from_psi0, to_psi0, divs_psi0)
        self.amplitudes_min = zeros((divs_r, divs_psi0))
        self.amplitudes_max = zeros((divs_r, divs_psi0))
        self.densities = zeros((divs_r, divs_psi0))
        self.iters = zeros((divs_r, divs_psi0))

        for i, r in enumerate(self.rvec):
            for j, psi0 in enumerate(self.psi0vec):
                if j == 0:
                    self.init('hex', 0.5, r, psi0)
                else:
                    self.y[:] = self.y[:] - self.psi0vec[j-1] + psi0
                    lib.fftw_execute(self.y_to_yp)
                    
                self.amplitudes_min[i, j], self.amplitudes_max[i, j], self.densities[i, j], self.iters[i, j] = self.run_aut()

                clf()
                imshow(self.amplitudes_max - self.densities, interpolation = 'none')
                colorbar()
                pause(0.1)
                
    def map_space_boundary(self, from_r, to_r, psi0_neg, psi0_pos, divs_r, divs_psi0):
        self.result = []
        
        for r in linspace(from_r, to_r, divs_r):
            psi0_boundary = - sqrt(abs(r) / 2.35)
            prev_psi0 = None
            print 'r = {0}, psi0_min = {1}, psi0_max = {2}\n'.format(r, psi0_boundary + psi0_pos, psi0_boundary - psi0_neg)
            stdout.flush()
            for psi0 in linspace(psi0_boundary + psi0_pos, psi0_boundary - psi0_neg, divs_psi0):
                if prev_psi0 == None:
                    self.init('hex', 0.5, r, psi0)
                else:
                    self.y[:] = self.y[:] - prev_psi0 + psi0
                    lib.fftw_execute(self.y_to_yp)
                    
                prev_psi0 = psi0

                self.result += [self.run_aut(r, psi0)]
                            
    def save_results(prefix):
        savetxt(prefix + '_amplitudes_min',self.amplitudes_min)
        savetxt(prefix + '_amplitudes_max',self.amplitudes_max)
        savetxt(prefix + '_densities',self.densities)
        savetxt(prefix + '_iters',self.iters)
        self.amplitudes_min.tofile(prefix + '_amplitudes_min_bin')
        self.amplitudes_max.tofile(prefix + '_amplitudes_max_bin')
        self.densities.tofile(prefix + '_densities_bin')

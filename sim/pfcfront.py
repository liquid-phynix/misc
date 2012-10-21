import numpy as n
from numpy import array, linspace, sin, cos, real, imag, prod, pi, sqrt, cosh
from numpy.fft import *
from pylab import *
import matplotlib.cm as cm
import fftw3
from fftw3.lib import lib
from sys import stdout

class P2Front(object):
    def __init__(self,r=-0.2, psi0=-0.5, nsx=10.0, nsy=10.0 ,nx=200, init = 'hex'):
        self.frame_counter = 0
        nsy /= 2
        self.sigma = 4.0 * pi / sqrt(3.0)
        self.r,self.psi0 = float(r),float(psi0)
        self.lx = (2.0 * pi * nsx) if (init == 'stripe') else (self.sigma * nsx)
        self.ly = self.sigma * nsy * sqrt(3.0)
        self.nx = nx
        self.ny = int(nx * self.ly / self.lx)
        vecx = linspace(0,self.lx,num=self.nx,endpoint=False)
        self.xvec = vecx
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

        sqrt3p2 = sqrt(3.0) / 2.0
        # init
        if init == 'ran':
            self.y[:] = self.psi0 + 0.001 * (rand(self.ny, self.nx)-0.5)[:]
        elif init == 'hex':
            self.y[:] = self.psi0 + (2.0 * cos(self.X[:] * sqrt3p2) * cos(self.Y[:] / 2.0) + cos(self.Y[:])) / 3.0
        elif init == 'hom':
            self.y[:] = self.psi0
        elif init == 'stripe':
            self.y[:] = self.psi0 + 0.5 * cos(self.X[:])
        elif init == 'stripe-exp':
            self.y[:] = self.psi0 + cos(self.X[:]) * sqrt(2.0 * (- self.r - 3.0 * self.psi0 ** 2.0) / 3.0) * exp(-0.0005 * (self.X[:] - self.lx / 2.0)**2.0)
        else:
            raise ValueError('not a proper initial condition')

        # envelope = 0.5 * (tanh(- 0.2 * (self.X - self.l / 2.0 - self.l / 15.0)) + 1.0) * 0.5 * (tanh(0.2 * (self.X - self.l / 2.0 + self.l / 15.0)) + 1.0)
        # self.y[:] = self.psi0 * (1.0 - envelope[:]) + envelope[:] * self.y[:]
        
        # for ix in xrange(self.y.shape[0]):
        #     for iy in xrange(self.y.shape[1]):
        #         if (ix-self.n/2)**2 + (iy-self.n/2)**2 > (self.n/8)**2:
        #             self.y[ix,iy] = self.psi0
 

        xlen, ylen = self.nx, self.ny
        kx = array([x for x in range(xlen//2+1)] + [x for x in range(-(xlen-xlen//2-1), 0)])
        ky = array([x for x in range(ylen//2+1)] + [x for x in range(-(ylen-ylen//2-1), 0)])

        self.lap = zeros_like(self.yp)
        lx, ly = self.lx, self.ly
        for iy in xrange(self.lap.shape[0]):
            for ix in xrange(self.lap.shape[1]):
                self.lap[iy,ix] = (-2*pi*1j*kx[ix]/lx)**2.0+(-2*pi*1j*ky[iy]/ly)**2.0
        self.dop = self.lap * (self.r + (1 + self.lap)**2.0)

        lib.fftw_execute(self.y_to_yp)
        
        #self.x = linspace(0, self.l, self.n, endpoint=False)
        #self.y = sin(2*2*pi/self.l*self.x)+sin(5*2*pi/self.l*self.x)+sin(3*2*pi/self.l*self.x)
        # tmp = rfft2(self.y)
        # mask = zeros_like(tmp)
        # mask[0:mask.shape[0]/2, 0:mask.shape[0]/2] = 1.0
        # self.y = irfft(tmp * mask)
        #self.y = self.psi0 + zeros_like(self.x)
        #        self.y = self.psi0 + zeros(self.n)#+ 0.000001 * randn(self.n)
        #        self.y = self.psi0+0.01*sin(2*pi*self.x)
        #        self.y = self.psi0+0.001*cos(2.0*pi*self.x)*exp(-2.0*(self.x-self.l/2.0)**2.0)
        # self.kmask=ones_like(self.yp)
        # self.kmask[:,:20]=0
        # self.kmask[:,-20:]=0
        # self.kmask[:20,:]=0
        # self.kmask[-20:,:]=0

        #        self.kmask = array([[exp(-0.08 * (k_x**2.0 + k_y**2.0)) for k_x in kx] for k_y in ky])
        self.kmask = array([[exp(-0.05 * (k_x**2.0 + (k_y ** 2.0 * self.lx / self.ly))) for k_x in kx] for k_y in ky])
        
        self.saved_curves = []


    def plot(self):
        # Plot 1
        subplot(311)
        imshow(real(self.y), interpolation='none', cmap=cm.spectral)
        # Plot 2
        subplot(312)
        #        imshow(abs(self.kmask))
        #        imshow(real(ifft2(self.yp * self.kmask)))
        plot(real(self.y[0,:]))
        #        plot(real(self.y.sum(axis = 0) / float(self.y.shape[0])))
        # Plot 3
        subplot(313)
        plot(real(ifft2(self.yp * self.kmask).sum(axis=0)) / float(self.y.shape[0]))
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
        self.frame_counter += 1
        if self.frame_counter % 10 == 0:
            self.saved_curves += [(dt * float(self.frame_counter),
                                   real(self.y[0,:]).copy())]
        
    def make_front(self, psiL):
        envelope = 0.5 * (tanh(- 0.2 * (self.X - self.lx / 2.0 - self.lx / 12.0 / 4.0)) + 1.0) * 0.5 * (tanh(0.2 * (self.X - self.lx / 2.0 + self.lx / 12.0 / 4.0)) + 1.0)
        enveloped = cosh(0.2*(-(1./2.-1./12./4.)*self.lx + self.X))**-2. * (1. - tanh(0.2*(-(1./2.+1./12./4.)*self.lx + self.X))) - cosh(0.2*(-(1./2.+1./12./4.)*self.lx + self.X))**-2. * (1. + tanh(0.2*(-(1./2.-1./12./4.)*self.lx + self.X)))

        #        envelope = amp(self.X, self.lx / 2.0 + self.lx / 12.0 / 4.0, 20) * amp(-self.X, - (self.lx / 2.0 - self.lx / 12.0 / 4.0), 20)

        #        envelope = cos(0.15730 * self.X)
        #        self.y[:] = psiL + envelope[:] * self.y[:]
        self.y[:] = psiL * (1.0 - envelope[:]) + envelope[:] * self.y[:]
        #        self.y[:] = self.y[:] + 0.05 * n.abs(enveloped)
        
        lib.fftw_execute(self.y_to_yp)

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

# def amp(x,x0,w):
#     if x > x0 + w:
#         return 0.
#     elif x > x0:
#         return cos(pi * (x - x0) / 2.0 / w)
#     else:
#         return 1.

def amp(x,x0,w):
    ret = n.ones_like(x)
    indices = x > x0
    ret[indices] = cos(pi * (x[indices] - x0) / 2.0 / w)
    ret[x > x0 + w] = 0.
    return ret

# from pfcfft2dfftw import *
# p=P2Front(n=200,l=100,psi0=-0.2)
# p.run()
# mencoder 'mf://3.*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o vid3.avi

# for i, frame in enumerate(data):
#      ...:     ax.cla()
#      ...:     ax.plot(frame[len(frame)//2:])
#      ...:     fig.savefig("movie2/%04d.png"%i)
# for i, frame in enumerate(data):
#      ...:     ax.cla()
#      ...:     ax.plot(frame[len(frame)//2:])
#      ...:     fname = "movie2/%04d.png"%i
#      ...:     fig.savefig(fname)
#      ...:     print 'saved ', fname
#      ...:     pause(0.1)

import numpy as n
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

def find_max(fun, start):
    dx = 0.1
    x = start
    fx = fun(x)
    while dx > 1e-5:
        fleft = fun(x - dx)
        fright = fun(x + dx)
        if fleft > fx:
            x = x - dx
            fx = fleft
        elif fright > fx:
            x = x + dx
            fx = fright
        else:
            dx = dx / 2.0
    return x

def loc_max(xvec, ar):
    indices = (ar > n.roll(ar,1)) & (ar > n.roll(ar,-1))
    return xvec[indices], ar[indices]

# def find_pos((xvec,ar)):
#     diff = ar-n.roll(ar,1)
#     maxpos = xvec[diff >= n.max(diff)][0]
#     funar = interp1d(xvec,diff,kind='quadratic')
#     return find_max(funar, maxpos)

def find_pos(xvec, ar):
    xvec, ar = loc_max(xvec, ar)
    indices = array(range(len(xvec)))
    inflexi = indices[(ar > 0.1 * ar.min()) & (ar < 0.9 * ar.max())][0]
    maxpos = xvec[inflexi]

    xvec = xvec[inflexi - 5 : inflexi + 5]
    ar = ar[inflexi - 5 : inflexi + 5]
    funar = interp1d(xvec, ar, kind = 'quadratic')
    xvec = linspace(xvec[0], xvec[-1], 10 * len(xvec))
    ar = array(map(funar, xvec))
    diff = ar - n.roll(ar, 1)
    funar = interp1d(xvec, diff, kind = 'quadratic')
    
    return find_max(funar, maxpos)


# map(lambda frame: find_pos(loc_max(frame)), data[11:350])

# import pickle
# with open("data.dat","w") as f:
#      pickle.dump(data, f)

# poses = []
# for i,frame in enumerate(data[12:]):
#     find_pos(loc_max(frame))
#     print 'done', i

#figure(); plot([pfcfront.find_pos(xvec, curve) for curve in p.saved_curves[200::10]])

from scipy.optimize import leastsq
def lin_fit(xv, yv):
    residuals = lambda ps, y, x: y - (ps[0] * x + ps[1])
    return leastsq(residuals, [1., 0.], args = (yv, xv))


# tvec = array([c[0] for c in p.saved_curves[200::100]])
# poses = array([pfcfront.find_pos(p.xvec, c[1]) for c in p.saved_curves[200::100]])
# pfcfront.lin_fit(tvec, poses)

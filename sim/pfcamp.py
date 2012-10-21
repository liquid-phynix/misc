import numpy         as     n
from   numpy         import array, linspace, sin, cos, real, imag, prod, pi, sqrt, cosh, dot, abs
from   numpy.fft     import *
from   pylab         import *
import matplotlib.cm as     cm
import fftw3
from   fftw3.lib     import lib
from   sys           import stdout

class P1Amp(object):

    def fftw_array(self, divs):
        return fftw3.create_aligned_array(divs, dtype=typeDict['complex'])

    def fftw_plan(self, src, dest, dir):
        #        nthreads = 1
        f = fftw3.fftw_flags
        return lib.fftw_plan_dft_1d(src.shape[0], src, dest, dir, f['measure'] | f['preserve input'])
    
    def __init__(self, eps = 0.2, psi0 = 0.25, lz = 10.0, nz = 200, theta = 0, init = 'ran', w = 10.0, falloff = 0.2, preamp = 1.0):
        self.frame_counter = 0

        self.eps, self.psi0, self.lz, self.nz = float(eps), float(psi0), float(lz), nz
        self.zvec = linspace(0, self.lz, num = self.nz, endpoint=False)
        self.kappa = self.eps - 3.0 * self.psi0 ** 2.0
        self.theta = theta
        print 'e-3p0^2=', self.kappa
        
        normal = array([sin(theta), cos(theta)])
        rotmat = array([[cos(theta), sin(theta)],
                        [-sin(theta), cos(theta)]])
        
        hex_base = array([[-0.5 * sqrt(3), -0.5],
                          [0, 1],
                          [0.5 * sqrt(3), -0.5]])
        self.etas = [dot(v, normal) for v in hex_base]

        rs = ('amp1', 'amp2', 'amp3', 'nlin')
        ks = ('amp1_f', 'amp2_f', 'amp3_f', 'nlin_f')
        for rspace, kspace in zip(rs, ks):
            r = self.__dict__[rspace] = self.fftw_array(self.nz)
            k = self.__dict__[kspace] = self.fftw_array(self.nz)
            self.__dict__[rspace + '_to_' + kspace] = self.fftw_plan(r, k, -1)
            self.__dict__[kspace + '_to_' + rspace] = self.fftw_plan(k, r, 1)

        # init
        if init == 'ran':
            self.amp1[:] = 0.01 * (rand(self.nz) - 0.5)[:] + 0.01j * (rand(self.nz) - 0.5)[:]
            self.amp2[:] = 0.01 * (rand(self.nz) - 0.5)[:] + 0.01j * (rand(self.nz) - 0.5)[:]
            self.amp3[:] = 0.01 * (rand(self.nz) - 0.5)[:] + 0.01j * (rand(self.nz) - 0.5)[:]
        elif init == 'exp2':
            self.amp1[:] = sqrt(2.0 * (self.eps - 3.0 * self.psi0 ** 2.0) / 3.0) * exp(-0.5 * (self.zvec[:] - self.lz / 2.0) ** 2.0)
        elif init == 'tanh':
            envelope = 0.5 * (tanh(- falloff * (self.zvec - self.lz / 2.0 - self.lz / w)) + 1.0) * \
              0.5 * (tanh(falloff * (self.zvec - self.lz / 2.0 + self.lz / w)) + 1.0)
              #            self.amp1[:] = - preamp * sqrt(2.0 * (self.eps - 3.0 * self.psi0 ** 2.0) / 3.0) * envelope[:]
              #            self.amp1[:] = - sign(self.kappa) * 0.1 * preamp * envelope[:]
            self.amp1[:] = - 0.1 * preamp * envelope[:]
            self.amp2[:] = self.amp1[:]
            self.amp3[:] = self.amp1[:]
        else:
            raise ValueError('not a proper initial condition')

        for r, k in zip(rs, ks)[:-1]:            
            lib.fftw_execute(self.__dict__[r + '_to_' + k])
            
        kz = array([k for k in range(self.nz // 2 + 1)] + [k for k in range(-(self.nz - self.nz // 2 - 1), 0)])
        der1 = - 2 * pi * 1j * kz / lz
        der2 = der1 * der1
        der3 = der2 * der1
        der4 = der3 * der1
        der5 = der4 * der1
        der6 = der5 * der1

        kappa = self.kappa
        # QDRG
        self.linops = \
          [ der6
            + 6.0j * eta * der5 \
            - (1.0 + 12.0 * eta ** 2.0) * der4 \
            - 4.0j * eta * (1.0 + 2.0 * eta ** 2.0) * der3 \
            + (4.0 * eta ** 2.0 - kappa) * der2 \
            - kappa * 2.0j * eta * der1 \
            + kappa \
            for eta in self.etas]
        
        # proto-RG
        # self.linops = \
        #   [ der6
        #     + 6.0j * eta * der5 \
        #     - (1.0 + 12.0 * eta ** 2.0) * der4 \
        #     - 4.0j * eta * (1.0 + 2.0 * eta ** 2.0) * der3 \
        #     + (4.0 * eta ** 2.0) * der2 \
        #     + kappa \
        #     for eta in self.etas]

        self.saved_curves = []

        xtoy = 16
        upperright = array([self.lz/xtoy/2, self.lz/2])
        lowerleft = -upperright
        self.X, self.Y = meshgrid(linspace(lowerleft[0], upperright[0], self.nz // xtoy), linspace(lowerleft[1], upperright[1], self.nz))
        self.Y = -self.Y
        self.X, self.Y = rotmat[0,0] * self.X + rotmat[0,1] * self.Y, rotmat[1,0] * self.X + rotmat[1,1] * self.Y
        
        for i, name in enumerate(('planew1', 'planew2', 'planew3')):
            pw = self.__dict__[name] = exp(1j * (hex_base[i][0] * self.X + hex_base[i][1] * self.Y))
            self.__dict__[name + '_conj'] = pw.conj()

    def step(self, dt):
        amps = (self.amp1, self.amp2, self.amp3)
        amps_f = (self.amp1_f, self.amp2_f, self.amp3_f)
        
        for i,j,l in ((0,1,2),(1,0,2),(2,0,1)):
            nlin = - 3.0 * amps[i] * (amps[i] * amps[i].conj() + 2.0 * amps[j] * amps[j].conj() + 2.0 * amps[l] * amps[l].conj()) - 6.0 * self.psi0 * amps[j].conj() * amps[l].conj()
            self.nlin[:] = nlin[:]
            lib.fftw_execute(self.nlin_to_nlin_f)
            
            amps_f[i][:] = (amps_f[i][:] + dt * self.nlin_f[:]) / (1.0 - dt * self.linops[i][:])
        
        lib.fftw_execute(self.amp1_f_to_amp1)
        lib.fftw_execute(self.amp2_f_to_amp2)
        lib.fftw_execute(self.amp3_f_to_amp3)

        self.amp1[:] /= float(self.nz)
        self.amp2[:] /= float(self.nz)
        self.amp3[:] /= float(self.nz)
        
        self.frame_counter += 1
        if self.frame_counter % 10 == 0:
            self.saved_curves += [(dt * float(self.frame_counter), # 0
                                   abs(self.amp1).copy(),          # 1
                                   abs(self.amp2).copy(),          # 2
                                   abs(self.amp3).copy())]         # 3

    def run(self, dt = 0.01, steps = 500, freq = 100, ptype = None):
        figure()
        self.plot(ptype)
        raw_input()
        while steps > 0:
            stdout.flush()
            self.step(dt)
            if steps % freq == 0:
                print 'remaining {0}'.format(steps)
                if len(self.saved_curves) > 100:
                    pnew = find_pos(self.zvec, self.saved_curves[-1][2])
                    pold = find_pos(self.zvec, self.saved_curves[-100][2])
                    tnew = self.saved_curves[-1][0]
                    told = self.saved_curves[-100][0]
                    print 'front velocity: {0}\n'.format((pold - pnew) / (tnew - told))
                stdout.flush()
                self.plot(ptype)                
            steps -= 1

    # def run_save(self,dt=0.01,steps=500):
    #     self.figs = [figure(figsize=(15,4)) for _ in range(3)]
    #     self.axes = [f.add_subplot(111) for f in self.figs]
    #     self.plot_save()
    #     while steps > 0:
    #         print steps,
    #         stdout.flush()
    #         self.step(dt)
    #         if steps%100==0:
    #             self.plot_save()
    #         steps-=1
    
    def plot_hex(self):
        #        figure()
        ar = \
          array([self.amp1]).transpose() * self.planew1 \
          + array([self.amp2]).transpose() * self.planew2 \
          + array([self.amp3]).transpose() * self.planew3 \
          + array([self.amp1.conj()]).transpose() * self.planew1_conj \
          + array([self.amp2.conj()]).transpose() * self.planew2_conj \
          + array([self.amp3.conj()]).transpose() * self.planew3_conj
        imshow(real(ar), interpolation = 'none', cmap=cm.PiYG)
        #        imshow(real(ar), interpolation = 'none', cmap=cm.spectral)
    
    def plot(self, ptype):
        if ptype == None:
            plot(self.zvec, abs(self.amp2))
            pause(0.1)
        elif ptype == 'all':
            clf()

            subplot2grid((6, 4), (0, 0), rowspan = 6)
            self.plot_hex()

            subplot2grid((6, 4), (0, 1), colspan = 3)
            plot(self.zvec, real(self.amp1))
            subplot2grid((6, 4), (1, 1), colspan = 3)
            plot(self.zvec, imag(self.amp1))

            subplot2grid((6, 4), (2, 1), colspan = 3)
            plot(self.zvec, real(self.amp2))
            subplot2grid((6, 4), (3, 1), colspan = 3)
            plot(self.zvec, imag(self.amp2))

            subplot2grid((6, 4), (4, 1), colspan = 3)
            plot(self.zvec, real(self.amp3))
            subplot2grid((6, 4), (5, 1), colspan = 3)
            plot(self.zvec, imag(self.amp3))

            tight_layout()
            pause(0.1)
        else:
            raise ValueError('not a proper plot type ')            
        
    # def plot(self, ptype):
    #     if ptype == None:
    #         plot(self.zvec, abs(self.amp2))
    #         pause(0.1)
    #     elif ptype == 'all':
    #         clf()

    #         subplot2grid((3, 4), (0, 0), rowspan = 3)
    #         self.plot_hex()

    #         subplot2grid((3, 4), (0, 1), colspan = 3)
    #         plot(self.zvec, abs(self.amp1))

    #         subplot2grid((3, 4), (1, 1), colspan = 3)
    #         plot(self.zvec, abs(self.amp2))

    #         subplot2grid((3, 4), (2, 1), colspan = 3)
    #         plot(self.zvec, abs(self.amp3))

    #         tight_layout()
    #         pause(0.1)
    #     else:
    #         raise ValueError('not a proper plot type ')            
                

    # def plot_save(self):
    #     for axis in self.axes: axis.cla()
    #     self.axes[0].imshow(real(self.y), interpolation='none', cmap=cm.spectral)
    #     self.axes[1].plot(real(self.y.sum(axis = 0) / float(self.y.shape[0])))
    #     self.axes[2].plot(real(ifft2(self.yp * self.kmask).sum(axis=0)) / float(self.y.shape[0]))
    #     for fig, prefix in zip(self.figs, ['1', '2', '3']):
    #         fig.savefig('movie/' + prefix + '.%05d.png' % self.frame_counter)
    #     pause(0.1)

# p4=pfc.P1Amp(psi0=0.28,lz=800,nz=1600,theta=0,init='tanh',w=7,preamp=0.8,falloff=0.02)
# [plot(curves[1]) for curves in p4.saved_curves[500::200]];

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

def find_pos(xvec, ar):
    diff = ar - n.roll(ar, 1)
    tfar = diff >= diff.max()
    maxindex = array(range(len(xvec)))[tfar][0]
    maxpos = xvec[tfar][0]
    mini = maxindex - 20
    maxi = maxindex + 20
    funar = interp1d(xvec[mini:maxi],diff[mini:maxi],kind='quadratic')
    return find_max(funar, maxpos)

# def find_pos(xvec, ar):
#     diff = (ar-n.roll(ar,1))
#     maxpos = xvec[diff >= n.max(diff)][0]
#     funar = interp1d(xvec,diff,kind='linear')
#     return find_max(funar, maxpos)

# map(lambda frame: find_pos(loc_max(frame)), data[11:350])

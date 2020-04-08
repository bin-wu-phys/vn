#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by Bin Wu on Tue Jun 11 19:33:26 2019
"""
from Code.eshow import *
from Code.kernel import *
from Code.observable import *
from Code.lattice import *
import numpy as np

class KinTran:
    """
    KinTran()

    A KinTran object solves kinetic transport.

    Parameters
    ----------
    t0 : float
        The initial time.
    FInit : np.array
        The initial distribution at t0.
    C : kernel object
    dt : float
        The time step.

    Attributes
    ----------
    F : ndarray
        The distribution F in kinetic transport in Fourier space.
    t : float
        The time variable.
    dt : float
        The time step.
    deltacs: list of floats
        The coefficients for adaptive time step taken from p. 717 of numerical recipies in c.

    """
    def __init__(self, t0=None, FInit=np.array([]), kern=None, dt=0.1, err=1.0e-3):
        self.t = t0
        self.F = FInit
        self.kern = kern
        self.dt = dt
        self.adaptive_dtQ = True

        self.deltacs = [37.0/378.0 - 2825.0/27648.0, 0.0, 250.0/621.0 - 18575.0/48384.0,
                        125.0/594.0 - 13525.0/55296.0, - 277.0/14336.0, 512.0/1771.0 - 0.25
                        ]

        self.err = err

    def set_adaptive_dtQ(self, ftQ):
        """
        Set whether to use the fixed time step
        :param ftQ: True or False
        """
        self.adaptive_dtQ = ftQ

    def __str__(self):
        """
        formulates the output of F.
        """
        return 'At t = {}, F = \n {}'.format(self.t, self.F)

    def setTimeStep(self, dt):
        """
        set the time step for time marching
        """
        self.dt = dt

    def setErr(self):
        self.err = err

    def RungeKutta(self):
        """
        Calculate F and delta using the 5th order Runge-Kutta method.

        returns
        ----------
        F: np.array
        The distribution F in Fourier space.
        err: float
        The estimate of the trancation error, defined as the average value of the values on all the grid points.
        """
        k1 = self.dt * self.kern.calc(self.t, self.F)
        k2 = self.dt * self.kern.calc(self.t + 0.2 * self.dt, self.F + 0.2 * k1)
        k3 = self.dt * self.kern.calc(self.t + self.dt * 0.3, self.F + k1 * 3.0 / 40.0 + k2 * 9.0 / 40.0)
        k4 = self.dt * self.kern.calc(self.t + self.dt * 0.6, self.F + k1 * 0.3 - k2 * 0.9 + k3 * 1.2)
        k5 = self.dt * self.kern.calc(self.t + self.dt,
                                      self.F - k1 * 11.0 / 54.0 + k2 * 2.5 - k3 * 70.0 / 27.0 + k4 * 35.0 / 27.0)
        k6 = self.dt * self.kern.calc(self.t + self.dt * 0.875, self.F + k1 * 1631.0 / 55296.0 + k2 * 175.0 / 512.0
                                      + k3 * 575.0 / 13824.0 + k4 * 44275.0 / 110592.0 + k5 * 253.0 / 4096.0)

        F = self.F + k1 * 37.0 / 378.0 + k3 * 250.0 / 621.0 + k4 * 125.0 / 594.0 + k6 * 512.0 / 1771.0

        delta = self.deltacs[0] * k1 + self.deltacs[1] * k2 + self.deltacs[2] * k3 \
                + self.deltacs[3] * k4 + self.deltacs[4] * k5 + self.deltacs[5] * k6

        dim = delta.shape
        nnonzero = 0.0
        Fmax = np.amax(np.abs(self.F))
        for n in range(dim[0]):
            for nv in range(dim[1]):
                for nr in range(dim[2]):
                    for nph in range(dim[3]):
                        if np.abs(self.F[n, nv, nr, nph])/Fmax > 1.0e-8:
                            delta[n, nv, nr, nph] = delta[n, nv, nr, nph]/self.F[n, nv, nr, nph]
                            nnonzero = nnonzero + 1.0
                        else:
                            delta[n, nv, nr, nph] = 0.0
        err = np.sum(np.abs(delta))/nnonzero
        print('At t = {}, average err = {}, max err = {} with dt = {}.'.format(self.t, err,
                                                                               np.amax(np.abs(delta)), self.dt))
        
        return F, err

    def nextTimeRK4(self):
        """
        calcuates F at nexttime step using the fourth-order Runge-Kuta.
        Here, we can either choose n or theta space. I choose theta space.
        """
        k1 = self.dt*self.kern.calc(self.t, self.F)
        k2 = self.dt*self.kern.calc(self.t + 0.5*self.dt,
                                 self.F + 0.5*k1)
        k3 = self.dt*self.kern.calc(self.t + 0.5*self.dt,
                                 self.F + 0.5*k2)
        k4 = self.dt*self.kern.calc(self.t + self.dt, self.F + k3)

        self.F = self.F + (k1 + k4)/6.0 + (k2 + k3)/3.0
        self.t = self.t + self.dt

    def nextTimeAdaptive(self):
        """
        calcuates F at nexttime step using the fourth-order Runge-Kuta.
        Here, we can either choose n or theta space. I choose theta space.
        """

        F, err = self.RungeKutta()

        while err > self.err:
            self.dt = 0.9 * self.dt * (self.err/err) ** 0.25
            F, err = self.RungeKutta()

        self.F = F
        self.t = self.t + self.dt
        self.dt = 0.9*self.dt*(self.err/err)**0.2
        # print('at t = {}, dt = {}\n'.format(self.t, self.dt))

    def calcvn(self):
        """
        Calculate vn's.
        """
        return self.kern.ob.calcvn(self.t, self.F)

    def nextTime(self):
        if self.adaptive_dtQ:
            self.nextTimeAdaptive()
        else:
            self.nextTimeRK4()

    def nextTimeTo(self, ts):
        if  self.t + self.dt - ts > 1e-3:
            self.dt = ts - self.t
        if self.adaptive_dtQ:
            self.nextTimeAdaptive()
        else:
            self.nextTimeRK4()

    #show results
    def showEnergyDensity(self):
        """
        Show the 3d and contour plots of the energy density.
        """
        #update the data
        self.kern.ob.setTime(self.t)
        self.kern.ob.interp(self.F)
        self.kern.ob.InterpRestFrame(self.kern.dr)
        
        #show the results
        rmax = self.t + 2.0
        dr = 0.1
        nr = int(rmax/dr)
        rLIST = np.linspace(0, nr*dr, nr+1)
        eLIST = np.array([[ self.kern.ob.e[nth](r)
        for nth in range(len(self.kern.ob.latt.lattice[self.kern.ob.idx_th]))
        ]
        for r in rLIST
        ]
        )
        ec = energyContour(rLIST, self.kern.ob.latt.lattice[self.kern.ob.idx_th], eLIST, self.t)
        e3D = energy3D(rLIST, self.kern.ob.latt.lattice[self.kern.ob.idx_th], eLIST, self.t)

        return rLIST, self.kern.ob.latt.lattice[self.kern.ob.idx_th], eLIST, ec, e3D

    #save and load
    def savekt(self, filename):
        np.savez_compressed(filename+'.kt', t=self.t, F=self.F, dt=self.dt)

    def save(self, filename):
        lattdim = np.array([self.kern.ob.idx_th, self.kern.ob.idx_v,
                            self.kern.ob.idx_r, self.kern.ob.idx_ph,
                            len(self.kern.ob.latt.lattice[self.kern.ob.idx_th]),
                            len(self.kern.ob.latt.lattice[self.kern.ob.idx_v]),
                            len(self.kern.ob.latt.lattice[self.kern.ob.idx_r]),
                            len(self.kern.ob.latt.lattice[self.kern.ob.idx_ph])
                            ])

        lattlim = np.array([self.kern.ob.latt.vzmin, self.kern.ob.latt.vzmax,
                            self.kern.ob.latt.lattice[self.kern.ob.idx_r][0],
                            self.kern.ob.latt.lattice[self.kern.ob.idx_r][-1],
                            ])
        np.savez_compressed(filename, t=self.t, F=self.F, dt=self.dt, err=self.err,
                            g=self.kern.g, dr=self.kern.dr, t0=self.kern.ob.t0,
                            lattdim=lattdim, lattlim=lattlim
                            )

    def loadkt(self, filename):
        loadkt = np.load(filename+'.kt.npz')
        self.t, self.F, self.dt= 1.0*loadkt['t'], loadkt['F'], 1.0*loadkt['dt']

    def load(self, filename):
        #load data
        if filename[-4:]=='.npz':
            loadkt = np.load(filename)
        else:
            loadkt = np.load(filename+'.npz')
        lattdim = loadkt['lattdim']
        lattlim = loadkt['lattlim']
        lattin = sorted(((lattdim[0], 'theta', lattdim[4]),
                        (lattdim[1], 'vz tilde', lattlim[0], lattlim[1], lattdim[5]), 
                        (lattdim[2], 'r tilde', lattlim[2], lattlim[3], lattdim[6]),
                        (lattdim[3], 'phir', lattdim[7]))
                        )

        #duplicate lattice
        latt = Lattice(lattin[0][1:], lattin[1][1:],
                       lattin[2][1:], lattin[3][1:])

        #duplicate observable
        ob = Observable(1.0*loadkt['t0'], latt)

        #load kern
        self.kern = Kernel(ob=ob, g=1.0*loadkt['g'], dr=1.0*loadkt['dr'])

        #load kt data
        self.t, self.F, self.dt = 1.0*loadkt['t'], loadkt['F'], 1.0*loadkt['dt']
        if 'err' in loadkt.files:
            self.err = 1.0 * loadkt['err']
        
        ob.setTime(self.t)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:45:42 2019

@author: biwu
"""

from math import sin, cos, sqrt, atan, pi, copysign, ceil
import numpy as np
from scipy import interpolate
from Code.proc import *
import os
import time
from Code.lattice import *


def arctan(x, y):
    value = None
    if x ==0 and y == 0:
        value = 0.0
    else:
        if x == 0.0:
            value = 0.5*copysign(1.0, y)*pi
        else:
            value = atan(y/x)
            if x < 0.0:
                if y == 0.0:
                    value = pi
                else:
                    value = copysign(1.0, y)*pi + value
    return value

def arctannp(x, y):
    if len(x) != len(y):
        print('x and y should have the same number of elements!\n')
        exit()
    return np.array([arctan(x[i], y[i]) for i in range(len(x))])

def gdd(xu, yu):
    """
    Returns g_{\mu\nu} x^\mu y^\nu.
    """
    res = np.zeros(len(xu[1,:]))
    if len(xu) != len(yu):
        raise  ValueError('The two vectors do not have the same len.')
    else:
        for i in range(len(res)):
            res[i] = np.dot(xu[1:, i], yu[1:, i]) - xu[0, i]*yu[0, i]
    return res

def heavisidetheta(x):
    #return 1.0*(np.real(x)>0)*(np.abs(np.imag(x))/np.real(x)<1.0e-6)
    return 1.0*(x>0.0)

def physicalQ(e, u):
    """
    Check whether the energy density or transverse pressure is physical.
    """
    return np.prod(np.imag(e)==0)*np.prod(np.imag(u)==0)

def setunphysical():
    return 0.0, np.array([1.0, 0.0, 0.0])

class Observable:
    """
    Calculate the energy density, vu, v_n.
    
    Attributes
    ----------
    idx_th : integer
        Index for $\theta$.
    idx_v : integer
        Index for $\tilde{v}_z$.
    idx_r : integer
        Index for $\tilde{r}$.
    idx_ph : integer
        Index for $\tilde{\phi}_r.
    phir : list
        phir as a function of $\tilde{v}_z, \tilde{r}$ and $\tilde{\phi}_r$.
    delta : list
        The quality: 1 - \tilde{v}_z^2 + (t_0 \tidle{v}_z/t)^2
    
    """
    #Comoving coordinate system
    
    def xu(self, xti, phi, vzti):
        vT2 = 1.0 - vzti*vzti
        return xti - cos(phi)*(self.t0 - self.t*sqrt(vT2 + (self.t0*vzti/self.t)**2))/sqrt(vT2)

    def yu(self, yti, phi, vzti):
        vT2 = 1.0 - vzti*vzti
        return yti - sin(phi)*(self.t0 - self.t*sqrt(vT2 + (self.t0*vzti/self.t)**2))/sqrt(vT2)

    def xtu(self, x, phi, vz):
        vT2 = 1.0 - vz*vz
        return x - cos(phi)*(self.t - self.t0*sqrt(vT2 + (self.t*vz/self.t0)**2))/sqrt(vT2)

    def ytu(self, y, phi, vz):
        vT2 = 1.0 - vz*vz
        return y - sin(phi)*(self.t - self.t0*sqrt(vT2 + (self.t*vz/self.t0)**2))/sqrt(vT2)

    def rtu(self, x, y, phi, vz):
        return sqrt(self.xtu(x, phi, vz)**2 + self.ytu(y, phi, vz)**2)

    def vzu(self, vzt):
        return self.t0*vzt/(self.t*sqrt(1.0 - vzt*vzt + (self.t0*vzt/self.t)**2))

    def vztu(self, vz):
        return self.t*vz/(self.t0*sqrt(1.0 - vz*vz + (self.t*vz/self.t0)**2))

    def thetatu(self, x, y, phi, vz):
        return arctan(self.xtu(x, phi, vz), self.ytu(y, phi, vz))

    def thetau(self, xt, yt, phi, vzt):
        return arctan(self.xu(xt, phi, vzt), self.yu(yt, phi, vzt))
    
    def ru(self, xt, yt, phi, vzt):
        return sqrt(self.xu(xt, phi, vzt)**2 + self.yu(yt, phi, vzt)**2)

    def calcphir(self):
        self.phir = np.array([[[[
                    np.exp(1.0j*n*(phit - self.thetau(rt, 0, phit, vzt)))
                    for phit in self.latt.lattice[self.idx_ph]
                    ]
                    for rt in self.latt.lattice[self.idx_r]
                    ]
                    for vzt in self.latt.lattice[self.idx_v]
                    ]
                    for n in range(len(self.latt.lattice[self.idx_th]))
                    ]
                    )
    
    def dphiIntegrator(self, f):
        """
        Integrate fi using the trapezoidal rule. fi is a periodic function of phi.
        """
        fb = 0.5*(f + np.concatenate((f[1:], [f[0]])))
        #print(fb)        
        return np.dot(fb, self.dphi)/(2.0*np.pi)

    def integratephi(self, f):
        """
        Integrate fi using the trapezoidal rule. fi is a periodic function of phi.
        """
        return  np.array([[[self.dphiIntegrator(f[:, nv, n, i])
                for i in range(6)]
                for n in range(self.dim_n)]
                for nv in range(len(self.latt.lattice[self.idx_v]))
                ])

    def dvIntegrator(self, f):
        """
        Integrate fi over vz using the trapezoidal rule.
        """
        fb = 0.5*(f + np.concatenate((f[1:], [f[-1]])))
        return np.dot(fb, self.dvz)

    def integratevz(self, f):
        """
        Integrate fi using the trapezoidal rule. fi is a periodic function of phi.
        """
        return  np.array([[self.dvIntegrator(f[:, n, i])
                for i in range(6)]
                for n in range(self.dim_n)
                ])

    def drIntegrator(self, f):
        rf = self.latt.lattice[self.idx_r]*f
        return np.dot(0.5*(rf[0:-1]+rf[1:]), self.dr)

    def initializetnn(self):
        """
        Calculate some component of v^\mu v^\nu.
        """
        self.vTt2 = 1.0 - self.latt.lattice[self.idx_v]*self.latt.lattice[self.idx_v]
        self.vTt = np.sqrt(self.vTt2)
        self.t0x = np.outer(self.vTt, np.cos(self.latt.lattice[self.idx_ph]))
        self.t0y = np.outer(self.vTt, np.sin(self.latt.lattice[self.idx_ph]))
        self.txx = np.outer(self.vTt2,
                            np.cos(self.latt.lattice[self.idx_ph])**2
                            )
        self.txy = np.outer(self.vTt2,
                            0.5*np.sin(2.0*self.latt.lattice[self.idx_ph])
                            )
        self.tyy = np.outer(self.vTt2,
                            np.sin(self.latt.lattice[self.idx_ph])**2
                            )

    def setLattice(self, latt):
        self.latt = latt
        if latt is None:
            self.idx_th, self.idx_v, self.idx_r, self.idx_ph = 0, 0, 0, 0
            self.vzT2 = self.vzT = []
            self.t0x = self.t0y = self.txx = self.txy = self.tyy = []
            self.dphi = self.dvz = self.dr = []
        else:
            self.idx_th, self.idx_v, self.idx_r, self.idx_ph = latt.indices()
            self.dim_n = len(self.latt.lattice[self.idx_th])//2 + 1
            self.dphi = np.concatenate((self.latt.lattice[self.idx_ph][1:], [2.0*np.pi])) - self.latt.lattice[self.idx_ph]
            dvz = self.latt.lattice[self.idx_v][1:] - self.latt.lattice[self.idx_v][0:-1]
            self.dvz = np.concatenate((dvz, [dvz[-1]]))
            self.dr = self.latt.lattice[self.idx_r][1:] - self.latt.lattice[self.idx_r][0:-1]
            self.initializetnn()
            self.e = [0 for n in range(len(self.latt.lattice[self.idx_th]))]
            self.u = [[] for n in range(len(self.latt.lattice[self.idx_th]))]

    def __init__(self, t0, latt=None, procs=1):
        self.t0 = t0
        self.t = t0
        self.Finterp = []
        self.e = []
        self.u = []
        self.setLattice(latt)
        self.phir = []
        self.warningQ = False
        self.delta = self.vT = []
        self.setProcs(procs)

    def setTime(self, t):
        self.t = t
    
    def update(self):
        """
        Update the time-dependent quantites. Run it after setTime.
        """
        self.delta = np.array([ self.vTt2[nvz] + (self.t0*self.latt.lattice[self.idx_v][nvz]/self.t)**2 for nvz in range(len(self.latt.lattice[self.idx_v]))])
        self.vT = self.vTt/np.sqrt(self.delta)

    def setWarning(Q):
        self.warningQ = Q

    def interpExt(self, F, rt, phi):
        """
        F is the distribution. rt is the list of \tilde{r}.
        phi is the list of \tilde(\phi) in [0, 2pi).
        """
        self.Finterp.clear()
        for n in range(F.shape[0]):
            self.Finterp.append([])
            for vz in range(F.shape[1]):
                self.Finterp[-1].append(
                        [interpolate.interp2d(
                                rt, phi, np.transpose(F[n, vz, :, :].real),
                                kind='linear', fill_value=0.0),
                        interpolate.interp2d(
                                rt, phi, np.transpose(F[n, vz, :, :].imag),
                                kind='linear', fill_value=0.0)]
                        )

    def interp(self, F):
        """
        Calcualte the interpolation function of F over \tilde{r} and \tilde{phi}_r.
        
        Parameters
        ----------
        F : np.array
            The distribution in n, \tilde{v}_z, \tilde{r} and \tilde{\phi}_r
            
        Returns
        ----------
            A 2-d list in n and \tilde{v}_z of the interpolation functions
            dependent of \tilde{r} \tilde(\phi) in [0, 2pi). Becaureful, 
            There is an additional square bracket in the result.
        """
        self.Finterp.clear()
        for n in range(F.shape[0]):
            self.Finterp.append([])
            for vz in range(F.shape[1]):
                self.Finterp[-1].append(
                        [interpolate.interp2d(
                                self.latt.lattice[self.idx_r],
                                self.latt.lattice[self.idx_ph],
                                np.transpose(F[n, vz, :, :].real),
                                kind='linear', fill_value=0.0),
                        interpolate.interp2d(
                                self.latt.lattice[self.idx_r],
                                self.latt.lattice[self.idx_ph],
                                np.transpose(F[n, vz, :, :].imag),
                                kind='linear', fill_value=0.0)]
                        )

    def TuuIntegrand(self, n, r, nv, nphi):
        """
        Calculate the components for T^{\mu\nu}.
        
        Parameters
        ----------
        n : an integer
            It is the index of theta in Fourier space.
        nv : an integer
            It is the the index of \tilde{v}_z.
        r : a float number
            It is the physical radius r.
        nphi : a float number in [0, 2\pi)
            It is the index for $\phi_r = \phi - \theta$.
            Here, we integrate over $phi_r$ but F is a function of $\tilde{\phi}_r$.
            
        Returns
        -------
            The integrands for 
            (t/t^0)*[T^{00}, T^{0x}, T^{0y}, T^{xx}, T^{xy}, T^{yy}].
        """
        vz = self.vzu(self.latt.lattice[self.idx_v][nv])
        m = sqrt(self.vTt2[nv] + (self.t0*self.latt.lattice[self.idx_v][nv]/self.t)**2)
        phit = self.latt.lattice[self.idx_ph][nphi] - self.thetatu(r, 0, 
                                self.latt.lattice[self.idx_ph][nphi],
                                vz)
        rt = self.rtu(r, 0.0,
                      self.latt.lattice[self.idx_ph][nphi], vz
                      )

        return np.array([m, self.t0x[nv, nphi], self.t0y[nv, nphi],
                        self.txx[nv, nphi]/m, self.txy[nv, nphi]/m,
                        self.tyy[nv, nphi]/m]
                        )*np.exp(1.0j*n*self.latt.lattice[self.idx_ph][nphi]
                        )*(self.Finterp[n][nv][0](rt, phit%(2.0*np.pi))[0]
                        +1.0j*self.Finterp[n][nv][1](rt, phit%(2.0*np.pi))[0])

    def Tuu(self, r, theta):
        """
        Evalute all the components of $T^{\mu\nu}$ with $\theta = 2\pi n/N$.
        """
        itgtab = np.array([[[self.TuuIntegrand(n, r, nv, nphi)
                for n in range(self.dim_n)]
                for nv in range(len(self.latt.lattice[self.idx_v]))]
                for nphi in range(len(self.latt.lattice[self.idx_ph]))
                ], dtype=np.complex)
        vzitg = self.integratephi(itgtab)
        tnuu = self.integratevz(vzitg)
        tnn = np.real(tnuu[0]+np.exp(1.0j*theta*(self.dim_n-1))*tnuu[-1])
        for n in range(1, self.dim_n-1):
            tnn += 2.0*np.real(np.exp(1.0j*n*theta)*tnuu[n])

        #return tnn
        return np.array([[tnn[0], tnn[1], tnn[2], 0.0],
                         [tnn[1], tnn[3], tnn[4], 0.0],
                         [tnn[2], tnn[4], tnn[5], 0.0],
                         [0.0, 0.0, 0.0, tnn[0] - tnn[3] - tnn[5]]
                         ])

    def RestFrame(self, r, theta):
        """
        Calucate the local energy density etc.
        """
        #print('({}, {})'.format(r, theta))
        e = 0
        u = np.array([])
        Tuu = self.Tuu(r, theta)
        tud = Tuu[0:-1, 0:-1]
        tud[:, 0] = -tud[:, 0]
        #print(tud)
        ee, uu = np.linalg.eig(tud)
        #print(ee)
        #print(uu)
        if physicalQ(ee, uu):
            udu = gdd(uu,uu)
            idx = np.nonzero(heavisidetheta(-udu))
            #print(len(idx),idx)
            if len(idx[0]) != 1:
                e, u = setunphysical()
                if self.warningQ:
                    print('Oops! One has {} energy eigenvalues!'.format(len(idex)))
                    print('Warning: at t = {} one has {} energy density at r = {} and theta = {}'.format(self.t, len(idex), r, theta))
            else:
                if ee[idx[0][0]] >= 0.0:
                    e, u = setunphysical()
                else:
                    e = -ee[idx[0][0]]*self.t0/self.t
                    u = uu[:, idx[0][0]]/sqrt(-udu[idx[0][0]])
                    u = u.reshape([3])
                    #print('u = ', u)
                    if u[0] < 0.0:
                        u = -u
        else:
            e, u = setunphysical()
            if self.warningQ:
                print('''Warning:
                        at t = {} energy density becomes unphysics r = {} and theta = {}
                        '''.format(self.t, r, theta)
                        )

        return e/len(self.latt.lattice[self.idx_th]), u

    def InterpRestFramen(self, n, rLIST):
        """
        Interpolate e and u.
        
        Parameters
        ----------
        dr : float
            r spacing for the interpolation.        
        """
#        print('{}: {}'.format(n, rLIST))
        #print(os.getpid())
        eLIST = np.zeros(len(rLIST))
        uLIST = np.zeros([len(rLIST), 3])

        theta = self.latt.lattice[self.idx_th][n]

        for nr in range(len(rLIST)):
            eLIST[nr], uLIST[nr, :] = self.RestFrame(rLIST[nr], theta)

        en = interpolate.interp1d(rLIST, eLIST,
                                          kind='linear', fill_value=0.0)

        un = [interpolate.interp1d(rLIST, uLIST[:, 0], kind='linear', fill_value=0.0),
                      interpolate.interp1d(rLIST, uLIST[:, 1], kind='linear', fill_value=0.0),
                      interpolate.interp1d(rLIST, uLIST[:, 2], kind='linear', fill_value=0.0),
                     ]
        return en, un

    def setProcs(self, procs):
        self.procs = procs

        if self.latt != None:
            # calculate the asignment table for the processes
            nmax = len(self.latt.lattice[self.idx_th])
            npp, res = nmax//self.procs, nmax%self.procs

            if res > 0:
                self.nrange = [range(0, npp + 1)]
            else:
                self.nrange = [range(0, npp)]

            for np in range(1, self.procs):
                if np < res:
                    self.nrange.append(range(self.nrange[np-1][-1] + 1, self.nrange[np-1][-1] + npp + 2))
                else:
                    self.nrange.append(range(self.nrange[np-1][-1] + 1, self.nrange[np-1][-1] + npp + 1))
        else:
            self.nrange = None

    def InterpRestFrame(self, dr):
        """
        Interpolate e and u. It does the job but it does not save any time.
        
        Parameters
        ----------
        dr : float
            r spacing for the interpolation.        
        """
        RMAX = self.latt.lattice[self.idx_r][-1] + (self.t - self.t0) + dr 
        nr = ceil(RMAX/dr)
        rLIST = np.linspace( 0, nr*dr, nr+1)
        if self.procs > 1:
            q = Queue()
            ps = []

            if self.nrange == None:
                self.setProcs(self.procs)

            for n in range(1, self.procs):
                p = Proc(self.nrange[n], rLIST, self, q)
                ps.append(p)
                p.start()

            pid = os.getpid()
            for n in self.nrange[0]:
                self.e[n], self.u[n] = self.InterpRestFramen(n, rLIST)
                # print('Main {} has done with {}.\n'.format(pid, n))

            nfp, nfps = 0, len(self.latt.lattice[self.idx_th]) - len(self.nrange[0])
            while nfp < nfps:
                res = q.get()
                # print('{} gets one!'.format(pid))
                self.e[res[0]], self.u[res[0]] = res[1], res[2]
                nfp = nfp + 1
        else:
            for n in range(len(self.latt.lattice[self.idx_th])):
                self.e[n], self.u[n] = self.InterpRestFramen(n, rLIST)


            
    def integrateFn(self, Fn):
        """
        Integrate out $\tilde{\phi}, \tidle r$ and $\tilde{v}_z$.
        """
        itg_ph = np.array([[ self.dphiIntegrator(Fn[nvz, nr, :])
                for nr in range(len(self.latt.lattice[self.idx_r]))
                ]
                for nvz in range(len(self.latt.lattice[self.idx_v]))
                ])
        itg_r = np.array([ self.drIntegrator(itg_ph[nvz, :])
                for nvz in range(len(self.latt.lattice[self.idx_v]))
                ])
        return self.dvIntegrator(self.vTt*itg_r)
            
    def calcvn(self, t, F):
        """
        Calculate vn's: cos(n\phi) and sin(n\phi)
        """
        vn = [t]
        v0 = self.integrateFn(F[0, :, :, :])

        for n in range(1, self.dim_n):
            rn = self.integrateFn(F[n, :, :, :])/v0
            vn.append(rn.real)
            vn.append(rn.imag)

        return vn

if __name__ == '__main__':

    xti, yti, phi, vzti, t, t0 = 1.0, 0.5, 0.3, 0.1, 2.0, 0.1
    latt = Lattice(('r tilde', 0.0, 3.0, 31),
                   ('vz tilde', 0.0, 10.0, 21),
                   ('phir', 40),
                   ('theta', 10))

    ob = Observable(t0=t0, latt=latt)
    ob.setTime(t)
    ob.setProcs(4)
    print(ob.nrange)
    #print(xti, phi, vzti, t, t0)
    #print(xu(xti, phi, vzti, t, t0), yu(xti, phi, vzti, t, t0))
    #print(xtu(xti, phi, vzti, t, t0), ytu(xti, phi, vzti, t, t0))
    print(ob.thetatu(xti, yti, phi, vzti))
    print(ob.thetau(xti, yti, phi, vzti))
    #print(ru(xti, yti, phi, vzti, t, t0))
    #print(ob.vzu(0.5, 3.0))
    #print(ob.vztu(0.5, 3.0))
    #for phi in 2.0*pi*np.arange(30)/30.0:
    #    print(arctan(cos(phi), sin(phi)))



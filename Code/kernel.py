#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 12:11:13 2019

@author: biwu
"""
import numpy as np
#from observable import *

class Kernel:
    """
    Kernel in the isotropization time approximation (ITA).
    
    Attributes
    ----------
    dr : float
        The r spacing, used for the energy density e and u interpolation. 
    """
    def __init__(self, g=1.0, ob=None, dr = 0.1):
        self.g = g
        self.ob = ob
        self.dr = dr

    def setObservable(self, ob):
        self.ob = ob

    def setCoupling(self, g):
        self.g = g

    def setdr(self, dr):
        """
        set the time step for time marching
        """
        self.dr = dr
        
    def Ftheta(self, F, nvt, nrt, nphit):
        """
        Calculate F in theta space.

        Parameters
        ----------
        nvt, nrt, nphit : integers
        indices respective for $\tilde{v}_z$, $\tilde{r}$ and $\tilde{\phi}_r$.
        
        Returns
        ----------
        Inverse fourier transform of $\tilde{f}_{(n)}$.
        """

        fn = np.concatenate((F[:, nvt, nrt, nphit],
                                           np.conjugate(F[-2:0:-1, nvt, nrt, nphit]))
                                         )

        return np.fft.ifft(self.ob.phir[:, nvt, nrt, nphit]*fn)

    def Ctheta(self, t, F, nvt, nrt, nphit):
        """
        Collision kernel in theta space.

        Parameters
        ----------
        t : float
            the time variable.
        F : list of floats
            the distribution function in theta space.
        nvt, nrt, nphit : integers
            indices respective for $\tilde{v}_z$, $\tilde{r}$ and $\tilde{\phi}_r$.
        
        Returns
            list of floats:
                collision kernel in theta space for given $\tilde{v}_z$, $\tilde{r}$ and $\tilde{\phi}_r$.
        ----------
                    
        """
        r = self.ob.ru(self.ob.latt.lattice[self.ob.idx_r][nrt], 0.0,
                       self.ob.latt.lattice[self.ob.idx_ph][nphit],
                       self.ob.latt.lattice[self.ob.idx_v][nvt]
                       )

        phir = self.ob.latt.lattice[self.ob.idx_ph][nphit] - self.ob.thetau(
                        self.ob.latt.lattice[self.ob.idx_r][nrt], 0.0,
                        self.ob.latt.lattice[self.ob.idx_ph][nphit],
                        self.ob.latt.lattice[self.ob.idx_v][nvt]
                        )

        v = np.array([1.0, self.ob.vT[nvt]*np.cos(phir), self.ob.vT[nvt]*np.sin(phir)])
 
        u = np.array([[ self.ob.u[nth][mu](r)
                        for mu in range(3)
                        ] 
                        for nth in range(len(self.ob.e))
                        ]
                    )
        vu = v[0]*u[:, 0] - v[1]*u[:, 1] - v[2]*u[:, 2]       
        
        e = np.array([self.ob.e[nth](r) for nth in range(len(self.ob.e))])
 
        Fth = self.Ftheta(F, nvt, nrt, nphit)
        
        return self.g*(e**1.25/((vu**3)*(self.ob.delta[nvt]**2)) - (e**0.25)*vu*Fth)

    def calc(self, t, F):
        """
        The collision kernel in Fourier space.

        Parameters
        ----------
        t : float
            the time variable.
        F : list of floats
            the distribution function in theta space.
        
        Returns
            list of floats:
                collision kernel in fourier space.        
        """
        self.ob.interp(F)
        self.ob.setTime(t)
        self.ob.update()
        self.ob.calcphir()
        self.ob.InterpRestFrame(self.dr)
        C = np.zeros(F.shape, dtype=np.complex)
        for nvt in range(len(self.ob.latt.lattice[self.ob.idx_v])):
            for nrt in range(len(self.ob.latt.lattice[self.ob.idx_r])):
                for nphit in range(len(self.ob.latt.lattice[self.ob.idx_ph])):
                    Cn = np.conjugate(self.ob.phir[:self.ob.dim_n, nvt, nrt, nphit]
                            )*(np.fft.fft(self.Ctheta(t, F, nvt, nrt, nphit))[:self.ob.dim_n])
                    for n in range(self.ob.dim_n):
                        C[n, nvt, nrt, nphit] = Cn[n]
        
        return C

    def save(self, filename):
        np.savez_compressed(filename+'.kern', g=self.g, dr=self.dr)

    def load(self, filename):
        loadkt = np.load(filename+'.kern.npz')
        self.g, self.dr= 1.0*loadkt['g'], 1.0*loadkt['dr']

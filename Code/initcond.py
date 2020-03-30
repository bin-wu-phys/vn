#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:46:10 2019

@author: biwu
"""

import numpy as np
#import lattice
from scipy import interpolate, integrate
from scipy.special import jv, jn_zeros


class InitCond:
    """
    An InitCond object sets up the initial condition.

    Parameters
    ----------
    type : string
        The type of the initial condition.
    t0 : float
        The initial time
    harmonics : tuple
        Its element is a tuple or list of the form (n, vn).

    Attributes
    ICtype : string
        The type of transverse profile: 'guassian' or 'hocky puck'.
    t0 : float
        The initial time
    vn : tuple of {int, float}
        The theta harmonics and their overall
    ----------
    """
    def Fbg(self, r, phir):
        """
        The background transverse profile.
        """
        return 2.0*np.exp(-r**2-self.t0**2+2.0*r*self.t0*np.cos(phir))/self.t0

    def J(self, n, l, r2):
        """
        Bessel-Gaussian parametrization of the radial profile.
        n : nth theta harmonic
        l: lth radial wavenumber
        r2: r^2
        :return: lth radial wavenumber of the nth theta harmonic.
        """
        rho = np.sqrt(1.0 - np.exp(-r2))

        return jv(n, jn_zeros(n, l)[-1]*rho)/self.t0


    def FTn(self, r, phir):
        res = np.ndarray(len(self.dn), dtype=np.complex)
        r2 = (r**2+self.t0**2-2.0*r*self.t0*np.cos(phir))
        for n in range(len(self.dn)):
            phase = np.exp(-1.0j*self.dn[n][0]*self.dn[n][2])
            if self.dn[n][0] == 1:
                res[n] = 0.5*phase*(r - self.t0*np.cos(phir)
                        - 1.0j*self.t0*np.sin(phir))
                if not self.ancientQ:
                    """
                    max e_1 = 0.413418
                    """
                    #d_2 = 3.375*epsilon_2
                    res[n] = 2.25*np.sqrt(np.pi)*res[n]*np.exp(-0.5*r2)
            elif self.dn[n][0] == 2:
                res[n] = 0.5*phase*(r*r - 2.0*r*self.t0*np.cos(phir)
                        + self.t0*self.t0*np.cos(2.0*phir) - 2.0j*self.t0*(r
                        - self.t0*np.cos(phir))*np.sin(phir))
                if not self.ancientQ:
                    #d_2 = 3.375*epsilon_2
                    res[n] = 3.375*res[n]*np.exp(-0.5*r2)
            elif self.dn[n][0] == 3:
                res[n] = 0.5*phase*(r**3 - 3.0*r*r*self.t0*np.cos(phir)
                        + 3.0*r*(self.t0**2)*np.cos(2.0*phir)
                        - (self.t0**3)*np.cos(3.0*phir) - 1.0j*self.t0*(
                        3.0*(r**2)*np.sin(phir)- 3.0*self.t0*r*np.sin(
                        2.0*phir) + (self.t0**2)*np.sin(3.0*phir)))
                if not self.ancientQ:
                    """
                    the maximum e_3 has a value smaller than 0.38.
                    """
                    # print(81.0*np.sqrt(np.pi)/64.0)
                    res[n] = res[n]*81.0*np.sqrt(np.pi)*np.exp(-0.5*r2)/64.0
            elif self.dn[n][0] == 4:
                res[n] = 0.5 *phase*(r ** 4 - 4.0 * (r**3) * self.t0 * np.cos(phir)
                                + 6.0 * (r**2) * (self.t0 ** 2) * np.cos(2.0 * phir)
                                - 4.0*r*(self.t0**3)*np.cos(3.0*phir)
                                + (self.t0 ** 4) * np.cos(4.0 * phir)
                                - 1.0j * self.t0 * (
                                                    4.0 * (r ** 3) * np.sin(phir)
                                                    - 6.0 * self.t0 * (r**2) * np.sin(2.0 * phir)
                                                    + 4.0*r*(self.t0 ** 2) * np.sin(3.0 * phir)
                                                    - (self.t0**3)*np.sin(4.0*phir)
                                                    )
                                )
                if not self.ancientQ:
                    """
                    the maximum e_4 has a value smaller than 0.36.
                    """
                    # print(81.0*np.sqrt(np.pi)/64.0)
                    res[n] = 81.0*res[n] * np.exp(-0.5 * r2)/64.0
            elif self.dn[n][0] == 5:
                res[n] = 0.5 *phase*(r ** 5 - 5.0 * (r**4) * self.t0 * np.cos(phir)
                                + 10.0 * (r**3) * (self.t0 ** 2) * np.cos(2.0 * phir)
                                - 10.0*(r**2)*(self.t0 ** 3) * np.cos(3.0 * phir)
                                + 5.0*r*(self.t0**4)*np.cos(4.0*phir)
                                - (self.t0**5)*np.cos(5.0*phir)
                                - 1.0j * self.t0 * (
                                                    5.0 * (r ** 4)* np.sin(phir)
                                                    - 10.0 * self.t0 * (r**3) * np.sin(2.0 * phir)
                                                    + 10.0*(r**2)*(self.t0 ** 2) * np.sin(3.0 * phir)
                                                    - 5.0*r*(self.t0**3)*np.sin(4.0*phir)
                                                    + (self.t0**4)*np.sin(5.0*phir)
                                                    )
                                )
                if not self.ancientQ:
                    """
                    the maximum e_5 has a value smaller than 0.345.
                    """
                    # print(81.0*np.sqrt(np.pi)/64.0)
                    res[n] = res[n] * 729.0 * np.sqrt(np.pi) * np.exp(-0.5 * r2) / 2048.0
            elif self.dn[n][0] == 6:
                res[n] = 0.5 * phase * (r ** 6
                                        - 6.0 * (r ** 5) * self.t0 * np.cos(phir)
                                        + 15.0 * (r ** 4) * (self.t0 ** 2) * np.cos(2.0 * phir)
                                        - 20.0 * (r ** 3) * (self.t0 ** 3) * np.cos(3.0 * phir)
                                        + 15.0 * (r**2) * (self.t0 ** 4) * np.cos(4.0 * phir)
                                        - 6.0*r*(self.t0 ** 5) * np.cos(5.0 * phir)
                                        + (self.t0**6)*np.cos(6.0 * phir)
                                        - 1.0j * self.t0 * (
                                                6.0 * (r ** 5) * np.sin(phir)
                                                - 15.0 * self.t0 * (r ** 4) * np.sin(2.0 * phir)
                                                + 20.0 * (self.t0**2)*(r ** 3) * np.sin(3.0 * phir)
                                                - 15.0 * (self.t0**3)*(r**2)  * np.sin(4.0 * phir)
                                                + 6.0*(self.t0 ** 4)*r * np.sin(5.0 * phir)
                                                - (self.t0 ** 5) * np.sin(6.0 * phir)
                                                )
                                        )
                if not self.ancientQ:
                    """
                    the maximum e_6 has a value smaller than 0.326544.
                    """
                    # print(81.0*np.sqrt(np.pi)/64.0)
                    res[n] = res[n] * 729.0 * np.exp(-0.5 * r2) / 2560.0
        return res

    def FTn_spatial(self, r, phir):
        res = np.ndarray(len(self.dn), dtype=np.complex)
        r2 = (r**2+self.t0**2-2.0*r*self.t0*np.cos(phir))
        for n in range(len(self.dn)):
            phase = np.exp(-1.0j*self.dn[n][0]*self.dn[n][2])
            if self.dn[n][0] == 1:
                res[n] = 0.5*phase*(r - self.t0*np.cos(phir)
                        - 1.0j*self.t0*np.sin(phir))
                if not self.ancientQ:
                    """
                    max e_1 = 0.413418
                    """
                    #d_2 = 3.375*epsilon_2
                    if self.dn[n][3] == 0:
                        res[n] = 2.25*np.sqrt(np.pi)*res[n]*np.exp(-0.5*r2)
            elif self.dn[n][0] == 2:
                res[n] = 0.5*phase*(r*r - 2.0*r*self.t0*np.cos(phir)
                        + self.t0*self.t0*np.cos(2.0*phir) - 2.0j*self.t0*(r
                        - self.t0*np.cos(phir))*np.sin(phir))
                if not self.ancientQ:
                    #d_2 = 3.375*epsilon_2
                    if self.dn[n][3] == 0:
                        res[n] = 3.375*res[n]*np.exp(-0.5*r2)
            elif self.dn[n][0] == 3:
                res[n] = 0.5*phase*(r**3 - 3.0*r*r*self.t0*np.cos(phir)
                        + 3.0*r*(self.t0**2)*np.cos(2.0*phir)
                        - (self.t0**3)*np.cos(3.0*phir) - 1.0j*self.t0*(
                        3.0*(r**2)*np.sin(phir)- 3.0*self.t0*r*np.sin(
                        2.0*phir) + (self.t0**2)*np.sin(3.0*phir)))
                if not self.ancientQ:
                    """
                    the maximum e_3 has a value smaller than 0.38.
                    """
                    # print(81.0*np.sqrt(np.pi)/64.0)
                    if self.dn[n][3] == 0:
                        res[n] = res[n]*81.0*np.sqrt(np.pi)*np.exp(-0.5*r2)/64.0
            elif self.dn[n][0] == 4:
                res[n] = 0.5 *phase*(r ** 4 - 4.0 * (r**3) * self.t0 * np.cos(phir)
                                + 6.0 * (r**2) * (self.t0 ** 2) * np.cos(2.0 * phir)
                                - 4.0*r*(self.t0**3)*np.cos(3.0*phir)
                                + (self.t0 ** 4) * np.cos(4.0 * phir)
                                - 1.0j * self.t0 * (
                                                    4.0 * (r ** 3) * np.sin(phir)
                                                    - 6.0 * self.t0 * (r**2) * np.sin(2.0 * phir)
                                                    + 4.0*r*(self.t0 ** 2) * np.sin(3.0 * phir)
                                                    - (self.t0**3)*np.sin(4.0*phir)
                                                    )
                                )
                if not self.ancientQ:
                    """
                    the maximum e_4 has a value smaller than 0.36.
                    """
                    # print(81.0*np.sqrt(np.pi)/64.0)
                    if self.dn[n][3] == 0:
                        res[n] = 81.0*res[n] * np.exp(-0.5 * r2)/64.0
            elif self.dn[n][0] == 5:
                res[n] = 0.5 *phase*(r ** 5 - 5.0 * (r**4) * self.t0 * np.cos(phir)
                                + 10.0 * (r**3) * (self.t0 ** 2) * np.cos(2.0 * phir)
                                - 10.0*(r**2)*(self.t0 ** 3) * np.cos(3.0 * phir)
                                + 5.0*r*(self.t0**4)*np.cos(4.0*phir)
                                - (self.t0**5)*np.cos(5.0*phir)
                                - 1.0j * self.t0 * (
                                                    5.0 * (r ** 4)* np.sin(phir)
                                                    - 10.0 * self.t0 * (r**3) * np.sin(2.0 * phir)
                                                    + 10.0*(r**2)*(self.t0 ** 2) * np.sin(3.0 * phir)
                                                    - 5.0*r*(self.t0**3)*np.sin(4.0*phir)
                                                    + (self.t0**4)*np.sin(5.0*phir)
                                                    )
                                )
                if not self.ancientQ:
                    """
                    the maximum e_5 has a value smaller than 0.345.
                    """
                    # print(81.0*np.sqrt(np.pi)/64.0)
                    if self.dn[n][3] == 0:
                        res[n] = res[n] * 729.0 * np.sqrt(np.pi) * np.exp(-0.5 * r2) / 2048.0
            elif self.dn[n][0] == 6:
                res[n] = 0.5 * phase * (r ** 6
                                        - 6.0 * (r ** 5) * self.t0 * np.cos(phir)
                                        + 15.0 * (r ** 4) * (self.t0 ** 2) * np.cos(2.0 * phir)
                                        - 20.0 * (r ** 3) * (self.t0 ** 3) * np.cos(3.0 * phir)
                                        + 15.0 * (r**2) * (self.t0 ** 4) * np.cos(4.0 * phir)
                                        - 6.0*r*(self.t0 ** 5) * np.cos(5.0 * phir)
                                        + (self.t0**6)*np.cos(6.0 * phir)
                                        - 1.0j * self.t0 * (
                                                6.0 * (r ** 5) * np.sin(phir)
                                                - 15.0 * self.t0 * (r ** 4) * np.sin(2.0 * phir)
                                                + 20.0 * (self.t0**2)*(r ** 3) * np.sin(3.0 * phir)
                                                - 15.0 * (self.t0**3)*(r**2)  * np.sin(4.0 * phir)
                                                + 6.0*(self.t0 ** 4)*r * np.sin(5.0 * phir)
                                                - (self.t0 ** 5) * np.sin(6.0 * phir)
                                                )
                                        )
                if not self.ancientQ:
                    """
                    the maximum e_6 has a value smaller than 0.326544.
                    """
                    # print(81.0*np.sqrt(np.pi)/64.0)
                    if self.dn[n][3] == 0:
                        res[n] = res[n] * 729.0 * np.exp(-0.5 * r2) / 2560.0

            # Spatial profile
            if self.dn[n][3] == 0:
                """
                Here, we put in the background such that it is the transverse profile just as the case self.dn[n][3] >0.
                """
                res[n] = res[n]*self.Fbg(r, phir)
            else:
                res[n] = self.J(self.dn[n][0], self.dn[n][3], r2)*res[n] / np.sqrt(r2)**self.dn[n][0]

        return res

    def Fn(self, theta, r, phir):
        """
        The nth-harmonic mode divided by Fbg at t = t0.
        """
        if self.dn[0][0] == 2:
            res = self.FTn(r, phir)*np.exp(self.dn[0][0]*1.0j*theta)
        else:
            print('The {}th harmonic has not been implemented!\n'.format
                  (self.dn[0][0]))
            raise ValueError('Unimplemented initial conditions.')
        return 2.0*res.real

    def gausstheta(self, latt):
        """
        gassian initial condition: FT*delta(vz)/(2pi) in theta without
        taking out exp(i n phi_r), that is, it is the distribution in
        ordinary coordinates.

        Parameters
        ----------
        latt : object of class Lattice
            Grids points for the initial conditon.
        """
        self.theta, self.vz, self.r, self.phi = latt.indices()
        F0 = np.zeros([latt.shape[self.theta], latt.shape[self.vz],
                       latt.shape[self.r], latt.shape[self.phi]],
                       dtype=np.complex)
        vz = 0
        for theta in range(latt.shape[self.theta]):
            for r in range(latt.shape[self.r]):  
                for phi in range(latt.shape[self.phi]):    
#                    print('theta = {}, r = {}, phi = {}, vz = {}'.
#                          format(theta, r, phi, vz))
                    F0[theta][vz][r][phi] = (1.0 + self.dn[0][1]
                        *self.Fn(latt.lattice[self.theta][theta], latt.lattice[self.r][r],
                                 latt.lattice[self.phi][phi])) * \
                        self.Fbg(latt.lattice[self.r][r], latt.lattice[self.phi][phi])\
                        / (latt.lattice[self.vz][1])

        return F0

    def fft(self, fun, latt):
        """
        initial condition from file.

        Parameters
        ----------
        latt : object of class Lattice
            Grids points for the initial conditon.
        """
        F0 = fun(latt)

        F0fft = np.empty(F0.shape, dtype=complex)
        for vzt in range(latt.shape[self.vz]):
            for rt in range(latt.shape[self.r]):
                for phirt in range(latt.shape[self.phi]):
                    F0fft[:, vzt, rt, phirt] = np.exp(-1.0j * np.array(range(latt.shape[self.theta])
                                                                       ) * latt.lattice[self.phi][phirt]
                                                    )*np.fft.fft(F0[:, vzt, rt, phirt])

        return F0fft[:latt.shape[self.theta] // 2 + 1, :, :, :]

    def gaussfft(self, latt):
        """
        gassian initial condition: FT*delta(vz)/(2pi) in theta.

        Parameters
        ----------
        latt : object of class Lattice
            Grids points for the initial conditon.
        """
        return self.fft(self.gausstheta, latt)

    def gauss(self, latt):
        """
        gassian initial condition: FT*delta(vz)/(2pi) in Fourier space.

        Parameters
        ----------
        latt : object of class Lattice
            Grids points for the initial conditon.
        """
        self.theta, self.vz, self.r, self.phi = latt.indices()
        F0 = np.zeros([latt.shape[self.theta]//2 + 1, latt.shape[self.vz],
                       latt.shape[self.r], latt.shape[self.phi]], dtype='complex')
            
        vz = 0

#        for r in range(latt.shape[self.r]):  
#            for phi in range(latt.shape[self.phi]):    
#                F0[0][vz][r][phi] = len(latt.lattice[self.theta])*self.Fbg(latt.lattice[self.r][r], latt.lattice[self.phi][phi])\
#                        / (latt.lattice[self.vz][1])
#                F0[self.dn[0][0]][vz][r][phi] = len(latt.lattice[self.theta])*self.dn[0][1]*self.FTn(latt.lattice[self.r][r],
#                                 latt.lattice[self.phi][phi]) * \
#                        self.Fbg(latt.lattice[self.r][r], latt.lattice[self.phi][phi])\
#                        / (latt.lattice[self.vz][1])
                        
        for r in range(latt.shape[self.r]):  
            for phi in range(latt.shape[self.phi]): 
                F0[0][vz][r][phi] = len(latt.lattice[self.theta])*self.Fbg(
                        latt.lattice[self.r][r], latt.lattice[self.phi][phi]
                        )/(latt.lattice[self.vz][1])
                fn = self.FTn_spatial(latt.lattice[self.r][r],
                              latt.lattice[self.phi][phi])
                for n in range(len(self.dn)):
                    #print(n, self.dn[n][0], self.dn[n][1])
                    F0[self.dn[n][0]][vz][r][phi] = np.exp(
                    -1.0j*self.dn[n][0]*latt.lattice[self.phi][phi]
                    )*self.dn[n][1]*fn[n]*len(latt.lattice[self.theta])/(latt.lattice[self.vz][1])

        return F0

    def hockypuck(self, latt):
        pass

    def ws(self, latt):
        pass

    def set_shiftQ(self, s):
        """
        set the value of shiftQ
        :param s: a bool
        """
        self.shiftQ = s

    def smoother(self, r2):
        """
        Smooth function.
        :param r2: a float or np.array, radial coordinate squaired
        :return: a float or np.array
        """
        return np.exp(-0.5*r2/self.R_smooth2)/(2.0*np.pi*self.R_smooth2)

    def read_density(self):
        """
        Read out transverse energy from fname.
        """
        f = open(self.fname)
        data = [l.replace('\n', '').split(' ') for l in f if ('#' not in l) and (l != '')]
        density = np.array([[float(d) for d in l] for l in data])
        xMax, dx = 10.0, 0.2
        xList = [-xMax + 0.5 * dx + n * dx for n in range(int(2.0 * xMax / dx))]
        yMax, dy = 10.0, 0.2
        yList = [-yMax + 0.5 * dy + n * dy for n in range(int(2.0 * yMax / dy))]

        if self.smoothQ:
            density = dx*dx*np.array(
                        [[ np.sum(
                            [[  density[ixp, iyp]*self.smoother(
                                    (xList[ix] - xList[ixp])**2 + (yList[iy] - yList[iyp])**2
                                )
                                for iyp in range(len(yList))
                            ]
                             for ixp in range(len(xList))]
                        )
                            for iy in range(len(yList))
                        ]
                        for ix in range(len(xList))
                        ]
                        )
        self.init_fun = interpolate.interp2d(xList, yList, density, kind='linear', fill_value=0.0)

        return density

    def read(self):
        """
        Read out transverse energy from fname.
        """
        self.read_density()

        #print('Calculating N ...\n')
        N = integrate.nquad(lambda x, y: self.init_fun(x, y)[0], [(-9.9, 9.9), (-9.9, 9.9)])
        if self.shiftQ:
            #print('Calculating xbar, ybar ...\n')
            self.xb = integrate.nquad(lambda x, y: self.init_fun(x, y)[0]*x, [(-9.9, 9.9), (-9.9, 9.9)])[0]/N[0]
            self.yb = integrate.nquad(lambda x, y: self.init_fun(x, y)[0]*y, [(-9.9, 9.9), (-9.9, 9.9)])[0]/N[0]
            #print('xbar = {}, ybar = {}.'.format(self.xb, self.yb))
        #print('Calculating R2 ...\n')
        R2 = integrate.nquad(lambda x, y: self.init_fun(x, y)[0]*((x - self.xb)**2 + (y - self.yb)**2), [(-9.9, 9.9), (-9.9, 9.9)])
        self.e0 = N[0]/np.pi
        self.R = np.sqrt(R2[0]/N[0])
        #print('R = {}.'.format(self.R))

    def FT_file(self, r, th):
        """
        Calculate FT from the interpolation function of the initial data taken from self.fname
        :param r: r tilde
        :param th: theta tilde
        :return: FT with normalizaiton FT(0)=2.0/self.t0
        """
        return 2.0*self.R*self.R*self.init_fun(self.xb + r*self.R*np.cos(th),
                                               self.yb + r*self.R*np.sin(th))[0]/(self.e0*self.t0)

    def from_file_theta(self, latt):
        self.theta, self.vz, self.r, self.phi = latt.indices()
        F0 = np.zeros([latt.shape[self.theta], latt.shape[self.vz],
                       latt.shape[self.r], latt.shape[self.phi]], dtype='complex')
        self.read()

        vz = 0

        for th in range(latt.shape[self.theta]):
            for r in range(latt.shape[self.r]):
                for phi in range(latt.shape[self.phi]):
                    F0[th][vz][r][phi] = self.FT_file(
                        latt.lattice[self.r][r], latt.lattice[self.theta][th]
                    ) / (latt.lattice[self.vz][1])

        return F0

    def from_file(self, latt):
        """
        initial condition from file.

        Parameters
        ----------
        latt : object of class Lattice
            Grids points for the initial conditon.
        """
        F0 = self.fft(self.from_file_theta, latt)
        if len(self.dn) == 1:
            if len(self.dn[0]) == 2: # n-perturbation has the priority.
                for n in range(1, len(F0)):
                    if n == self.dn[0][0]:
                        F0[n, :, :, :] = self.dn[0][1]*F0[n, :, :, :]
                    else:
                        F0[n, :, :, :] = self.dos*F0[n, :, :, :]
            elif len(self.dn[0]) == 1:
                for n in range(self.dn[0][0] + 1, len(F0)):
                    F0[n, :, :, :] = self.dos * F0[n, :, :, :]
        return F0

    lookup = {'gaussian': gauss, 'hocky puck': hockypuck, ' Woods-Saxo':ws, 'file':from_file}

    def __init__(self, ICtype, t0=0.1, *harmonics):
        """
        InitCond()

        A InitCond object evalute the initial condition for a KinTran object.

        Parameters
        ----------
        t0 : float
            The initial time.
        dos: float
            The amplitude (delta) of other modes.
        Attributes
        ----------
        smoothQ : bool
            If true, a guassian smooth will apply to the initial transverse profile in coordinate space.
        fname : basestring
            Filename for the initial transverse energy. It is only constrained to Farid's file style.
        R_smooths : float
            the radius squared for the smooth function: smoother
        init_fun : smoothed initial energy density if smoothQ is True. It is, otherwise, the same as input_fun.
        """
        self.ancientQ = False

        #parameters for input initial profile
        self.fname = None
        self.dos = 0.0
        self.smoothQ = False
        self.R_smooth2 = 1.0
        self.init_fun = None

        self.xb = 0.0 # average value of x
        self.yb = 0.0 # average value of y
        self.R = 1.0 # the root mean square radius
        self.shiftQ = True # whether shif the origin or not
        self.e0 = 1.0 # variable used for normalization

        if ICtype in self.lookup:
            self.ICtype = ICtype
            self.t0 = t0
            self.dn = harmonics
            self.theta = self.r = self.phi = self.vz = None
            #print(self.dn)
            #print(len(self.dn))
        else:
            print("""Sorry! So far we have only codes the following types
                  of initial conditions:\n
                  """)
            print(self.lookup.keys(), '\n')
            raise ValueError('Unimplemented initial conditions!')

    def initialize(self, latt):
        """
        F on the grids latt at t = t0.
        """
        self.theta, self.vz, self.r, self.phi = latt.indices()

        return self.lookup[self.ICtype](self, latt)

    def setAncient(self, aQ):
        self.ancientQ = aQ

    def setFile(self, fn):
        self.fname = fn

    def setSmooth(self, sQ):
        """
        evaluate self.smoothQ
        :param sQ: bool
        """
        self.smoothQ = sQ

    def setR_smooth(self, Rs):
        """
        :param Rs: float, smearing radius
        """
        self.R_smooth2 = Rs**2

    def set_dos(self, dos):
        """
        Set the maplitude for the other modes.
        :param dos: float
        """
        self.dos = dos
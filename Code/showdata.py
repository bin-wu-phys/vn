import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate, integrate
import scipy

from Code.kintran import *

class ShowData:
    """
    ShowData()

    A ShowData object collects all the information from the filename.

    Parameters
    ----------
    filename : string
        The filename for loading the data.

    Attributes
    ----------
    kt : a KinTran object
         It can be used to pick all the information of the save KinTran information at a given t.

    vn : ndarray

    """
    def read_vn(self):
        """
        Readout self.v, the v_n's as a function of v[0]=t.
        """
        if self.folder == '':
            file = self.filename + '/' + self.filename + '.vn.npz'
        else:
            file = self.folder + '/' + self.filename + '/' + self.filename + '.vn.npz'

        vn = np.load(file)['vn']
        self.t = vn[:, 0]
        print()
        self.v = np.array([
                            np.concatenate([[vn[t, 0]],
                                            [ vn[t, 2*n+1] + 1.0j*vn[t, 2*n+2]
                                              for n in range((len(vn[t,:])-1)//2)
                                            ]
                                            ]
                                           )
                                            for t in range(len(vn))
                        ])

    def e(self, n, comment=False, check_origin_shifted=False):
        """
        Calculate en from the distribution function.
        If check_origin_shifted == True, for n = 1, calculate <F_1>/<F_0> instead.
        """
        # Read F0 if necessary
        self.read_F0()

        if n == 1:
            if check_origin_shifted:
                pn = 2
            else:
                pn = 3
        else:
            pn = n + 1

        fbg = interpolate.interp1d(self.r, self.F0[0], kind='linear',
                                   fill_value=0.0, bounds_error=False)
        fn = interpolate.interp1d(self.r, self.F0[n], kind='linear',
                                  fill_value=0.0, bounds_error=False)
        ebg = integrate.nquad(lambda r: r ** pn * scipy.real(fbg(r)),
                              [(self.r[0], self.r[-1])]
                              )[0]
        res = integrate.nquad(lambda r: r ** pn * scipy.real(fn(r)),
                              [(self.r[0], self.r[-1])]
                              )[0] + \
              1.0j * integrate.nquad(lambda r: r ** pn * scipy.imag(fn(r)),
                                     [(self.r[0], self.r[-1])]
                                     )[0]

        resnum = np.inner(self.r ** pn, self.F0[n]) / np.inner(self.r ** pn, self.F0[0])


        if comment:
            plt.scatter(self.r, self.r ** pn * self.F0[n].real)
            plt.scatter(self.r, self.r ** pn * self.F0[n].imag)

            rs = np.arange(0, 3.0, 0.05)
            fintep = fn(rs)
            plt.plot(rs, rs ** pn * fintep.real)
            plt.plot(rs, rs ** pn * fintep.imag)

        return resnum, res / ebg

    def set_file(self, filename, folder=''):
        """
        Set self.filename and self.folder.
        :param filename: string
        :param folder: string
        """
        self.filename = filename
        self.folder = folder

        self.load()

    def load(self):
        """
        Load all the info from the data
        """
        self.read_vn()
        self.F0 = None
        self.r = None

    def read_F0(self):
        # Find out the initial distribution
        if self.F0 is None:
            if self.folder == '':
                fd = self.filename
            else:
                fd = self.folder + '/' + self.filename

            files = [f for f in os.listdir(fd) if '.kt.' in f]

            f_init = files[0]
            for i in range(1, len(files)):
                if files[i] < f_init:
                    f_init = files[i]
            f_init = fd + '/' + f_init

            phi, vz = 0, 0
            self.kt.load(f_init)
            self.r = self.kt.kern.ob.latt.lattice[0]
            self.F0 = self.kt.F[:, vz, :, phi]

    def show_F0(self, n, ax=None, lw=4, linestyle='-'):
        """
        Show F0[n] as a function of r.
        :param n: integer
                the harmonic number (in theta).
        :param ax: an Axies object
        """
        # Read F0 if necessary
        self.read_F0()

        returnQ = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
            returnQ = True

        ax.plot(self.r, self.F0[n].real, linestyle, lw=lw,
                    label=r'Re $F_{}$'.format(n)
                )
        ax.plot(self.r, self.F0[n].imag, linestyle, lw=lw,
                label=r'Im $F_{}$'.format(n)
                )

        ax.set_xlabel(r'$r/R$', fontsize=14)
        ax.legend()

        if returnQ:
            return fig, ax

    def __init__(self, filename='', folder=''):
        self.filename = filename
        self.folder = folder
        self.kt = KinTran()

        self.F0 = None
        self.r = None

        if filename != '':
            self.load()
        else:
            self.v = None

    def show_v(self, n, enQ=True, ax=None, lw=4, label='', linestyle='-'):
        """
        """
        returnQ = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
            returnQ = True
        if label == '':
            label = self.filename

        if enQ:
            e2 = self.e(n)[0]
            ax.plot(self.v[:, 0].real, np.abs(self.v[:, n] / e2), linestyle, lw=lw,
                    label=label)
            ax.set_ylabel(r'$v_{}/\epsilon_{}$'.format(n, n), fontsize=14)
        else:
            ax.plot(self.v[:, 0].real, np.abs(self.v[:, n]), linestyle, lw=lw,
                    label=label)
            ax.set_ylabel(r'$v_{}$'.format(n), fontsize=14)

        ax.set_xlabel(r'$\tau/R$', fontsize=14)
        ax.legend()

        if returnQ:
            return fig, ax
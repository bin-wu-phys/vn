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

    def e(self, n, comment=False):
        """
        Calculate en from the distribution function.
        """
        # Find out the initial distribution
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

        self.kt.load(f_init)

        if n == 1:
            pn = 3
        else:
            pn = n + 1
        phi, vz = 0, 0
        fbg = interpolate.interp1d(self.kt.kern.ob.latt.lattice[0], self.kt.F[0, vz, :, phi], kind='linear',
                                   fill_value=0.0, bounds_error=False)
        fn = interpolate.interp1d(self.kt.kern.ob.latt.lattice[0], self.kt.F[n, vz, :, phi], kind='linear',
                                  fill_value=0.0, bounds_error=False)
        ebg = integrate.nquad(lambda r: r ** pn * scipy.real(fbg(r)),
                              [(self.kt.kern.ob.latt.lattice[0][0], self.kt.kern.ob.latt.lattice[0][-1])]
                              )[0]
        res = integrate.nquad(lambda r: r ** pn * scipy.real(fn(r)),
                              [(self.kt.kern.ob.latt.lattice[0][0], self.kt.kern.ob.latt.lattice[0][-1])]
                              )[0] + \
              1.0j * integrate.nquad(lambda r: r ** pn * scipy.imag(fn(r)),
                                     [(self.kt.kern.ob.latt.lattice[0][0], self.kt.kern.ob.latt.lattice[0][-1])]
                                     )[0]

        resnum = np.inner(self.kt.kern.ob.latt.lattice[0] ** pn, self.kt.F[n, vz, :, phi]) / np.inner(
            self.kt.kern.ob.latt.lattice[0] ** pn, self.kt.F[0, vz, :, phi])


        if comment:
            plt.scatter(self.kt.kern.ob.latt.lattice[0],
                        self.kt.kern.ob.latt.lattice[0] ** pn * self.kt.F[n, vz, :, phi].real)
            plt.scatter(self.kt.kern.ob.latt.lattice[0],
                        self.kt.kern.ob.latt.lattice[0] ** pn * self.kt.F[n, vz, :, phi].imag)
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

    def __init__(self, filename='', folder=''):

        self.filename = filename
        self.folder = folder
        self.kt = KinTran()

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
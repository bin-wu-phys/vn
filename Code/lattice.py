#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 12:11:13 2019

@author: biwu
"""
import numpy as np

class Lattice:
    """
    A Lattice object sets up the grid for solving kinetic transport.

    Parameters
    ----------
    shape : tuple of lists or tuples
        an element is of the form:
            (name of coordinate, start, end, number of grid points).
        Here, the coordinate inludes:
                'theta', 'vz tilde', 'r tilde', 'phir tilde'

    Attributes
    ----------
    lattice : list of np arrays
        lists of grid points in theat, r tilde, vz tilde and phir FOR NOW.
    shape : list of ints
        The list of the number of grid points in each dimension.
    """
    #Grids
    
    def theta(rng):
        return 2.0*np.pi*np.arange(rng[0])/rng[0]

    def r_tilde(rng):
        return np.linspace(rng[0], rng[1], rng[2])

    def phir(rng):
        return 2.0*np.pi*np.arange(rng[0])/rng[0]+1.0e-7

    def vz_tilde(rng):
        lnz = np.linspace(rng[0], rng[1], rng[2])
        return (np.exp(lnz) - 1.0)/(np.exp(lnz) + 1.0)

    def __str__(self):
        """
        formulates the output of F.
        """
        return 'Description: {}\nShape: {}'.format(
                self.description, self.shape)

    lookup = {'theta': theta, 'r tilde': r_tilde, 'vz tilde': vz_tilde,
              'phir': phir}

    def __init__(self, *shape):
        self.lattice = []
        self.description = []
        for coor in shape:
            try:
                self.lattice.append(self.lookup[coor[0]](coor[1:]))
                self.description.append(coor[0])
                if coor[0] == 'vz tilde':
                    self.vzmin = coor[1]
                    self.vzmax = coor[2]
            except KeyError:
                print("""Please make sure that the coordinate is one
                      of the following ones:\n""")
                print(self.lookup.keys(), '\n')
                raise KeyError('coorindate ' + coor[0] + ' is not defined!')

        self.shape = [row[-1] for row in shape]

    def indices(self):
        """
        Returns the index of each coordinate in self.lattice.
        """
        theta = self.description.index('theta')
        vz = self.description.index('vz tilde')
        r = self.description.index('r tilde')
        phi = self.description.index('phir')

        return theta, vz, r, phi

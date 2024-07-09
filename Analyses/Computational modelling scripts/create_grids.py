# -*- coding: utf-8 -*-
"""
create reward distribution for the grids
"""

import numpy as np


def bivariate(W, L, l = 2, max_r = 100, min_r = 0):
    """
    Function to generate a 2D reward distribution
    
    Parameters
    ----------
    W : int
        width of the grid.
    L : int
        length of the grid.
    l : float, optional
        the strength of the spatial correlations (given by the kRBF) 
        between the rewards. The default is 2.
        
    max_r : float, optional
        maximum reward. The default is 100.
    min_r : float, optional
        minimum reward. The default is 0.

    Returns
    -------
    numpy array with WxL rewards
    
    """
    kernel_matrix = k_matrix(l, W, L)
   
    bandit = np.random.multivariate_normal([0]*L*W, kernel_matrix)
    bandit = (bandit - min(bandit))/(max(bandit)-min(bandit)) #to normalize
    bandit = bandit * max_r #to let it vary between 0 and max_r
    bandit = np.array([value + min_r for value in bandit])
   
        
    return bandit


def k_matrix(l, W, L):
    """
    Function to generate a correlation matrix
    Needed to generate mulitvariate normal correlated rewards
    
    Parameters
    ----------
    l : float
        the strength of the spatial correlations (given by the kRBF)
        between the rewards.
    W : int
        width of the grid.
    L : int
        length of the grid.

    Returns
    -------
    numpy array with WxL x WxL correlations 
    one correlation for every cell-cell combination

    """
    spatial_vector_list = np.array([np.array([i, j]) for i in range(L) for j in range(W)])
   
    return [[kRBF(spatial_vector_list[i], spatial_vector_list[j], l) 
             for i in range(0, W*L)] for j in range(0, W*L)]


def kRBF(v1, v2, l):
    """
    the Radial Basis Function kernel
    this kernel translates the spatial distances
    to correlations with smoothness lam
    Parameters
    ----------
    v1 : 2D numpy array of ints
        the coordinates of one cell.
    v2 : 2D numpy array of ints
        the coordinates of one cell.
    l : float
        the strength of the spatial correlations (given by the kRBF)
        between the rewards.

    Returns
    -------
    the correlation between two cells,
    based on their euclidian distance
    and the smoothness

    """
    
    return np.exp(-(np.linalg.norm(v1-v2))**2/(2*l**2))


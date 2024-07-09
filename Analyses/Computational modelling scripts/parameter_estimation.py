# -*- coding: utf-8 -*-
"""
script to estimate the model parameters of one (behavioural/computational) participant
"""
import math
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
import scipy
from scipy.optimize import minimize


"""
Prior with which we want to estimate the model parameters:
"""
#prior values of model parameters, used to start the search and to bias the optimization
x0 = [10, 0.5, 0.01]
#std of the priors on the model parameters, used to bias the optimization
xstd = [8, 0.2, 0.2] #actually is this the variance, just wrong name but correct values and use
#covariance matrix of the priors on the model parameters
xcov = np.diag(xstd)

bounds = np.exp(np.array([(-5, 6), (-5, 3), (-5, 3)]))
localized = True #to change if we want to use the localized model or not

def estimate(W, L, nr_trials, nr_blocks, data):
    """
    estimates the model parameters of the data
    using leave-one-out method where one round is left out
    for cross validation
    
    
    then use the median of the estimations as the final result
    and the quartile deviation as a measurement of variance
    
    Parameters
    ----------
    W : int
        width of the grid.
    L : int
        length of the grid.
    nr_trials : int
        nr of trials per round.
    nr_blocks : int
        nr of blocks per experiment.
    data : dataframe
        all the participant's data.

    Returns
    -------
    the estimated parameters and the NLL belonging to that estimation.

    """
    
    #make usefull lists of the data
    datalast = data.query('trial_nr == {}'.format(nr_trials-1))  #this is a dataframe that only contains the data of the last trial of every grid
    #lists with info per grid
    initial_opened = [value for value in datalast.initial_opened]
    
    l_fit_est_list = [0]*nr_blocks
    beta_est_list = [0]*nr_blocks
    tau_est_list = [0]*nr_blocks
    NLL_list = [0]*nr_blocks
    
    for out_index in range(nr_blocks):
        data_partial = data.query('block_nr != {}'.format(out_index))
        args = (W, L, nr_trials, initial_opened[:out_index] + initial_opened[out_index+1:],
                [value for value in data_partial.selected_choice],
                [value for value in data_partial.reward],
                [value for value in data_partial.average_reward])
        
        res = minimize(fun=wrapper, x0 = np.log(x0), args=args, method='SLSQP', bounds=np.log(np.array(bounds)))
        (l_fit_est, beta_est, tau_est) = np.exp(res.x)
        
        l_fit_est_list[out_index] = l_fit_est
        beta_est_list[out_index] = beta_est
        tau_est_list[out_index] = tau_est
        
        #cross validation
        data_partial = data.query('block_nr == {}'.format(out_index))
        cross_val = NLL(W, L, nr_trials, l_fit_est, beta_est, tau_est, 
                        [initial_opened[out_index]], 
                        [value for value in data_partial.selected_choice],
                        [value for value in data_partial.reward], 
                        [value for value in data_partial.average_reward])
        
        NLL_list[out_index] = cross_val
        
    print("we estimated a median of:")    
    print(r'l_fit = %.3f +- %.3f, beta = %.3f +- %.3f, tau = %.3f +- %.3f'
          % (np.median(l_fit_est_list), (np.percentile(l_fit_est_list, 75) - np.percentile(l_fit_est_list, 25))/2, 
             np.median(beta_est_list), (np.percentile(beta_est_list, 75) - np.percentile(beta_est_list, 25))/2, 
                  np.median(tau_est_list), (np.percentile(tau_est_list, 75) - np.percentile(tau_est_list, 25))/2 ))
    print("NLL = %.3f\n" % np.sum(NLL_list))
    print("\n")
    resx = (np.median(l_fit_est_list), np.median(beta_est_list), np.median(tau_est_list))
    resfun = np.sum(NLL_list)
    
        
    return (resx, resfun)


def wrapper(par, *args):
    """
    this function is just a wrapper function that brings the NLL function into the correct format
    to be called by the optimization function
    
    used for the 1 l, 1 beta and 1 tau per data set
    
    Parameters
    ----------
    par : array of floats
        the array of the model parameters that need to be optimized
    *args : tuple of ...
        all the other arguments needed for NLL

    Returns
    -------
    returns the function NLL
    """
    l_fit = np.exp(par[0])
    if (l_fit == 0):
        l_fit = 10e-8
    
    beta = np.exp(par[1])
    if (beta == 0):
        beta = 10e-8
    
    tau = np.exp(par[-1])
    if (tau == 0):
        tau = 10e-8    

    (W, L, nr_trials, initial_opened, selected_choice, reward, average_reward) = args
    
    '''
    The probability that the data is generated 
    by a model with the given model parameters
    is 
    prob = L * bias
    so 
    NLprob = NLL + NLbias = NLL - Lbias
    '''
    term1 = NLL(W, L, nr_trials, l_fit, beta, tau, initial_opened, selected_choice, reward, average_reward)
    term2 = NLbias((l_fit, beta, tau), x0, xcov)
    NLprob = term1 + term2
    return NLprob

    

'''
functions calculating the probability of each choice given a parameter set
'''
def NLL(W, L, nr_trials, l_fit, beta, tau, initial_opened, selected_choice, reward, average_reward):
    """
    NLL (par, observations)
    This function calculates the negative log likelihood of a set of observations (all the choices) given a parameter set

    Parameters
    ----------
    l_fit : float
        model parameter defining the generalization strenght
    beta : float
        model parameter defining the exploration bonus
    tau : float
        model parameter defining the softmax temperature
    
    initial_opened : list of ints
        per grid, the first opened cell
    selected_choice : list of ints
        all chosen cells
    reward : list of floats
        the observed rewards per chosen cell
    average_reward : list of floats
        the average reward per choice, needed to calculate the reward of the initially opened cell

    Returns
    -------
    returns a float, this will be given to the optimizer

    """
    LL = 0
    
    for round_nr in range(0, len(initial_opened)):
        opened_cells = [initial_opened[round_nr]]
        #we need to change the format to a list of coordinates
        opened_cells2D = [[math.floor(value/W), value%W] for value in opened_cells]
        
        first_observation = 2*average_reward[round_nr*nr_trials] - reward[round_nr*nr_trials]
        observations = [first_observation]
        
        for trial_nr in range(0, nr_trials):
            choice = selected_choice[round_nr*nr_trials + trial_nr]
            LL += log_probability(W, L, choice, opened_cells2D, observations, l_fit, beta, tau)
            
            #update the observations before moving on to the next 
            opened_cells.append(choice)
            opened_cells2D = [[math.floor(value/W), value%W] for value in opened_cells]
            observations.append(reward[round_nr*nr_trials + trial_nr])
       
    return -LL


'''
functions needed to calculate the NLL
'''
#more stable if we use the log_softmax function of scipy instead of calculing P our selves and then logging it
def log_probability(W, L, choice, opened_cells2D, observations, l_fit, beta, tau):
    m, s = GP(observations, W, L, l_fit, opened_cells2D)
    UCB = [m[i] + beta * s[i] for i in range(0, W*L)]
    if localized == True:
        #UCB gets weighted
        prior_choice = [(tile_number-(tile_number%W-1))//W, tile_number%W]
        UCB = localizer(UCB, prior_choice, W, L)
    UCBtau = [value/tau for value in UCB]
    log_P = scipy.special.log_softmax(UCBtau)
    
    return log_P[int(choice)]
    

def GP(observations, W, L, l_fit, opened_cells2D, noise=True, epsilon = 0.0001):   
    '''
    This function starts from the observations
    since the GP assumes that not observed states have as default a value of 0,
    we rescale the reward to have a mean 0 and vary between max bounds of -0.5, 0.5
    
    Parameters
    ----------
    observations : the observations
    W : width of the bandit.
    L : lenght of the bandit.
    l_fit: the generalization strength with which the participant smooths out the observed rewards
    opened_cells2D : list of which cells are opened.
    noise: True if the participant assumes noise in their observations
    epsilon: the assumed noise, 
            measured as the variance of the reward around the mean
            std is 1 up 100, 
            var is 0.01**2 = 0.0001
             
    
    Returns
    -------
    Returns the mean function m(x) and the uncertainty function s(x) per tile.
    '''
    
    cells =[[i,j] for i in range(0, W) for j in range(0, L)]    
    
    kernel = RBF(l_fit, "fixed")
    if noise:
        kernel += WhiteKernel(epsilon, "fixed")
        
    gp = GaussianProcessRegressor(kernel=kernel)
    gp.fit(opened_cells2D, [(reward-40)/90 for reward in observations]) 
    
    mlist, sigmalist = gp.predict(cells, return_std=True)
    
    return mlist, sigmalist



def NLbias(par, x0, xcov):
    """
    function to calculate the negative log of the bias term
    this term comes from a multivariate normal distribution
    (the prior of the model parameters)

    Parameters
    ----------
    par: (l_fit, beta, tau)
    x0 : mean of the prior
    xcov : covariance matrix of the prior

    Returns
    -------
    return the NL of the bias term

    """
    NLbias = 0
    for i in range(len(par)):
        NLbias += (par[i]-x0[i])**2/xcov[i][i]
   
    return NLbias/2

def localizer(UCB, prior_choice, W, L):
    cells =[[i,j] for i in range(0, W) for j in range(0, L)]  
    IMD =  [1/(abs(x-prior_choice[0]) + abs(y-prior_choice[1])) if [x,y]!=prior_choice else 1 for [x,y] in cells]
    UCBloc = [IMD[i]*UCB[i] for i in range(W*L)]
    
    return UCBloc
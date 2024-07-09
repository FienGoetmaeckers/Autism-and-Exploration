# -*- coding: utf-8 -*-
"""
parameter recovery code
"""
from csv import DictWriter
import pandas as pd
import math
import numpy as np
import random
from parameter_estimation import estimate, localizer
from solving_models import GP, softmax
from bandits11_l2 import bandits

#bounds defined by Tuckeys Fence on our data estimates
bounds_lfit =  (0.006737946999085467, 28.1967)
bounds_beta =  (0.006737946999085467, 1.01493)
bounds_tau =  (0.006737946999085467, 0.11889)


#experiment settings to simulate
nr_trials = 25
nr_blocks = 10
W = L = 11
l = 2
localized = True

    
"""
step 1: generate data with given parameters
to do this, we follow exactly the same steps 
as in the Learnsolver function of the solving_models script
but simultaneously, we write out data in the format of behavioural data
"""
#we need to make a dataframe similar to the behavioural data
data = pd.DataFrame(columns=["block_nr", "trial_nr", "initial_opened", "selected_choice", "reward", "average_reward"])
    
#sample model parameters for the agent
beta = random.uniform(bounds_beta[0], bounds_beta[1])
tau = random.uniform(bounds_tau[0], bounds_tau[1])
l_fit = random.uniform(bounds_lfit[0], bounds_lfit[1])

for block_nr in range(0, nr_blocks):
    bandit = bandits[block_nr]  
    hidden = np.array([True]*W*L) 
    observed_bandit = [[] for value in bandit]
        
    '''
    step 1: one random cells is revealed at the begin of the experiment
    '''
    tile_number = random.randint(0, W*L-1)
    reward = random.normalvariate(bandit[tile_number], 1)
	    
    #save that this cell has been opened before and save the history in observed_bandit
    hidden[tile_number] = False #save that this cell has been opened
    observed_bandit[tile_number].append(reward)
    
    #save for dataframe
    initial_opened = tile_number
    rewardlist = [reward] #to calculate the average accumulated reward per trial
       
        
    for trial_nr in range(0, nr_trials):
        '''
        per choice to make, 
        the first step is to learn from the prior rewards
        and make predictions about the cells of the grid
        this is done via Gaussian Process regression
        '''
        m, s = GP(observed_bandit, W, L, l_fit, hidden)
         
        '''
        to translate the expectations to a probability to select a cell,
        we use UCB and a softmax rule
        '''
        UCB = [m[i] + beta * s[i] for i in range(0, W*L)]
        if localized == True:
            #UCB gets weighted
            prior_choice = [(tile_number-(tile_number%W-1))//W, tile_number%W]
            UCB = localizer(UCB, prior_choice, W, L)


        P = softmax(UCB, W, L, tau)
        
        #the agent will choose from the tiles, for which the probabilites are given by P 
        tile_number = random.choices(np.arange(0, W*L), weights=P)[0]
        reward = random.normalvariate(bandit[tile_number], 1)
         
        #save that this cell has been opened before and save the history in observed_bandit
        hidden[tile_number] = False #save that this cell has been opened
        observed_bandit[tile_number].append(reward)
        
        rewardlist.append(reward)
           
        result_trial = {"block_nr": block_nr, "trial_nr": trial_nr, "initial_opened": initial_opened, "selected_choice": tile_number, "reward": reward, "average_reward": np.mean(rewardlist)} 
        data = data.append(result_trial, ignore_index=True)
 
    
"""
step 2: parameter estimation
"""
est = estimate(W, L, nr_trials, nr_blocks, data)
       
"""
step 3: save generating and estimated parameter combo
"""
results = {"l_fit_g": l_fit, "l_fit_e": est[0][0], "beta_g": beta, "beta_e": est[0][1], "tau_g": tau, "tau_e": est[0][2], "optimal NLL": est[1]}

file_name = 'estimated_parameters_{}x{}_{}_l{}.csv'.format(W,L, nr_trials, l)
field_names=["l_fit_g", "l_fit_e", "beta_g", "beta_e", "tau_g", "tau_e", "optimal NLL"]
#open CSV file in append mode
with open(file_name, 'a') as f_object:
	dictwriter_object = DictWriter(f_object, fieldnames=field_names)
	dictwriter_object.writerow(results)
	f_object.close()
              
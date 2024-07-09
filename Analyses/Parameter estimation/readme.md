## The Gaussian Process Regression - Upper Confidence Bound model
* **solving_models.py**  
    contains the computational steps for the GP-UCB model. Also, the localized version of the model is included here.

## Estimating model parameters
* **parameter_estimation_behavioural.py**  
    is the main script to estimate the three model parameters of one participant. It reads in the data of the one participant (participant nr p_index, with p_index a real number between 0 and 660 of choice), calls the estimate function from parameter_estimation.py and writes the estimated model parameters, the out-of-sample NLL and the AIC as a new line in csv file.
* **parameter_estimation.py**  
    contains the estimate function, called upon in parameter_estimation_behavioural.py, to estimate the model parameters of one block of data (one participant). The estimate function uses leave-one-*(round-)*out cross-validation to estimate model parameters iteratively over 9 out of 10 rounds, and cross-validates the estimation on the left-out round to assess model fits.


## Creating reward distributions
* **create_grids.py**  
    used to create the reward distributions used in the behavioural task. 

    

"""
script to estimate the model parameters of one participant
"""
import sys
from csv import DictWriter
import pandas as pd
from parameter_estimation import estimate
import os


#specific variables for this data file
nr_blocks = 10 
nr_trials = 25
nr_participants = 660
W = L = 11
#name of the file to read the data from
data_name = "../../Data/data"

"""
read in the data
"""
data = pd.read_csv(data_name + '.csv', delimiter=',')
#select one participant to estimate in this script

p_index = 0 #change for next participant estimation

participant = data.subjectID.unique()[p_index]
print("For participant {}".format(participant))
data_p = data.query('subjectID == "{}"'.format(str(participant)))

"""
estimate the model parameters of this participant
"""
est = estimate(W, L, nr_trials, nr_blocks, data_p)

"""
save the output
"""

condition = 0

resultsM2      = {"Participant": data_p["subjectID"].values[0], "l_fit": est[0][0], "beta": est[0][1], 
	     	     "tau": est[0][2], "condition": condition, "NLL": est[1], "AIC": 2*5 + 2 * est[1]}

with open("estimates.csv", 'a') as f_object:
    field_names = ["Participant", "l_fit", "beta", "tau", "condition", "NLL", "AIC"]
    dictwriter_object = DictWriter(f_object, fieldnames=field_names)
    dictwriter_object.writerow(resultsM2)
    f_object.close()


del data_p

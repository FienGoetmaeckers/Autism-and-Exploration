# Autism-and-Exploration
All scripts for the project "Different exploration strategies along the autism spectrum: Diverging effects of autism diagnosis and autism traits"
* link to preprint https://osf.io/preprints/psyarxiv/3mvxk 
* preregistration https://osf.io/mz3h5 

## Behavioural experiment
All needed code to run the online behavioural experiment described in the paper
* MAB.html will run the behavioural experiment
* grid_file.js contains reward distributions
  
The experiment's code is written in JavaScript and uses the jsPsych framework (https://www.jspsych.org/latest/). All needed scripts, packages and plugins are included in the *Behavioural experiment* folder.

Also, if you want to try out the truffle task, you can access it online here: https://braemlab.ugent.be/Fien/Correlation_N/MAB.html 

This link serves as a demo to understand the mechanisms of the behavioural task. This version does not save your data and can be run on Chrome internet browsers on laptops and desktops. Completing the full task should take about 15 minutes. 


## Analysis
### Model and parameter estimation
Contains all scripts to read in behavioural data and estimate the model parameters per participant. These scripts use Python programming language. We used Python version 3.9.


### Statistical analyses
* Analysis_script.R : R script for conducting statistical analyses on behavioural and modelling data, includes code for creation of figures. We used R version 4.4.1. 

## Data
Includes all needed data to replicate findings. 

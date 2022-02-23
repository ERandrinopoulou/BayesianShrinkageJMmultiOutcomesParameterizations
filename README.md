# BayesianShrinkageJMmultiOutcomesParameterizations
## Bayesian shrinkage approach for a joint model of longitudinal and survival outcomes assuming different association structures

This repository includes the code for a Bayesian shrinkage approach that models multiple longitudinal (2 continuous and 1 binary) outcomes together with a time-to-event assuming different association structures (value, slope and area). 

Specifically:
* "**ModelJAGS**": includes the joint model for jags.
* "**PrepareData**": includes the preparation of the data (pbc2 from the survival package).
* "**Fit**": includes the main code. Specifically, it loads the packages and the data (pbc2 from the survival package), runs the functions, creates the model, runs the model in JAGS and saves the results.
* All aforementioned files use the pbc2 (from the survival package) data to run the model. The folder "**Simulate data**" consist of code on how to simulate data sets from the joint model assuming the value, slope and area parameterization.

How does it work:
* Download all files and place them in one folder.
* Set as working directory in R this folder.
* Run the code in "Fit" for fitting the joint model.


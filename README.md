# BayesianShrinkageJMmultiOutcomesParameterizations
## Bayesian shrinkage approach for a joint model of longitudinal and survival outcomes assuming different association structures

This repository includes the code for a Bayesian shinkage approach that models multiple longitudinal (2 continuous and 1 binary) outcomes together with a time-to-event assuming different association structures. 

Specifically:
* "**ModelJAGS**": includes the joint model for jags.
* "**PrepareData**": includes the preparation of the data.
* "**Fit**": includes the main code. Specifically, it loads the packages and the data, runs the functions, creates the model, runs the model in JAGS and saves the results.

How does it work:
* Download all files and place them in one folder.
* Set as working directory in R this folder.
* Run the code in "Fit" for fitting the joint model.


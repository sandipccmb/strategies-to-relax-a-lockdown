# strategies-to-relax-a-lockdown
Combining serology with case-detection, to allow the easing of restrictions against SARS-CoV-2: a modelling-based study in India 

Sandip Mandal1, Hemanshu Das2, Sarang Deo2, Nimalan Arinaminpathy3*

1 Independent Consultant, New Delhi, India; 2 Indian School of Business, Hyderabad, India; 3 MRC Centre for Global Infectious Disease Analysis, School of Public Health, Imperial College London, London, UK

Copyright (C) <2020>, Sandip Mandal et. al. All rights reserved. Released under the GNU General Public License (GPL)

This repository contains codes and data used to simulate and analyze SARS-CoV-2 transmission to examine how seroprevalence data could guide a test-and-isolate strategy, for lifting a lockdown.
The model code is written in MATLAB and results are saved as MATLAB data files (extension .mat) and to make plots the following MATLAB codes are used:
1. Simulate_Fig1.m
2. Simulate_Fig2.m
3. Simulate_Fig3.m
4. Simulate_Fig4.m

OS requirements
The codes developed here are tested on Windows and Mac OSX operating system. 

Codes & their functionality
1.	Setup_model_India2.m is the main source file where all the variables are defined, input parameters values are assigned and ‘transmission rate’ is estimated corresponding to a given R0 value. 
2.	get_addresses.m is a function to define the variables in structural form. Here all variables are structured in three age groups. 
3.	make_model4.m is a function which captures the transition between different variables and gives linear and non-linear parts of the model separately.
4.	goveqs_basis3.m This captures the full model equations and returns state of each compartment at next step.
5.	get_objective.m gives the output of the model in Pre-lockdown and post-lockdown scenarios.
6.	goveqs_scaleup.m this function scaleup the model parameters linearly with assigned time period. 
7.	find_R0.m calculates reproduction number R0 for given transmission parameter 'beta' and other model parameters.
8.	linspecer.m This function can be used to plot lots of lines with distinguishable and nice looking colors.
9.	Make_dots.m this is used to create the dots needed for Figure 3, 4.
10.	Simulate_Fig1.m Plot of figure1 in the paper.
11.	Simulate_Fig2.m Plot of figure2 in the paper.
12.	Simulate_Fig3.m Plot of figure3 in the paper.
13.	Simulate_Fig4.m Plot of figure4 in the paper.

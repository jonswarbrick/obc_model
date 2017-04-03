## Synopsis

Code associated to paper 'Endogenous financial crises in a model of occasionally binding constrained borrowing' to calibrate, simulate time series and impulse response functions, report moments and plot results.

## Details of Files and Models

For solving and simulating, printing results and plotting figures:
- run.m - solve and simulate model(s)
- multiple_run.m - solve and simulate model(s) looping over range of settings/parametrisations
- print_moments.m - compute and display simulated moments
- paperMomentsTable.m - display table of simulated moments for paper
- plot_IRF.m - plot IRF simulations
- multipleFigurePlot_IRF.m - plot IRFs for large number of models/setups
- plot_SIM.m - plot simulated time series

For calibration:
- calibration.m - to calibrate multiple models (calling models/calibrate_parameters.m)
- models/calibrate_parameters.m - to set up calibration (which parameters, intial values, step size, tolerance etc)

Main model files in folder 'model':
models/rbc.mod
models/gk.mod
models/obc.mod
models/nk.mod
models/nkobc.mod

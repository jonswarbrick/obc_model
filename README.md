## Synopsis

Code associated to paper 'Credit crunches from occasionally binding bank borrowing constraints' to calibrate, simulate time series and impulse response functions, report moments and plot results.

## Details of Files and Models

For solving and simulating, printing results and plotting figures:
- run_simulation.m - solve and simulate model(s)
- analysis.m - compute and display simulated moments
- plots.m - plot IRF simulations

For calibration:
- calibration.m - to calibrate multiple models (calling models/calibrate_parameters.m)
- models/calibrate_parameters.m - to set up calibration (which parameters, intial values, step size, tolerance etc)

Main model files in folder 'model':

Data:
- fetch_data.m - fetches data on real economy from FRED
- run_moments.m - report empirical moments for paper and plot equity issuance time series
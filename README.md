## Synopsis

Code associated to paper 'Credit crunches from occasionally binding bank borrowing constraints' to simulate time series and impulse response functions, report moments and plot results.

## Details of Files and Models

For solving and simulating, printing results and plotting figures:
- run_simulation.m - solve and simulate model(s), using files in directory 'model'
- analysis.m - compute and display simulated moments
- plots.m - print paper figures of empirical and model data

For calibration:
- calibration.m - calibrate model shock parameters

Data:
- data/fetch_data.m - fetches data on real economy from FRED
- data/run_moments.m - report empirical moments for paper
-
## Additional software

In combination with the supplied code, the following software was used to generate the results in the paper:
- Matlab R2018b (win)
- Dynare 4.5.6
- dynareOBC (https://github.com/tholden/dynareOBC version available 28/11/18) and dependencies
- Gurobi 8.1.0

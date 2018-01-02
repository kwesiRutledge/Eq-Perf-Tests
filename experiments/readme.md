# Experiments

The experiments that were performed, related to Equalized Recovery, are contained in this government.

## Usage

For example, to run Experiment Set #13 (observer_comparison13.m) , while defining the 'verbosity','time_horizon','M1', and 'M2'.

```
results = observer_tests(13,{'verbosity,time_horizon,M1,M2'})
```

## Experiment List

13. The objective of this experiment is to design a finite horizon, affine estimator (fhae) that is robust against the possibility of 1 observation missing in the entire sequence of length T.

14. The objective of this set of experiments was to attempt to use the new robust optimization method to synthesize a finite horizon, affine estimator (assuming that all data is available or at most 1 piece of data missing during the time horizon). Unfortunately, for the default experiment set, a feasible estimator exists ONLY if we exclude the possibility of missing observations at t=4 and t=5.

15. The objective of this set of experiments are:
	* Consider the effect of changing the constant M1 on the minimal M2 that can be found for Equalized Recovery (assumes that the time horizon T is 6, by default). Note: A line search could be implemented to find minimal M1 AND minimial M2 using this method, but it is not considered here.
	* Simulate when an estimator is defined for one time horizon, T, but the controller is run for an arbirary time horizon, T2.
	


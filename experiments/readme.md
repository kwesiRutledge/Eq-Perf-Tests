# Experiments

The experiments that were performed, related to Equalized Recovery, are contained in this government.

## Usage

For example, to run Experiment Set #13 (observer_comparison13.m) , while defining the 'verbosity','time_horizon','M1', and 'M2'.

```
results = observer_tests(13,{'verbosity,time_horizon,M1,M2'})
```

## Experiment List

6. The objective of this experiment was to observe how small the norm ||xi(T)||{\infty} could become for the worst case process and measurement noise for various choices of time horizon T. In this case, we considered the choice of using a full state or reduced order observer/estimator and compared the results in a plot.

7.
8.
9.
10.
11.
12.

13. The objective of this experiment is to design a finite horizon, affine estimator (fhae) that is robust against the possibility of 1 observation missing in the entire sequence of length T.

14. The objective of this set of experiments was to attempt to use the new robust optimization method to synthesize a finite horizon, affine estimator (assuming that all data is available or at most 1 piece of data missing during the time horizon). Unfortunately, for the default experiment set, a feasible estimator exists ONLY if we exclude the possibility of missing observations at t=4 and t=5.

15. The objective of this set of experiments are:
	* Consider the effect of changing the constant M1 on the minimal M2 that can be found for Equalized Recovery (assumes that the time horizon T is 6, by default). Note: A line search could be implemented to find minimal M1 AND minimial M2 using this method, but it is not considered here.
	* Simulate when an estimator is defined for one time horizon, T, but the controller is run for an arbirary time horizon, T2.

16. The objective of this set of experiments was to duplicate the results of Experiment 14, when the proper 'E' matrix was used. Started before experiment 15 was completed because an error (no consideration of the E matrix) was observed in the design of experiments 14 and 15.

17. Reran the counterexample provided in the paper. Here we compare the correctness of multiple synthesis methods, where the goal is to correctly identify when an estimator (aka gains L or (F,u0)) is feasible. We defined feasibility as "feasible" if an optimal objective was less than M and "infeasible" otherwise. The 3 objectives considered were:
	* ||A+LC||{infty} M + ||L||{\infty} \eta_v + ||E||{\infty} \eta_w
	* ||(A+LC)\xi(0)||{\infty} + ||Lv(0)||{\infty} + ||Ew(0)||{\infty}
	* ||(A+FC)\xi(0) + Fv(0) + u_0 + Ew(0)||{\infty}

18. Created a system that could achieve equalized recovery with a condition of "any 1 observation in the time horizon T" can be missing.

19. Tests a function that generates fhae.
	* Repeats the tests from 13,14,16
	* function should be robust enough to allow for E/B_w and G/C_v matrices

20. Tests the ROO solution to the ACC problem to see if it can achieve equalized performance as well as plot how it might achieve equalized recovery. 

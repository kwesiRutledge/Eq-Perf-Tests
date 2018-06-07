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
	* Figure 1 of this script became Figure 1 in the paper.

17. Reran the counterexample provided in the paper. Here we compare the correctness of multiple synthesis methods, where the goal is to correctly identify when an estimator (aka gains L or (F,u0)) is feasible. We defined feasibility as "feasible" if an optimal objective was less than M and "infeasible" otherwise. The 3 objectives considered were:
	* ||A+LC||{infty} M + ||L||{\infty} \eta_v + ||E||{\infty} \eta_w
	* ||(A+LC)\xi(0)||{\infty} + ||Lv(0)||{\infty} + ||Ew(0)||{\infty}
	* ||(A+FC)\xi(0) + Fv(0) + u_0 + Ew(0)||{\infty}

18. Created a system that could achieve equalized recovery with a condition of "any 1 observation in the time horizon T" can be missing.

19. Tests a function that generates fhae.
	* Repeats the tests from 13,14,16
	* function should be robust enough to allow for E/B_w and G/C_v matrices

20. Tests the ROO solution to the ACC problem to see if it can achieve equalized performance as well as plot how it might achieve equalized recovery. 
21.
22.
23. Can we robustify against languages with more than one word?
	* The answer is that our math, as of April 12th did not correctly handle more than one word.
	* In this experiment, I compare one of our earlier experiments (#15) with the performance of a "worst case" language that we have proven to be correctly handled in the current framework.
	* Plots compare the use of our old method for robustifying to languages with more than one word with a "worst case" language design. Plots 3 and 4 specifically display each estimators error when the following words are applied:
		* 101111
		* 110111
		* 111011
		* 111101

24. Recreating the results of ROO ACC with "worst case" language L-star
	* Figure 3 in the ADHS 2018 Paper comes from this experiment.
25. Solving the language problem with the gain switching controller concept. Proof of concept. Still need theory for this.
26. Observing the effect of underconstraining the optimization variable Q when data is missing.
	* When data is missing, the C-bar matrix contains zeros in certain rows which (when placed into the prodcuts that form our constraints) should lead to an under constraining of Q.
	* Formerly, the Q was constrained because we did not force v(t) = 0 when data was missing.
	* Figure 2 in the ADHS 2018 paper, comes from this experiment.
27.
28.	
29. Set of experiments that are required to implement a FHAE_pb design change where L is saved locally as a cell matrix.
	* Introduced Free Recovery Problem and Function.
	* Testing Language Generation for Filter Synthesis.
30. Set of experiments for synthesizing filters that are "modular" in a sense.

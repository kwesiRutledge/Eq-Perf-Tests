# Eq-Perf-Tests

## Description

This repository contains the code used to generate the results of:
* [Optimization-Based Design of Bounded-Error Estimators Robust to Missing Data](https://doi.org/10.1016/j.ifacol.2018.08.027)
* [Prefix-based Bounded-error Estimation with Intermittent Observations](https://doi.org/10.23919/ACC.2019.8814707)
* [Finite horizon constrained control and bounded-error estimation in the presence of missing data](https://doi.org/10.1016/j.nahs.2020.100854)
* [Belief-prefix Control for Autonomously Dodging Switching Disturbances](https://doi.org/10.23919/ECC51009.2020.9143990)

Most of the reproducible code should be found in the 'results' directory.

## Reproducing Results

### HSCC 2021

(In progress.)

## Requirements

Note that this repository also makes use of the following external libraries:
* [Multi-Parametric Toolbox 3](http://mpt3.org) (which utilizes [YALMIP](http://yalmip.github.io) and optimization solvers such as [gurobi](http://gurobi.com))
* [Mosek](http://mosek.com) (An optimization library which should be optional.)
* [Gurobi](http://gurobi.com) (This is a commercial optimization solver.)

One source of errors is that one of the above required libraries is not installed.

## Usage

By including the contents of the 'functions' directory, one should have access to the functions necessary to reproduce the results.

For an introduction into how to use the code, refer to [this repository's wiki](https://github.com/kwesiRutledge/Eq-Perf-Tests/wiki). If you have any questions/would like for something to be added to this list, then please let me know.

## FAQ

### I've previously installed the above libraries, but when I try to run the code I get an error that says they cannot be found (i.e. sdpvar() is not recognized or Polyhedron() is not recognized). How do I fix this?

This is usually because YALMIP or MPT3 is not on the path. To add these back onto the path, find the directory where they are saved and from MATLAB you can right click on them and select `Add Selected Folders and Subfolders to Path`. You can also add the directories to the path programmatically using `addpath()` and `genpath()`.

### Whenever I try to run the bilinear optimization solvers with YALMIP, I get an error dealing with array bounds. How do I fix this?

An example of this type of error:
```
Index exceeds array bounds.

Error in yalmip2gurobi (line 220)
                Qi(map(monomials(k),1),map(monomials(k),2)) = Qi(map(monomials(k),1),map(monomials(k),2)) + di(k)/2;

Error in callgurobi (line 12)
model = yalmip2gurobi(interfacedata);

Error in solvesdp (line 368)
    eval(['output = ' solver.call '(interfacedata);']);

Error in optimize (line 31)
[varargout{1:nargout}] = solvesdp(varargin{:});

Error in LCSAS/FindConsistentBeliefController (line 303)
		optim0 = optimize(optimization_constraints,[],ops);
```

This typically is because YALMIP's installer has chosen one of the "stable" versions of YALMIP and not the one that correctly support's Gurobi's bilinear optimization tools. To fix this replace `yalmip2gurobi.m` in your YALMIP installation with the "bleeding edge" version of it [here](https://github.com/yalmip/YALMIP/blob/develop/solvers/yalmip2gurobi.m).
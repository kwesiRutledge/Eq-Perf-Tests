# Eq-Perf-Tests

## Description

This repository contains the code used to generate the results of:
* [Optimization-Based Design of Bounded-Error Estimators Robust to Missing Data](https://doi.org/10.1016/j.ifacol.2018.08.027)
* [Prefix-based Bounded-error Estimation with Intermittent Observations](https://doi.org/10.23919/ACC.2019.8814707)
* [Finite horizon constrained control and bounded-error estimation in the presence of missing data](https://doi.org/10.1016/j.nahs.2020.100854)
* [Belief-prefix Control for Autonomously Dodging Switching Disturbances](https://doi.org/10.23919/ECC51009.2020.9143990)

Most of the reproducible code should be found in the 'results' directory.

## Usage

By including the contents of the 'functions' directory, one should have access to _most_ of the functions necessary to reproduce the results.

Note that this repository also makes use of the following external libraries:
* [Multi-Parametric Toolbox 3](http://mpt3.org) (which utilizes [YALMIP](http://yalmip.github.io) and optimization solvers such as [gurobi](http://gurobi.com))
* [Mosek](http://mosek.com) (An optimization library which should be optional.)
* [Gurobi](http://gurobi.com) (This is a commercial optimization solver.)
* [PENBMI](http://www.penopt.com/penbmi.html) (Optional. This is one of the few bilinear optimization solvers that exist.)

One source of errors is that one of the above required libraries is not installed.

## Warnings

This code is designed to run with all directories automatically added to path.

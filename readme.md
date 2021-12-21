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

## Requirements

Note that this repository also makes use of the following external libraries:
* [Multi-Parametric Toolbox 3](http://mpt3.org) (which utilizes [YALMIP](http://yalmip.github.io) and optimization solvers such as [gurobi](http://gurobi.com))
* [Mosek](http://mosek.com) (An optimization library which should be optional.)
* [Gurobi](http://gurobi.com) (This is a commercial optimization solver.)
* [PENBMI](http://www.penopt.com/penbmi.html) (Optional. This is one of the few bilinear optimization solvers that exist.)

One source of errors is that one of the above required libraries is not installed.

## Usage

By including the contents of the 'functions' directory, one should have access to _most_ of the functions necessary to reproduce the results.

This section will include some small introductions into how to use the code. If you have any questions/would like for something to be added to this list, then please let me know.

### Affine Systems

The fundamental control system used in this library is the discrete-time, affine system. In math, this is a system whose state at time t is $$x_t$$ updates over time according to the equation:

$$x_{t+1} = Ax_t + Bu_t + f$$
where there is also an input $$u_t$$ that the controller or user can provide to impact what the next state will be.
To create a system of that form in MATLAB one would write:
```
A = [ %... your code here ... % ]
B = [ %... your code here ... % ]
f = [ %... your code here ... % ]
system1 = Aff_Dyn(A,B,f,zeros(n))
```
where n is the dimension of the system.

We also consider systems with output feedback:
$$x_{t+1} = Ax_t + Bu_t + f,$$
$$y_{t} = C x_t$$
To create this system, one would write:
```
A = [ %... your code here ... % ]
B = [ %... your code here ... % ]
f = [ %... your code here ... % ]
C = [ %... your code here ... % ]
system2 = Aff_Dyn(A,B,f,C)
```
And in the most general case, we also consider systems where there are also disturbances that impact how the state evolves and how the measurements are created. That is there is a disturbance which impacts the next state $$x_{t+1}$$ and measurements $$y_t$$. Such a system will look like this:
$$x_{t+1} = Ax_t + Bu_t + B_w w_t + f,$$
$$y_{t} = C x_t + C_v v_t$$
where the process disturbance at time t $$w_t$$ and the measurmeent disturbance at time t $$v_t$$. These disturbances are assumed to come from a type of set known as a polytope (i.e. $$w_t \in \mathcal{W}$$ and similarly $$v_t \in \mathcal{V}$$). We use the toolbox [MPT3](http://mpt3.org) to easily create these sets. To see how to create such a discrete-time affine system consider the following example:
```
g = 9.8;
dt = 0.01;
n_x = 2; %Dimension of x

A1 = [ 1 , dt ; 0 , 1 ];
B1 = [ 0 ; dt/(m1) ];
f1 = [ 0 ; -dt * g ];
C1 = eye(n_x)

eta_w = 0.1;
W1 = Polyhedron('lb',-eta_w,'ub',eta_w);
V1 = Polyhedron('lb',-0.05*ones(1,n_x),'ub',0.05*ones(1,n_x));

ad1 = Aff_Dyn(	A1 , B1 , f1 , C1, ...
				W1 , V1 , ...
				B1 , C1 );
```

### POB_Feedback

The result of the paper [Belief-prefix Control for Autonomously Dodging Switching Disturbances](https://doi.org/10.23919/ECC51009.2020.9143990) is a special class of controller that is based on the idea of a BeliefGraph. In order to generate this controller to solve a reachability problem you can call `the solve_reachability_problem()` member function:

```
your_lcsas = LCSAS( %... Your code here ... % )
[fb_law,opt_out] = lcsas.solve_reachability_problem(P_target)
```
Note that normally at least a target region is needed (in the form of a Polyhedron() object from MPT3.0).

## Warnings

There are some functions in this code which manipulate the path only if you are using an Ozay Lab computer. In most cases, you will have to modify your own path in order to use this toolbox. Be prepared to include YALMIP, MPT3, etc. in your path before running these results.

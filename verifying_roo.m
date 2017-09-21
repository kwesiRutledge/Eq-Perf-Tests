% verifying_roo.m
% 	This script verifies that the ACC Observer System

clear all;
close all;
clc;

%% Constants

ki = 0.2890; % already divided by m
f_0 = 0.6558; % already divided by m
f_1 = 0.0058; % already divided by m
T = 1/25; % controller period

%% Continuous matrices
A = [-f_1, 0, 0;
     -1, 0, 1;
      0, 0, 0;];
B = [ki;
     0;
     0;];
F = [-f_0;
      0;
      0;];
E = [0;
     0;
     1;];
C = [1, 0, 0;
     0, 1, 0];

%% convert to discret time
A_s = @(s) expm(s*A);
Ad = A_s(T);
Bd = integral(A_s, 0, T, 'ArrayValued', true) * B;
Fd = integral(A_s, 0, T, 'ArrayValued', true) * F;
Ed = integral(A_s, 0, T, 'ArrayValued', true) * E;

% Ad = [Aaa, Aau;
%       Aua, Auu];
Aaa = Ad(1:2,1:2);
Aau = Ad(1:2,3);
Aua = Ad(3,1:2);
Auu = Ad(3,3);

Ba = Bd(1:2,1);
Bu = Bd(3,1);
Fa = Fd(1:2,1);
Fu = Fd(3,1);

Ea = Ed(1:2,1);
Eu = Ed(3,1);

%% Optimize to find the optimal L

% System Disturbance Constraints
d = 0.1;
m = 0.05;
perf_level = 1;

%Define Optimization Variables
L = sdpvar(size(Auu,1),size(Aau,1),'full');
err0 = sdpvar(size(Ad,1),1,'full');

delta_t0 = sdpvar(1,1,'full'); 
mu_t0 = sdpvar(size(C,1),1,'full');
mu_t1 = sdpvar(size(C,1),1,'full');

alpha1 = sdpvar(1,1,'full');

% Define Error in the next step for our estimate of lead car velocity
next_step_err = norm( (Auu - L*Aau) * err0 + Eu * delta_t0 - (Aua - L*Aaa)*mu_t0 - L*Ea*delta_t0 - L*mu_t1,Inf);

% Define Constraints
constrs = [ -perf_level <= err0 <= perf_level , uncertain(err0) ];
constrs = constrs + [ -m <= [ mu_t0 ; mu_t1 ] <= m , uncertain(mu_t0) , uncertain(mu_t1) ];
constrs = constrs + [ -d <= [ delta_t0 ] <= d, uncertain(delta_t0) ];
constrs = constrs + [ next_step_err <= alpha1 ];

%Optimize
sol = optimize(constrs,alpha1);

disp(' ')
disp(['The optimal error for the ACC System with ROO is: ' num2str(value(alpha1))])




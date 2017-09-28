% verifying_roo.m
% 	This script verifies that the ACC Observer System

clear all;
close all;
clc;

%% Constants

addpath('./OzayTest1/')

con = constants_normal; 

rmpath('./OzayTest1/')

kappa = con.f1_bar/con.mass;
ekt = exp(-kappa*con.dT);

A = [ ekt 				0 		0 ; 
      (ekt-1)/kappa 	1 		con.dT;
      0 				0 		1]; 
B = (1/(con.f1_bar^2))*[ con.f1_bar*(1-ekt) ;
                         con.mass*(1-ekt) - con.dT*con.f1_bar ;
                         0 ];
B_cond_number = max(abs(B));
B = B/B_cond_number;

E = [0 ; con.dT^2/2; con.dT];
K = [ (con.f0_bar/con.f1_bar)*(ekt-1) ;
          (con.f0_bar/(con.f1_bar^2) )*(con.dT*con.f1_bar + con.mass*(ekt-1)) ;
          0 ];
      
% Ad = [Aaa, Aau;
%       Aua, Auu];
Aaa = A(1:2,1:2); 
Aau = A(1:2,3); %becomes C
Aua = A(3,1:2); 
Auu = A(3,3); %becomes A

Ba = B(1:2,1);
Bu = B(3,1);
Fa = K(1:2,1);
Fu = K(3,1);

Ea = E(1:2,1);
Eu = E(3,1);

%% Optimize to find the optimal L

% System Disturbance Constraints
d = 0.1;
m = 0.05;
perf_level = 1;

%Define Optimization Variables
L = sdpvar(size(Auu,1),size(Aau,1),'full');
err0 = sdpvar(size(Auu,1),1,'full');

delta_t0 = sdpvar(1,1,'full'); 
mu_t0 = sdpvar(size(Aaa,1),1,'full');
mu_t1 = sdpvar(size(Aaa,1),1,'full');

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
disp(sol.info)



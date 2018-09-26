function [results] = observer_comparison35(varargin)
%observer_comparison35.m
%   This set of experiments is meant to show the multiple ways that we can
%   represent the final, guaranteed set to which an estimate must lie
%   within after doing our Equalized (or Free) Recovery Design.
%   It is well known that our 

%% Constants

L = [ 1, 1, 0 ];

eta_w = 0.1;
eta_v = 0.2;

M1 = 1;

expm_under_test = [2];

%% Creating a few interesting systems to test

A = [ 1 0 0;
      1 1 0;
      -1 5 -2];
  
B = [ 1; 1; 0 ];
C = [0 1 0]; % v1(:,[1 3])';

K_sub = place(A([1 2],[1 2]),B([1 2],:),[-0.5 0.5]);

ad1 = Aff_Dyn(A,B,zeros(size(A,1),1),C,...
                eta_w,eta_v);

J = [ -1.5 0; 0 -0.5];
Q = [ [ 1 ; 2 ]*sqrt(1/5) [ -2 ; 3 ]*sqrt(1/13) ];
A = Q*J*Q';

C = [0 1];

ad2 = Aff_Dyn(A,eye(2),zeros(size(A,1),1),C,...
                eta_w, eta_v);

%% Designate System Under Test
ad = ad2;
            
%% Plot the Initial Set of States that the system can be in. 

n = size(ad.A,1);
p = size(ad.C,1);

ic_poly = Polyhedron('ub',[M1,M1],'lb',-[M1,M1]);

figure;
plot(ic_poly)
axis([-2 2 -2 2])

%% Synthesize The Optimal Guarantee that can be made at step t_0+1
[opt_out,estima] = ad.free_rec_design_pb('Min_M3',M1,10,L(1),'verbosity',0);

exp1.opt_out = opt_out;
exp1.estima = estima;
xi1_poly = Polyhedron('ub',opt_out.M3*ones(1,n),'lb',-opt_out.M3*ones(1,n));

figure;
hold on;
plot(ic_poly)
plot(xi1_poly,'color','cyan')
axis([-2 2 -2 2])

exp1.xi1_poly = xi1_poly;

%% Use Set Based Operations to Determine the True Set
F_00 = estima.F_set{1};
A_tilde = ...
    [ ad.A + F_00*ad.C, eye(n) , F_00, -eye(n) ;
      -(ad.A + F_00*ad.C), -eye(n) , -F_00, eye(n) ;
      eye(n), zeros(n,n+p+n);
      -eye(n), zeros(n,n+p+n);
      zeros(n), eye(n), zeros(n,p+n);
      zeros(n), -eye(n), zeros(n,p+n);
      zeros(p,n+n), eye(p), zeros(p,n);
      zeros(p,n+n), -eye(p), zeros(p,n)];
  
u0_0 = estima.u0_set{1};
b_tilde = ...
    [ -u0_0 ;
      u0_0 ;
      M1*ones(2*n,1);
      ad.eta_w*ones(2*n,1);
      ad.eta_v*ones(2*p,1)];
  
P1 = Polyhedron('A',A_tilde,'b',b_tilde);

F1 = P1.projection(n+n+p+[1:n]);

exp1.P1 = P1;
exp1.F1 = F1;

%% Experiment 2, Looking at the Structure of Some F matrices, when there are repeated patterns.
[opt_out,estima] = ad.free_rec_design_pb('Min_M3',M1,10,[1,1,0,1,0,1,1,0,1,1],'verbosity',0);

exp2.estima = estima;
exp2.opt_out = opt_out;

%% Assigning Results

results.exp1 = exp1;
results.exp2 = exp2;

figure;
hold on;
plot(ic_poly)
plot(xi1_poly,'color','cyan')
plot(F1,'color','magenta')
axis([-2 2 -2 2])

end
%experiment42.m
%% Description
%	The purpose of this experiment is to:
%	- test some of the new modifications to the Aff_Dyn() object.

%% Constants
h = 0.1; %Discretization step
epsil0 = 0.6; %Coefficient of restitution

%Define the Piecewise Affine System
dim = 2;

A1 = [ 1, h; 0 1 ];
B1 = [0; (h/m)];
C1 = eye(dim);
f1 = zeros(dim,1);

eta_v = 0; eta_w = 0;
Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

A2 = [1,0;0,-epsil0];

ad2 = Aff_Dyn(A2,B1,f1,C1,Pw1,Pv1);

%% Testing the mpc matrix generator
[H1,S1,~,J1,~] = get_mpc_matrices([ad1,ad2],[1,2]);
[H,S,~,f_bar] = get_mpc_matrices([ad1,ad2],[1,1,1,1,2,2,2,2]);
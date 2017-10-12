% observer_comparison4.m

clear all;
close all;
clc;

%% Constants

%Using ACC System
load('data/system_examples/acc_p.mat');

T = 2; %Time horizon of the problem we will solve in Skaf-land

pl = 0.5;
verbosity = 1;

n = size(acc.A,1);

alpha1 = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve problem with constant bound on ONLY disturbance %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('================================================')
disp('Experiment 0: Standard Bounds (ONLY Disturbance)')
disp('================================================')
disp(' ')

% Experimental Params
%--------------------
verbosity = 0;

% Optimization
%-------------

acc_e = acc;
acc_e.B = eye(size(acc.A,1));	%For error system, the input we are interested in
								% isn't modified by a B matrix
acc_e.x0 = Inf;		%The initial condition needs to be set due to my rigid,
					% old way of doing things. It will not be used when pl is given.

sol_dist_bounds = generate_skaf_controller( acc_e , T , verbosity , 'PL' , pl );

disp(['The optimal bound on ||e(' num2str(T) ')|| is ' num2str(sol_dist_bounds.opt_obj) '.' ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve problem with constant bounds on disturbance and input %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===================================================')
disp('Experiment 1: Constant Bounds on Input, Disturbance')
disp('===================================================')
disp(' ')

% Experiment Constants
%---------------------

u_bound = 0.25;
verbosity = 0;

% Set Up Optimization Variables
%------------------------------

l_diag_constr = [];
robust_constr = [];
epi_constr = [];

acc_e.x0 = sdpvar(n,1,'full');

[G,H,Cm,x0m] = create_skaf_n_boyd_matrices(acc_e,T);

if verbosity >= 1
	disp('Created Skaf and Boyd Matrices.')
end

% Perform Robust Optimization Using YALMIP "uncertain"
%-----------------------------------------------------
if verbosity >= 1
	disp('YALMIP Robust Optimization: ')
end

Q = sdpvar(size(H,2),size(Cm,1),'full');
r = sdpvar(size(H,2),1,'full');

%Disturbance vectors
w = sdpvar(n*T,1,'full');
v = sdpvar(size(acc_e.C,1)*T,1,'full');

if verbosity >= 1
	disp('- Created Optimization Variables.')
end

%Create Expressions Containing Q,r for Optimization
Pxw = (eye(n*(T+1))+H*Q*Cm)*G;
Pxv = H*Q;
x_tilde = (eye(n*(T+1)) + H*Q*Cm)*x0m + H*r;

Puw = Q*Cm*G;
Puv = Q;
u_tilde = Q*Cm*x0m + r;

u = Puw * w + Puv * v + u_tilde;

%Create Objective
R = [zeros(n,n*T) eye(n)]; %Create the standard selection matrix

objective = norm( R*(x_tilde + Pxw * w + Pxv * v) , Inf );

if verbosity >= 1
	disp('- Created Objective.')
end

%Create Constraints

%Q is lower diagonal.
for bl_row_num = 1 : T-1
	l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
											[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
end

if verbosity >= 2
	l_diag_constr
end

%Robustifying against w and v
robust_constr = robust_constr + [ -(acc_e.m) <= v <= (acc_e.m) , uncertain(v) ];
robust_constr = robust_constr + [ -acc_e.d <= w <= acc_e.d , uncertain(w) ];
robust_constr = robust_constr + [ -pl <= acc_e.x0 <= pl , uncertain(acc_e.x0) ];

if verbosity >= 2
	robust_constr
end

%Epigraph Constraints
alpha0 = sdpvar(1,1,'full');
epi_constr = [ objective <= alpha0 ];

if verbosity >=2 
	epi_constr
end

if verbosity >= 1
	disp('- Created Constraints.')

	%Solve Optimization
	disp('Solving YALMIP Robust Optimization...')
end

op_num = verbosity;
ops = sdpsettings('verbose',op_num);
sol_const_bounds = optimize(l_diag_constr+robust_constr+epi_constr + [ -u_bound <= Puw * w + Puv * v + u_tilde <= u_bound ], ...
							alpha0, ...
							ops);

%

if verbosity >= 1
	if sol_const_bounds.problem == 0
		disp(['YALMIP Robust Optimization Solved'])
	else
		%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
		disp(['YALMIP Robust Optimization NOT Solved.'])
	end
end

sol_const_bounds.Q = value(Q);
sol_const_bounds.r = value(r);
sol_const_bounds.opt_obj = value(alpha0);

sol_const_bounds.F = value( (pinv(value(eye(size(Q,1)) + Q*Cm*H)) ) * Q);
sol_const_bounds.u0 = value((eye(size(sol_const_bounds.F,1)) + sol_const_bounds.F*Cm*H) * r);

disp(['With a constraint, ||u|| <= ' num2str(u_bound) ','])
disp(['the optimal bound on ||e(' num2str(T) ')|| is ' num2str(sol_const_bounds.opt_obj) ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve with input dependent bounds on disturbance %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===================================================')
disp('Experiment 2: Input Dependent bounds on Disturbance')
disp('===================================================')
disp(' ')

% Experiment Constants
%---------------------

verbosity = 0;

clear l_diag_constr robust_constr epi_constr

% Set Up Optimization Variables
%------------------------------

l_diag_constr = [];
robust_constr = [];
epi_constr = [];

acc_e.x0 = sdpvar(n,1,'full');

[G,H,Cm,x0m] = create_skaf_n_boyd_matrices(acc_e,T);

if verbosity >= 1
	disp('Created Skaf and Boyd Matrices.')
end

% Perform Robust Optimization Using YALMIP "uncertain"
%-----------------------------------------------------
if verbosity >= 1
	disp('YALMIP Robust Optimization: ')
end

Q = sdpvar(size(H,2),size(Cm,1),'full');
r = sdpvar(size(H,2),1,'full');

%Disturbance vectors
w = sdpvar(n*T,1,'full');
v = sdpvar(size(acc_e.C,1)*T,1,'full');

if verbosity >= 1
	disp('- Created Optimization Variables.')
end

%Create Expressions Containing Q,r for Optimization
Pxw = (eye(n*(T+1))+H*Q*Cm)*G;
Pxv = H*Q;
x_tilde = (eye(n*(T+1)) + H*Q*Cm)*x0m + H*r;

Puw = Q*Cm*G;
Puv = Q;
u_tilde = Q*Cm*x0m + r;

u = Puw * w + Puv * v + u_tilde;

%Create Objective
R = [zeros(n,n*T) eye(n)]; %Create the standard selection matrix

objective = norm( R*(x_tilde + Pxw * w + Pxv * v) , Inf );

if verbosity >= 1
	disp('- Created Objective.')
end

%Create Constraints

%Q is lower diagonal.
for bl_row_num = 1 : T-1
	l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
											[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
end

if verbosity >= 2
	l_diag_constr
end

%Robustifying against w and v
robust_constr = robust_constr + [ -(acc_e.m) <= v <= (acc_e.m) , uncertain(v) ];
robust_constr = robust_constr + [ -acc_e.d <= w+alpha1*(Puw*w+Puv*v+u_tilde) <= acc_e.d , uncertain(w) ];
robust_constr = robust_constr + [ -pl <= acc_e.x0 <= pl , uncertain(acc_e.x0) ];

if verbosity >= 2
	robust_constr
end

%Epigraph Constraints
alpha0 = sdpvar(1,1,'full');
epi_constr = [ objective <= alpha0 ];

if verbosity >=2 
	epi_constr
end

if verbosity >= 1
	disp('- Created Constraints.')

	%Solve Optimization
	disp('Solving YALMIP Robust Optimization...')
end

op_num = verbosity;
ops = sdpsettings('verbose',op_num);
sol_i2d_bounds = optimize(l_diag_constr+robust_constr+epi_constr,alpha0,ops);

%

if verbosity >= 1
	if sol_i2d_bounds.problem == 0
		disp(['YALMIP Robust Optimization Solved'])
	else
		%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
		disp(['YALMIP Robust Optimization NOT Solved.'])
	end
end

sol_i2d_bounds.Q = value(Q);
sol_i2d_bounds.r = value(r);
sol_i2d_bounds.opt_obj = value(alpha0);

sol_i2d_bounds.F = value( (pinv(value(eye(size(Q,1)) + Q*Cm*H)) ) * Q);
sol_i2d_bounds.u0 = value((eye(size(sol_i2d_bounds.F,1)) + sol_i2d_bounds.F*Cm*H) * r);

disp('With constraints:')
disp('    - || w + alpha*u || <= d')
disp(['The optimal bound on ||e(' num2str(T) ')|| is ' num2str(sol_i2d_bounds.opt_obj) '.'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve with an input-dependent bound on x %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('==============================================')
disp('Experiment 3: Input-Dependent Bound on State x')
disp('==============================================')


% Experiment Constants
%---------------------

verbosity = 1;

clear l_diag_constr robust_constr epi_constr

% Set Up Optimization Variables
%------------------------------

l_diag_constr = [];
robust_constr = [];
epi_constr = [];

acc_e.x0 = sdpvar(n,1,'full');

[G,H,Cm,x0m] = create_skaf_n_boyd_matrices(acc_e,T);

if verbosity >= 1
	disp('Created Skaf and Boyd Matrices.')
end

% Perform Robust Optimization Using YALMIP "uncertain"
%-----------------------------------------------------
if verbosity >= 1
	disp('YALMIP Robust Optimization: ')
end

Q = sdpvar(size(H,2),size(Cm,1),'full');
r = sdpvar(size(H,2),1,'full');

%Disturbance vectors
w = sdpvar(n*T,1,'full');
v = sdpvar(size(acc_e.C,1)*T,1,'full');

if verbosity >= 1
	disp('- Created Optimization Variables.')
end

%Create Expressions Containing Q,r for Optimization
Pxw = (eye(n*(T+1))+H*Q*Cm)*G;
Pxv = H*Q;
x_tilde = (eye(n*(T+1)) + H*Q*Cm)*x0m + H*r;

Puw = Q*Cm*G;
Puv = Q;
u_tilde = Q*Cm*x0m + r;

x = (x_tilde + Pxw * w + Pxv * v);
u = Puw * w + Puv * v + u_tilde;

%Create Objective
R = [zeros(n,n*T) eye(n)]; %Create the standard selection matrix

objective = norm( R*(x_tilde + Pxw * w + Pxv * v) , Inf );

if verbosity >= 1
	disp('- Created Objective.')
end

%Create Constraints

%Q is lower diagonal.
for bl_row_num = 1 : T-1
	l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(acc_e.B,2)+1:bl_row_num*size(acc_e.B,2)], ...
											[bl_row_num*size(acc_e.C,1)+1:end] ) == 0 ];
end

if verbosity >= 2
	l_diag_constr
end

%Robustifying against w and v
robust_constr = robust_constr + [ -(acc_e.m) <= v <= (acc_e.m) , uncertain(v) ];
robust_constr = robust_constr + [ -acc_e.d <= w <= acc_e.d , uncertain(w) ];
robust_constr = robust_constr + [ -pl <= acc_e.x0 <= pl , uncertain(acc_e.x0) ];

if verbosity >= 2
	robust_constr
end

%Epigraph Constraints
alpha0 = sdpvar(1,1,'full');
epi_constr = [ objective <= alpha0 ];

if verbosity >=2 
	epi_constr
end

if verbosity >= 1
	disp('- Created Constraints.')

	%Solve Optimization
	disp('Solving YALMIP Robust Optimization...')
end

op_num = verbosity;
ops = sdpsettings('verbose',op_num);
sol_i2d_bounds = optimize(	l_diag_constr+robust_constr+epi_constr+[ u <= x([4:9]) ], ...
							alpha0, ...
							ops);

%

if verbosity >= 1
	if sol_i2d_bounds.problem == 0
		disp(['YALMIP Robust Optimization Solved'])
	else
		%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
		disp(['YALMIP Robust Optimization NOT Solved.'])
	end
end

sol_i2d_bounds.Q = value(Q);
sol_i2d_bounds.r = value(r);
sol_i2d_bounds.opt_obj = value(alpha0);

sol_i2d_bounds.F = value( (pinv(value(eye(size(Q,1)) + Q*Cm*H)) ) * Q);
sol_i2d_bounds.u0 = value((eye(size(sol_i2d_bounds.F,1)) + sol_i2d_bounds.F*Cm*H) * r);

disp('With constraints:')
disp('    - || w + alpha*u || <= d')
disp(['The optimal bound on ||e(' num2str(T) ')|| is ' num2str(sol_i2d_bounds.opt_obj) '.'])










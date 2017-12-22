function [ controller , optim_data ] = achieve_eq_recovery_for( varargin )
% achieve_eq_recovery_for.m
%   The purpose of this function is to take a discrete linear system with process and measurement noise
%   and synthesize a controller (F,u_0) such that given ||xi(0)|| <= M1, T, we can guarantee 
%	||\xi(t)|| <= M2 for t in [1,T-1] and ||\xi(T)||\leq M1.
%	Important Distinction: We will take the consntants M1,M2, and T as being given
%		
%	Model Assumed in this Code:
%
%     x[k+1] = Ax[k] + Bu[k] + w[k]
%     y[k]   = Cx[k] + v[k]
%
%	Objective Assumed in this Code:
%
%
%							
%	Potential Usage:
%		- achieve_eq_recovery_for( sys , t_horizon , M1 , M2 , verbosity )
%
%   Inputs:
%		sys:		- This is a struct containing information about the time horizon of the model, as well as
%				  	  some information on bounds of the process and measurement noises.
%					- Members: .A,.B,.C,.m,.d,
%
%		t_horizon:	- Time horizon that we are using for our objective/design.
%
%		verbosity:	- This is a number (currently from 0 to 2), that determines how many messages are
%					 	  sent to the user during operation. (The higher the number the more messages sent.)
%
%		(Qualifiers)
%
%	Outputs:
%		controller 	- Struct containing the final F and u0 matrices.
%		optim_data 	- Struct containing the YALMIP results, and important matrix values. (In one variable so that
%				  	  this can be easily saved.)

%%%%%%%%%%%%%%%%%%%
%% Process Input %%
%%%%%%%%%%%%%%%%%%%

if nargin ~= 5
	%The first expression (containing mod() ) reflects that we expect to have 3 + 2*n number of arguments (where n=0,1,2,...)

	%If the function is improperly called, tell the user.
	error(['Improper number of arguments given or improper qualifiers given. (Received ' num2str(nargin) ')']);
end

sys 		= varargin{1};
t_horizon 	= varargin{2};
M1 			= varargin{3};
M2 			= varargin{4};
verbosity 	= varargin{5};

%%%%%%%%%%%%%%%
%% Constants %%
%%%%%%%%%%%%%%%

n = size(sys.A,1);
p = size(sys.C,1);
d_u = size(sys.B,2);

%Select matrix
select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synthesize Controller %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables
delta 	= sdpvar(n*t_horizon,1,'full');
mu 		= sdpvar(p*t_horizon,1,'full');
sys.x0 	= sdpvar(n,1,'full');

alpha_l 	= sdpvar(t_horizon+1,1,'full');

% Create Trajectory Matrices
[S,~,Cm,xi0m] = create_skaf_n_boyd_matrices(sys,t_horizon);

% Create Objective (Epigraph constraint)
Q = sdpvar(size(S,2),size(Cm,1),'full');
r = sdpvar(size(S,2),1,'full');

Pxd = (eye(n*(t_horizon+1))+S*Q*Cm)*S ;
Pxm = S*Q;
xi_tilde = (eye(n*(t_horizon+1)) + S*Q*Cm)*xi0m + S*r;

objective = norm( select_m(t_horizon,t_horizon)*(xi_tilde + Pxd * delta + Pxm * mu) , Inf );
interm_norms = norm([ eye(n*t_horizon) zeros(n*t_horizon,n) ]*( xi_tilde + Pxd * delta + Pxm * mu ) , Inf);
% interm_vals = [ eye(n*T_available) zeros(n*T_available,n) ]*( xi_tilde + Pxd * delta(n*T_missing+1:end,1) + Pxm * mu(p*T_missing+1:end,1) );

epi_constr2 = [ objective <= M1 , interm_norms <= M2 ];

% Feasibility Constraints
% +++++++++++++++++++++++

feasib_constrs = [];
feasib_constrs = feasib_constrs + [-M2 <= alpha_l <= M2];

% Create Robustification Constraints
% ++++++++++++++++++++++++++++++++++

robust_constrs = [];
robust_constrs = robust_constrs + [ -sys.d <= delta <= sys.d , uncertain(delta) ];
robust_constrs = robust_constrs + [ -sys.m <= mu <= sys.m , uncertain(mu)];
robust_constrs = robust_constrs + [ -M1 <= sys.x0 <= M1 , uncertain(sys.x0) ];

% Create Causality (Lower Diagonal) Constraint
% ++++++++++++++++++++++++++++++++++++++++++++
l_diag_constr = [];
for bl_row_num = 1 : t_horizon-1
	l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*d_u+1:bl_row_num*d_u], ...
											[bl_row_num*p+1:end] ) == 0 ];
end

% Create Graph Constraints
% ++++++++++++++++++++++++
graph_constrs = [];
for i = 0:t_horizon
	graph_constrs = graph_constrs + [ norm(select_m(i,t_horizon)*( xi_tilde + Pxd * delta + Pxm * mu ) , Inf) <= alpha_l(i+1) ];
end

% Optimize!
ops = sdpsettings('verbose',verbosity);
optim1 = optimize(epi_constr2+robust_constrs+l_diag_constr+graph_constrs ,sum(alpha_l),ops);

if verbosity > 0
	if optim1.problem ~= 0
		error('Optimization was NOT successfully solved.')
	end
end
%%%%%%%%%%%%%%%%%%%%
%% Saving Results %%
%%%%%%%%%%%%%%%%%%%%

optim_data.sol = optim1;
optim_data.max_err_norms = value(alpha_l);
optim_data.Q = value(Q);
optim_data.r = value(r);

controller.F = value( (inv(value(eye(size(Q,1)) + Q*Cm*S)) ) * Q);
controller.u0 = value( inv( value(eye(size(Q,1)) + Q*Cm*S)) * r );


end
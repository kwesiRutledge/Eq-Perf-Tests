function [ results ] =  generate_yong_controller( varargin )
%	generate_yong_controller.m
%		Description:
%			This function should create the matrices necessary for a dynamic
%			output feedback controller, as defined by ASU Prof. Sze Zheng Yong.
%
%		Potential Usage:
%			generate_yong_controller( sys , t_horizon , verbosity )
%			generate_yong_controller( ... , 's' , dim_s )
% 			generate_yong_controller( ... , 'PL' , perf_level)
%
%		Inputs:
%			sys			- Struct containing system matrices and other valuable system
%				  		  properties
%						- Must contain: .A,.B,.C,.m,.d,.x0
%
%			t_horizon	- Time horizon that we are using for our objective/design.			
%
%			verbosity	- This is a number (currently from 0 to 2), that determines how many messages are
%					 	  sent to the user during operation. (The higher the number the more messages sent.)
%
%		Outputs:
%			results	- Struct containing results as well as other intermediate values involved in the design.
%					  (e.g. Optimization solved/not solved,value of optimal optimization variables, etc.)

% Input Processing
%-----------------

sys = varargin{1};
t_horizon = varargin{2};
verbosity = varargin{3};

feasible_strs = {'s','PL'};

%The first expression (containing mod() ) reflects that we expect to have 3 + 2*n number of arguments (where n=0,1,2,...)
if mod(nargin-3,2) 
	%If the function is improperly called, tell the user.
	error('Improper number of arguments given.')
end


%The second condition expresses that if the expression has more than 3 arguments, we expect one of those arguments to be
%	one of our qualifiers (e.g. 'PL' or 'R')
for ind = (3+1):2:nargin
	if (~any(strcmp(feasible_strs,varargin{ind})) )
		error('Unrecognized string in input.')
	end
end

%Check the fields of sys
sys_fields = { 'A','B','C','m','d','x0'};
for ind = 1 : length(sys_fields)
	if (~isfield(sys,sys_fields{ind}))
		error(['Missing the ' sys_fields{ind} ' field of input system.'])
	end
end

% Constants
%----------

n = size(sys.A,1);

%Set s to default or to user's request
if ~any(strcmp(varargin,'s'))
	dim_s = n;
else
	dim_s = varargin{ find(strcmp(varargin,'s')) + 1 };
end

if any(strcmp(varargin,'PL'))
	%Remove x0 field from system.
	rmfield(sys,'x0');
	sys.x0 = sdpvar(n,1,'full');
	perf_level = varargin{ find(strcmp(varargin,'PL')) + 1 };
end

% Creating Skaf Matrices with Yong's Modifications
%-------------------------------------------------

dyn_obs_sys = dyn_obs_ify(sys,dim_s);

[G,H,Cm,x0m] = create_skaf_n_boyd_matrices(dyn_obs_sys,t_horizon);

if verbosity >= 1
	disp('Created Skaf Constants')
end

% Create YALMIP Optimization Variables
%-------------------------------------

if verbosity >= 1
	disp('YALMIP Robust Optimization: ')
end

Q = sdpvar(size(H,2),size(Cm,1),'full');
r = sdpvar(size(H,2),1,'full');

%Disturbance vectors								LOOK INTO THE SIZE/STRUCTURE OF THIS LATER IF THERE ARE STILL ERRRORS
% w0 = sdpvar(n*t_horizon,1,'full');
% v0 = sdpvar(size(sys.C,1)*t_horizon,1,'full');

% w = []; v = [];
% for t = 1 : t_horizon
% 	w = [ w ; w0((t-1)*n+1:t*n) ; zeros(dim_s,1)];
% 	v = [ v ; zeros(dim_s,1) ; v0((t-1)*size(sys.C,1) + 1:t*size(sys.C,1))];
% end

w_y = sdpvar(n,t_horizon,'full');
v_y = sdpvar(size(sys.C,1),t_horizon,'full');

%Make w and v into vectors
w = []; v = [];
for t = 1 : t_horizon
	w = [ w ; w_y(:,t) ; zeros(dim_s,1) ];
	v = [ v ; zeros(dim_s,1) ; v_y(:,t) ];
end

%whos('Q','r','w','v');

if verbosity >= 1
	disp('- Variables Created')
end

% Create Compound Expressions with YALMIP Variables
%--------------------------------------------------
Pxw = (eye((n+dim_s)*(t_horizon+1))+H*Q*Cm)*G;
Pxv = H*Q;
x_tilde = (eye((n+dim_s)*(t_horizon+1)) + H*Q*Cm)*x0m + H*r;

R = [zeros(n,(n+dim_s)*t_horizon) eye(n) zeros(n,dim_s)];

objective = norm( R*(x_tilde + Pxw * w + Pxv * v) , Inf );

if verbosity >= 1
	disp('- Objective Created')
end

% Create YALMIP Optimization's Constraints
%-----------------------------------------

l_diag_constr = []; robust_constr = []; epi_constr = []; disturb_constrs = [];

%Lower Diagonal Constraint

for bl_row_num = 1 : t_horizon-1
	l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(dyn_obs_sys.B,2)+1:bl_row_num*size(dyn_obs_sys.B,2)], ...
											[bl_row_num*size(dyn_obs_sys.C,1)+1:end] ) == 0 ];
end

if verbosity >= 2
	l_diag_constr
end

%Robustifying against w0 and v0
for t = 1 : t_horizon
	robust_constr = robust_constr + [ -sys.m <= v_y(:,t) <= sys.m , uncertain(v_y(:,t)) ];
	robust_constr = robust_constr + [ -sys.d <= w_y(:,t) <= sys.d , uncertain(w_y(:,t)) ];
end

if any(strcmp(varargin,'PL'))
	robust_constr = robust_constr + [ -perf_level <= dyn_obs_sys.x0 <= perf_level , uncertain(dyn_obs_sys.x0) ];
end

if verbosity >= 2
	robust_constr
end

%Epigraph Constraint
alpha0 = sdpvar(1,1,'full');
epi_constr = [ objective <= alpha0 ];

if verbosity >=2 
	epi_constr
end

%Disturbance Constraints (Used because the disturbances has some forced zeros)
% for t = 1 : t_horizon
% 	disturb_constrs = disturb_constrs + [ w((t-1)*(n+dim_s)+n+1:t*(n+dim_s)) == 0 ];
% 	disturb_constrs = disturb_constrs + [ v((t-1)*(size(sys.C,1)+dim_s)+1:(t-1)*((size(sys.C,1)+dim_s))+dim_s ) == 0 ];
% end

% if verbosity >= 2
% 	disturb_constrs
% end

if verbosity >= 1
	disp('- Constraints Created')
end

% YALMIP Optimization
%--------------------

ops = sdpsettings('verbose',verbosity);
results.sol_robust = optimize(l_diag_constr+robust_constr+epi_constr,alpha0,ops);

if verbosity >= 1
	if results.sol_robust.problem == 0
		disp(['YALMIP Robust Optimization Solved'])
	else
		error(['YALMIP Robust Optimization NOT Solved.'])
	end
end

% Saving Results
%---------------
results.Q = value(Q);
results.Q( isnan(value(Q)) ) = 0;
results.r = value(r);
results.r( isnan(value(r)) ) = 0;
results.opt_obj = value(alpha0);
value(alpha0);

if any(isnan(value(Q))) | any(isnan(value(r))) | any(isnan(value(alpha0)))
	%Warn the user
	warning('There are NaNs in Q,r,or opt_obj. This can be because of bad inputs or because some part of the optimization variables are unused. Setting NaNs to zero.')
end

results.F = value( (pinv(value(eye(size(results.Q,1)) + results.Q*Cm*H)) ) * results.Q);
results.u0 = value((eye(size(results.F,1)) + results.F*Cm*H) * results.r);

results.G = G;
results.H = H;
results.Cm = Cm;
results.x0m = x0m;


end
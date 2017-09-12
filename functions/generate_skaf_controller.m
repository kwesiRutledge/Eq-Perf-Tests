function [ results ] = generate_skaf_controller( varargin )
% generate_skaf_controller.m
%   The purpose of this function is to take a discrete linear system with process and measurement noise
%   and minimize the objective ||x[T]||. Its output includes the result of YALMIP Optimization as well as
%   the factors (Q,r) of the feedback matrix.
%
%	Model Assumed in this Code:
%
%     x[k+1] = Ax[k] + Bu[k] + w[k]
%     y[k]   = Cx[k] + v[k]
%
%	Objective Assumed in this Code:
%
%		||R*x||_{\infty},	where R can be arbitrarily chosen
%							
%	Potential Usage:
%		- generate_skaf_controller( sys , t_horizon , verbosity )
%		- generate_skaf_controller( sys , t_horizon , verbosity , 'PL' , perf_level )
%		- generate_skaf_controller( sys , t_horizon , verbosity , 'R' , R )
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
%		perf_level:	- Must make the previous input before this the string 'PL'.
%					- Positive real number
%
%		R:			- The selection matrix for the trajectory x
%
%	Outputs:
%		results - Struct containing the YALMIP results, and important matrix values. (In one variable so that
%				  this can be easily saved.)

% Process Inputs
%---------------

if any( strcmp(varargin,'PL') )
	%Option for Stochastic x0 is given

	%Find Performance Level
	PL_loc = find( strcmp(varargin,'PL') );

	%Save Performance Level
	perf_level = varargin{PL_loc+1};
	stochastic_x0 = true;
else
	stochastic_x0 = false;
end

sys 		= varargin{1};
t_horizon 	= varargin{2};
verbosity 	= varargin{3};

% Constants
%----------

n = size(sys.A,1);

% Constraints 
%------------

l_diag_constr = [];
robust_constr = [];
epi_constr = [];

% Modify x0 if necessary:
if stochastic_x0
	%Make x0 an optimization Variable
	sys.x0 = sdpvar(size(sys.A,1),1,'full');

	%Create constraints on x0 (based on performance level)
	robust_constr = robust_constr + [ -perf_level <= sys.x0 <= perf_level , uncertain(sys.x0) ];
end

% Create Constant Matrices Based On Model
%----------------------------------------

[G,H,Cm,x0m] = create_skaf_n_boyd_matrices(sys,t_horizon);

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
w = sdpvar(n*t_horizon,1,'full');
v = sdpvar(size(sys.C,1)*t_horizon,1,'full');

if verbosity >= 1
	disp('- Created Optimization Variables.')
end

%Create Expressions Containing Q,r for Optimization
Pxw = (eye(n*(t_horizon+1))+H*Q*Cm)*G;
Pxv = H*Q;
x_tilde = (eye(n*(t_horizon+1)) + H*Q*Cm)*x0m + H*r;

%Create Objective
R = [zeros(n,n*t_horizon) eye(n)]; %Create selection matrix
objective = norm( R*(x_tilde + Pxw * w + Pxv * v) , Inf );

if verbosity >= 1
	disp('- Created Objective.')
end

%Create Constraints

%Q is lower diagonal.
for bl_row_num = 1 : t_horizon-1
	l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(sys.B,2)+1:bl_row_num*size(sys.B,2)], ...
											[bl_row_num*size(sys.C,1)+1:end] ) == 0 ];
end

if verbosity >= 2
	l_diag_constr
end

%Robustifying against w and v
robust_constr = robust_constr + [ -sys.m <= v <= sys.m , uncertain(v) ];
robust_constr = robust_constr + [ -sys.d <= w <= sys.d , uncertain(w) ];

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
results.sol_robust = optimize(l_diag_constr+robust_constr+epi_constr,alpha0,ops);

if verbosity >= 1
	if results.sol_robust.problem == 0
		disp(['YALMIP Robust Optimization Solved'])
	else
		%error(['YALMIP Robust Optimization #' num2str(T) ' NOT Solved.'])
		disp(['YALMIP Robust Optimization NOT Solved.'])
	end
end

%Save results
results.Q = value(Q);
results.r = value(r);
results.opt_obj = value(alpha0);

results.F = value( (pinv(value(eye(size(Q,1)) + Q*Cm*H)) ) * Q);
results.u0 = value((eye(size(results.F,1)) + results.F*Cm*H) * r);
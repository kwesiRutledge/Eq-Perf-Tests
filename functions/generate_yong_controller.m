function [ results ] =  generate_yong_controller( varargin )
%	generate_yong_controller.m
%		Description:
%			This function should create the matrices necessary for a dynamic
%			output feedback controller, as defined by ASU Prof. Sze Zheng Yong.
%
%		Potential Usage:
%			generate_yong_controller( sys , t_horizon , verbosity )
%			generate_yong_controller( sys , t_horizon , verbosity , 's' , dim_s )
%			generate_yong_controller( sys , t_horizon , verbosity , 'solver' , solver_name )
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

feasible_strs = {'s','solver'};

%The first expression (containing mod() ) reflects that we expect to have 3 + 2*n number of arguments (where n=0,1,2,...)
if mod(nargin-3,2) 
	%If the function is improperly called, tell the user.
	error('Improper number of arguments given.')
end


%The second condition expresses that if the expression has more than 3 arguments, we expect one of those arguments to be
%	one of our qualifiers (e.g. 'PL' or 'R')
for ind = (3+1):2:nargin
	if (~strcmp(varargin{ind},'PL')) & ( ~strcmp(varargin{ind},'R') )
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

%% Constants

if ~any(strcmp(varargin,'solver'))
	%Default Solver will be Yalmip
	solver_name = 'yalmip';
else
	solver_name = varargin{ find( strcmp(varargin,'solver')) + 1 };
end

% CVX Solution
%-------------

if strcmp(solver_name,'yalmip')
	% Create Modified System
	%~~~~~~~~~~~~~~~~~~~~~
	dyn_obs_sys = dyn_obs_ify(sys)
else
	

end
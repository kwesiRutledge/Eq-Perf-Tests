function [ sys_w_dyn_obs ] = dyn_obs_ify( varargin )
	% convert_to_dyn_obs_sys.m
	%	This function will use the A,B,C matrices of input_sys to
	% 	create the modified system (the output) with a dynamic output
	%	based on input_sys.
	%
	%	Usage:
	%		- convert_to_dyn_obs_sys(input_sys)
	%		- convert_to_dyn_obs_sys(input_sys,s)
	%		- convert_to_dyn_obs_sys(input_sys,s,phi0)
	%
	%	Requirements:
	%		- input_sys is a struct containing members: .A,.B,.C,.x0
	%
	%	Input:
	%		input_sys 	- NOT OPTIONAL
	%				  	- Struct containing system matrices as members.
	%				      (Must have .A,.B,.C defined along with initial
	%					  condition .x0)
	%
	%		s 			- Optional Input
	%					- Size of the "output state", related to memory of error
	%
	%		phi0 		- Optional Input
	%					- The initial condition of the output state.
	%					- Size should match s
	% 		
	%	Output:
	%		- sys_w_dyn_obs (the only output) will be a struct with
	%		  the following members defined: .A,.B,.C,.x0

	% Input Handling
	%---------------

	input_sys = varargin{1};

	% Checking the fields of input_sys

	if (~isfield(input_sys,'A')) | (~isfield(input_sys,'B')) | (~isfield(input_sys,'C')) | (~isfield(input_sys,'x0'))
		error('input_sys was not properly defined.')
	end

	if nargin > 3
		error('This function accepts a maximum of 2 arguments.')
	end

	if nargin == 3
		if size(varargin{3},1) ~= varargin{2}
			error('The size of initial condition phi0 (input 3) does not match s (input 2).')
		end
	end

	% Constants
	%----------

	n = size(input_sys.A,1);

	if nargin == 2
		%the size of the hidden state is given
		s = varargin{2};
	elseif nargin == 1
		s = n;
	else
		error('Unexpected number of arguments.')
	end
		
	% Creating Matrices
	%------------------
	sys_w_dyn_obs.A = [ input_sys.A zeros(n,s) ; zeros(s,n) zeros(s) ];
	sys_w_dyn_obs.B = [ zeros(n,s) input_sys.B ; eye(s) zeros(s,size(input_sys.B,2)) ];
	sys_w_dyn_obs.C = [ zeros(s,n) eye(s) ; input_sys.C zeros(size(input_sys.C,1),s) ];

	% Creating Initial Condition
	%---------------------------

	if nargin <= 2
		sys_w_dyn_obs.x0 = [ input_sys.x0 ; zeros(s,1) ];
	else
		sys_w_dyn_obs.x0 = [ input_sys.x0 ; varargin{3} ];
	end

end


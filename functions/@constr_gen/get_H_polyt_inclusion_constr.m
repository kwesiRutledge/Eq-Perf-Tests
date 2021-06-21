function [varargout] = get_H_polyt_inclusion_constr( varargin )
	%Description:
	%	Provides the constraints and optimization variables needed to define the constraint: \mathbb{X} \subseteq \mathbb{Y}
	%	where \mathbb{X} = \{ x | H_x * x <= h_x } and \mathbb{Y} = \{ y | H_y * y <= h_y }
	%
	%Usage:
	%	[ Lambda , constrs ] = cg.get_H_polyt_inclusion_constr( H_x , h_x , H_y, h_y )
	%	[ Lambda , eq_constraint , le_constraint ] = cg.get_H_polyt_inclusion_constr( H_x , h_x , H_y, h_y , 'ReturnValueType' , 'Symbolic' )
	%
	%Inputs:
	%	H_x 	- Hyperplane representation of the x-polytope's normal vectors
	%	h_x 	- Hyperplane representaiton of the X Polytope's offset values
	%	H_y 	- Etc.
	%	h_y 	- Etc.
	%	eq_constraint 	- an array of symbolic variables representing the equality constraint Lambda * H_x - H_y = 0 (represents only the left hand side)
	%	ineq_constraint - an array of symbolic variables representing the equality constraint Lambda * h_x - h_y <= 0 (represents only the left hand side)

	%% Accept Inputs

	[ obj , H_x , h_x , H_y, h_y , options ] = get_H_polyt_inclusion_constr_input_processing( varargin{:} );

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size( H_x , 2 );

	q_x = size(H_x,1);
	q_y = size(H_y,1);

	if n ~= size(H_y,2)
		error(['The two polytopes appear to be in different dimension. X in dimension ' num2str(n) ', while Y in dimension ' num2str(size(H_y,2))])
	end

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	switch options.ReturnValueType
		case 'YALMIP'
			Lambda = sdpvar(q_y,q_x,'full');

			constrs = [];
			constrs = constrs + [ Lambda >= 0 ]; %Lambda must be nonnegative
			constrs = constrs + [ Lambda * H_x == H_y ];
			constrs = constrs + [ Lambda * h_x <= h_y ];

			varargout{1} = Lambda;
			varargout{2} = constrs;

		case 'Symbolic'
			Lambda = sym('Lambda',[q_y,q_x]);

			eq_constraint = Lambda * H_x - H_y; % == 0
			le_constraint = Lambda * h_x - h_y; % <= 0

			varargout = { Lambda, eq_constraint , le_constraint };

		otherwise
			error(['Unexpected ReturnValueType ' options.ReturnValueType '.' ])
	end

	

end

function [ obj , H_x , h_x , H_y, h_y , options ] = get_H_polyt_inclusion_constr_input_processing( varargin )

	% Input Processing

	obj = varargin{1};
	H_x = varargin{2};
	h_x = varargin{3};
	H_y = varargin{4};
	h_y = varargin{5};

	% Defaults

	options = [];
	options.ReturnValueType = 'YALMIP';

	% Input Processing

	if nargin > 5
		argidx = 6;
		while argidx <= nargin
			switch varargin{argidx}
				case 'ReturnValueType'
					options.ReturnValueType = varargin{argidx+1};
					argidx = argidx + 2;
				otherwise
					error(['Unrecognized input to the function: ' varargin{argidx} ])
			end
		end
	end

end
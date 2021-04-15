function [tf] = p_subseteq_pu( varargin )
	% Description:
	%	Determines if the polyhedron p_in is a subset of the Polyhedron pu_in.
	%	Uses optimization/MILP.
	% Usage:
	%	[tf] = p_subseteq_pu( p_in , pu_in )
	%	[tf] = p_subseteq_pu( p_in , pu_in , 'verbosity' , 0 )

	%% Input Processing

	[ p_in , pu_in , verbosity ] = input_processing_p_subseteq_pu( varargin{:} );

	%% Constants

	dim = p_in.Dim;
	union_count = pu_in.Num;

	%% Algorithm

	x = sdpvar(dim,1,'full');
	b = binvar(union_count,1);

	% Constrain x to be in p_in
	x_in_p_constraint = [ p_in.A * x <= p_in.b ];
	if ~isempty(p_in.Ae)
		x_in_p_constraint = x_in_p_constraint + [ p_in.Ae * x == p_in.be ];
	end

	% Constrain x so that if it is in any individual element of the polyunion
	% pu_in, then the binary flag is triggered.

	x_in_pu_constraint = [];
	for set_index = 1:union_count

		temp_set = pu_in.Set(set_index);

		x_in_pu_constraint = x_in_pu_constraint + [ iff( temp_set.A * x <= temp_set.b , b(set_index) == 1 ) ];
		if ~isempty(temp_set.Ae)
			x_in_pu_constraint = x_in_pu_constraint + [ iff( temp_set.Ae * x == temp_set.be , b(set_index) == 1 ) ];
		end
	end

	% Create objective
	objective = sum(b);

	% Solve Optimization Flag
	ops = sdpsettings('verbose',1,'debug',1);
	ops = sdpsettings(ops,'solver','gurobi');

	optim0 = optimize(x_in_p_constraint+x_in_pu_constraint,objective,ops);

	if value(objective) == 0
		tf = false;
	else
		tf = true;
	end

end

function [ p_in , pu_in , verbosity ] = input_processing_p_subseteq_pu( varargin )
	%Description:
	%	Process Inputs

	%% Input Checking

	p_in = varargin{1};
	pu_in = varargin{2};

	if ~isa(p_in,'Polyhedron')
		error(['Expected p_in to be a Polyhedron object. Received ' class(p_in) '.' ])
	end

	if ~isa(pu_in,'PolyUnion')
		error(['Expected pu_in to be a PolyUnion object. Received ' class(pu_in) '.' ])
	end

	if p_in.Dim ~= pu_in.Dim
		error(['The dimension of p_in (' num2str(p_in.Dim) ') is different from the dimension of pu_in (' num2str(pu_in.Dim) ').' ])
	end

	%% Setting Defaults

	verbosity = 0;

	%% Adjusting For Inputs
	varargin_idx = 3;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case 'verbosity'
				verbosity = varargin{varargin_idx+1};
				varargin_idx = varargin_idx + 2;
			otherwise
				error(['Unexpected input to p_subseteq_pu: ' varargin{varargin_idx} ])
		end
	end

end
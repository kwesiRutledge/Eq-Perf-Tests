function [ z , mccormick_constraints , mccormick_binvars ] = create_mccormick_envelope( varargin )
	%create_mccormick_envelope
	%Description:
	%	Creates a McCormick Envelope based upon the inputs x and y which can be either:
	%	- a scalar and a scalar, respectively
	%	- a vector and a vector, respectively, or
	%	- a matrix and a vector, respectively.
	%	The McCormick envelope creates a new variable z which is approximately equal to the product x * y.
	%
	%Usage:
	%	[ z , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope(x,y,grid_def_in)
	%	[ z , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope(x,y,grid_def_in, 'eta_z_bounds' , eta_z_bound_array )
	%	[ z , mccormick_constraints , mccormick_binvars ] = cg.create_mccormick_envelope(x,y,grid_def_in, 'debug_flag' , 1 )
	%
	%Inputs:
	%	- grid_definition: A struct with the following members: .x_lb, .x_ub, .y_lb, .y_ub, .NumberOfRegionsInX, .NumberOfRegionsInY
	%					   The size of each member should match the corresponding input variable (e.g. if x is a m x n matrix, then 
	%					   .x_lb, .x_ub, and .NumberOfRegionsInX should all be m x n matrices.)

	%%%%%%%%%%%%%%%%%%
	%% Check Inputs %%
	%%%%%%%%%%%%%%%%%%

	if nargin < 4
		error('Require at least 4 input arguments.')
	end

	cg = varargin{1};
	x = varargin{2};
	y = varargin{3};
	grid_definition = varargin{4};

	check_mccormick_envelope_inputs( x , y , grid_definition );

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	sizeX = size(x);
	sizeY = size(y);

	x_is_scalar = all(sizeX == [1,1]);
	y_is_scalar = all(sizeY == [1,1]);

	x_is_vector = any(sizeX == 1) && ~x_is_scalar;
	y_is_vector = any(sizeY == 1) && ~y_is_scalar;

	x_is_matrix = all(sizeX ~= 1);

	dim = length(y);

	mccormick_envelope_settings = handle_create_mccormick_envelope_arguments(varargin{:});

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	mccormick_binvars = {};

	if x_is_scalar && y_is_scalar
		[ z , mccormick_constraints , mccormick_binvars ] = create_mccormick_envelope_for_scalar_product( x , y , grid_definition );
	elseif x_is_vector && y_is_vector
		[ z , mccormick_constraints , mccormick_binvars ] = create_mccormick_envelope_for_vector_product( x , y , grid_definition , mccormick_envelope_settings.eta_z_bounds );
	elseif x_is_matrix && y_is_vector

		temp_mccormick_constraints = {};
		temp_mccormick_binvars = {};
		
		mccormick_constraints = [];
		z = sdpvar(sizeX(1),1,'full');

		for row_index = 1:sizeX(1)
			% Redefine Grid w.r.t. rows of the matrix.
			temp_grid_def_i = struct( ...
				'x_lb',grid_definition.x_lb(row_index,:)', ...
				'x_ub',grid_definition.x_lb(row_index,:)', ...
				'y_lb',grid_definition.y_lb, ...
				'y_ub',grid_definition.y_ub, ...
				'NumberOfRegionsInX',grid_definition.NumberOfRegionsInX(row_index,:)', ...
				'NumberOfRegionsInY',grid_definition.NumberOfRegionsInY );

			[ temp_zi(row_index) , temp_mccormick_constraints{row_index} , temp_mccormick_binvars{row_index} ] = create_mccormick_envelope_for_vector_product( ...
				x(row_index,:) , y , ...
				temp_grid_def_i , mccormick_envelope_settings.eta_z_bounds );
			
			% Create output variables.
			mccormick_constraints = mccormick_constraints + temp_mccormick_constraints{row_index};
			mccormick_constraints = mccormick_constraints + [ z(row_index) == temp_zi(row_index) ];

			mccormick_binvars = {mccormick_binvars{:},temp_mccormick_binvars{row_index}{:}};

			if mccormick_envelope_settings.debug_flag > 1
				disp(['[create_mccormick_envelope.m]: Completed ' num2str(row_index) '/' num2str(sizeX(1)) ' mccormick envelope components in matrix equality.' ])
			end

		end

		mccormick_binvars = temp_mccormick_constraints;
	else
		error('We should never reach this part. There is a bug in the code!')
	end


end

function check_mccormick_envelope_inputs( x , y , grid_definition )
	%Description:
	%	This function checks to make sure that the McCormick Envelope function's inputs are consistent.
	%	If they are not, it produces an error.

	%% Constants

	sizeX = size(x);
	sizeY = size(y);

	grid_def_fieldnames = {'x_lb','x_ub','y_lb','y_ub','NumberOfRegionsInX','NumberOfRegionsInY'};

	%% Algorithm

	x_is_scalar = all(sizeX == [1,1]);
	y_is_scalar = all(sizeY == [1,1]);

	x_is_vector = any(sizeX == 1) && ~x_is_scalar;
	y_is_vector = any(sizeY == 1) && ~y_is_scalar;

	x_is_matrix = all(sizeX ~= 1);

	% Check that the inputs x and y are compatible generally.
	if ~(x_is_scalar && y_is_scalar) && ~(x_is_vector && y_is_vector) && ~(x_is_matrix && y_is_vector)
		error('It must be the case that x and y are of compatible dimensions. Make sure that x,y are both scalars, both vectors, or a matrix and a vector.')
	end

	% Check that the grid_definition object has all of its necessary fields.
	for fieldname_index = 1:length(grid_def_fieldnames)

		if ~isfield( grid_definition , grid_def_fieldnames{fieldname_index} )
			error([ 'The field name ' grid_def_fieldnames{fieldname_index} ' is not part of the input grid_definition.' ])
		end

	end

	% Scalar-Scalar checks
	if x_is_scalar && y_is_scalar



	end

	% Vector-Vector checks
	if x_is_vector && y_is_vector

		if length(x) ~= length(y)
			error(['Expected the lengths of x and y to be the same. Instead, length(x) = ' num2str(length(x)) ' and length(y) = ' num2str(length(y)) '.'  ])
		end

		% Checking grid_definition
		size_x_lb = size(grid_definition.x_lb);
		size_x_ub = size(grid_definition.x_ub);
		size_y_lb = size(grid_definition.y_lb);
		size_y_ub = size(grid_definition.y_ub);

		if ~all(size_x_lb == size(x))
			error(['The x_lb field (size ' num2str(size_x_lb) ') of grid_definition is not of the same size as x (' num2str(size(x)) ').'])
		end

		if ~all(size_x_ub == size(x))
			error(['The x_ub field (size ' num2str(size_x_ub) ') of grid_definition is not of the same size as x (' num2str(size(x)) ').'])
		end

		if ~all(size_y_lb == size(y))
			error(['The y_lb field (size ' num2str(size_y_lb) ') of grid_definition is not of the same size as y (' num2str(size(y)) ').'])
		end

		if ~all(size_y_ub == size(y))
			error(['The y_ub field (size ' num2str(size_y_ub) ') of grid_definition is not of the same size as y (' num2str(size(y)) ').'])
		end

	end

	% Check that all objects are sdpvars

	if ~strcmp(class(x),'sdpvar')
		error(['x is expected to be an sdpvar, but instead its class is ' class(x) '.' ])
	end

	if ~strcmp(class(y),'sdpvar')
		error(['y is expected to be an sdpvar, but instead its class is ' class(y) '.' ])
	end


end

function [ z , mccormick_constraints , mccormick_binvars ] = create_mccormick_envelope_for_scalar_product( x , y , grid_definition )
	%create_mccormick_envelope_for_scalar_product
	%Description:
	%	Creates a mccormick envelope to represent the product of scalar YALMIP variables x * y.
	%
	%Inputs:
	%	grid_definition - A struct with members .x_lb, .x_ub, .y_lb, .y_ub, .NumberOfRegionsInX, .NumberOfRegionsInY
	%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	num_mccormick_envelopes = grid_definition.NumberOfRegionsInX * grid_definition.NumberOfRegionsInY;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%% Create Optimization Variables
	mccormick_binvars = binvar(num_mccormick_envelopes,1,'full');
	z = sdpvar(1,1,'full');

	%Create Constraints
	mccormick_constraints = [];
	for envelope_idx = 1:num_mccormick_envelopes
		[ x_bounds , y_bounds ] = get_envelope_bounds_helper_lbs( ...
			envelope_idx , ...
			[ grid_definition.NumberOfRegionsInX, grid_definition.NumberOfRegionsInY] , ...
			grid_definition.x_lb , grid_definition.x_ub , ...
			grid_definition.y_lb , grid_definition.y_ub );

		x_L = x_bounds(1); x_U = x_bounds(2);
		y_L = y_bounds(1); y_U = y_bounds(2);

		% Incorporate bounds as fixed in a new constraint on K and z
		mccormick_constraint_antecedent = ...
			[ x_L*y + x*y_L - x_L*y_L <= z ] + [ x_U*y + x*y_U - x_U*y_U <= z ] + ...
			[ x_U*y + x*y_L - x_U*y_L >= z ] + [ x_L*y + x*y_U - x_L*y_U >= z ];

		mccormick_constraints = mccormick_constraints + ...
			implies(mccormick_binvars(envelope_idx), mccormick_constraint_antecedent );

	end

	mccormick_constraints = mccormick_constraints + [ sum(mccormick_binvars) == 1 ];

end

function [ z , mccormick_constraints , mccormick_binvars ] = create_mccormick_envelope_for_vector_product( x , y , grid_definition , eta_z_components )
	%create_mccormick_envelope_for_vector_product
	%Description:
	%	Creates a mccormick envelope to represent the product of vector variables x'*y.
	%
	%Inputs:
	%	grid_definition - A struct with members .x_lb, .x_ub, .y_lb, .y_ub, .NumberOfRegionsInX, .NumberOfRegionsInY
	%					  Where every member is a vector of the same dimension.
	%					  Where Number of RegionsInX is a vector of the same dimension of X that determines how many grid cells exist
	%					  along each direction of X's dimensions.

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing / Checking %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if size(x,1) == 1
		x = x';
	end


	if size(y,1) == 1
		y = y';
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = size(x,1);

	if ~exist('eta_z_components')
		eta_z_components = NaN;
	end

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%
	mccormick_constraints = [];
	mccormick_binvars = {};

	for pair_index = 1:dim
		temp_grid_def = struct( ...
			'x_lb',grid_definition.x_lb(pair_index), ...
			'x_ub',grid_definition.x_ub(pair_index), ...
			'y_lb',grid_definition.y_lb(pair_index), ...
			'y_ub',grid_definition.y_ub(pair_index), ...
			'NumberOfRegionsInX',grid_definition.NumberOfRegionsInX(pair_index), ...
			'NumberOfRegionsInY',grid_definition.NumberOfRegionsInY(pair_index) );

		[ z_summand , temp_constraints , mccormick_binvars{pair_index} ] = create_mccormick_envelope_for_scalar_product( x(pair_index) , y(pair_index) , temp_grid_def );

		% Append constraints to mccormick_constraints
		mccormick_constraints = mccormick_constraints + temp_constraints;

		%Add z_summand to z
		if pair_index == 1
			z = z_summand;
		else
			z = z + z_summand;
		end
		
		%If there was a z bound provided, then apply it and add it to the mccormick constraints.
		if ~isnan(eta_z_components)
			mccormick_constraints = mccormick_constraints + [ -eta_z_components <= z_summand <= eta_z_components ];
		end

	end

	% % z is the sum of all the dot product terms
	% z = ones(1,dim)'*z_elt;

end

function [ settings_struct ] = handle_create_mccormick_envelope_arguments( varargin )
	%Description:
	%	Creates the settings struct based on the inputs to the function.

	%% Set Defaults

	settings_struct.eta_z_bounds = NaN;
	settings_struct.debug_flag = 0; %Send no debugging messages.

	%% Algorithms

	if nargin >= 4
		
		argidx = 5;
		while argidx <= nargin
			switch varargin{argidx}
				case 'eta_z_bounds'
					settings_struct.eta_z_bounds = varargin{argidx+1};
					argidx = argidx + 2;
				case 'debug_flag'
					settings_struct.debug_flag = varargin{argidx+1};
					argidx = argidx + 2;
				otherwise
					error('Unexpected input.')
			end
		end

	end

end

function [ x_bounds , y_bounds ] = get_envelope_bounds_helper_lbs( envelope_idx , grid_dims , x_lb , x_ub , y_lb , y_ub )
	%Description:
	%	Creates the "envelope_idx"-th envelope for a bilinear function of x * xy where there
	%	are "num_mccormick_envelopes" number of envelopes created and the bounds of x and y
	%	are provided.
	%
	%Notes:
	%	For a grid broken up into grid_dims = [ num_x1_ticks , num_x2_ticks ], there should be num_x1_ticks * num_x2_ticks = prod(grid_dims)
	%	envelopes.
	%
	%Usage:
	%	[ K_bounds , y_bounds ] = get_envelope_bounds_for_number( evelope_idx , grid_dims , 'LowerBounds', K_lb , K_ub , y_lb , y_ub )

	% Input Processing

	% Constants

	% Algorithm

	x_Dimension = length(x_lb);
	y_Dimension = length(y_lb);

	if (x_Dimension + y_Dimension) ~= length(grid_dims)
		error(['The number of entries in grid_dims does not partition the full space of (x,y) pairs with ' num2str(x_Dimension+y_Dimension) ' dimension.'] )
	end

	x_ranges = x_ub - x_lb;
	y_ranges = y_ub - y_lb;

	num_ticks_x = grid_dims([1:x_Dimension]);
	num_ticks_y = grid_dims([x_Dimension+1:x_Dimension+y_Dimension]);

	xy_tick_indices = cell(1,x_Dimension+y_Dimension);

	[xy_tick_indices{:}] = ind2sub( grid_dims , envelope_idx ); 

	x_tick_indices = xy_tick_indices{[1:x_Dimension]};
	y_tick_indices = xy_tick_indices{[x_Dimension+1:x_Dimension+y_Dimension]};

	% Algorithm

	env_x_lb = x_lb + (x_tick_indices-1).*( x_ranges ./ num_ticks_x);
	env_x_ub = x_lb + x_tick_indices.*( x_ranges ./ num_ticks_x);

	env_y_lb = y_lb + (y_tick_indices-1).*( y_ranges ./ num_ticks_y);
	env_y_ub = y_lb + y_tick_indices.*( y_ranges ./ num_ticks_y);

	x_bounds = [ env_x_lb , env_x_ub ];
	y_bounds = [ env_y_lb , env_y_ub ];

end
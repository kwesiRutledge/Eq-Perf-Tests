function [ varargout ] = find_hypothesis_generating_disturbances( varargin )
	%Description:
	%	Finds the set of all disturbances which 	
	%Usage:
	%	[ PwT ] = lcsas0.find_hypothesis_generating_disturbances( path1 )
	%	[ PwT_array ] = lcsas0.find_hypothesis_generating_disturbances( PathMatrix2 )


	%% Input Processing

	[ lcsas_in , KnowledgePaths , Px0 , Pu , fb_type ] = find_hypothesis_generating_disturbances_input_processing(varargin{:});

	%% Constants

	[ T , num_paths ] = size(KnowledgePaths);

	L_end = KnowledgePaths(end,1);
	card_L_end = L_end.cardinality();
	[ n_x , n_u , n_y , n_w , n_v ] = lcsas_in.Dimensions();

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	if num_paths > 1
		PwT_array = {};
		for path_index = 1:num_paths
			PwT_array{path_index} = lcsas_in.find_hypothesis_generating_disturbances( KnowledgePaths(:,path_index) );
		end
		varargout{1} = PwT_array; %Assign to output.
		return;
	end

	% If we moved past the last "if" statement, then there is only one KnowledgePath
	KnowledgePath = KnowledgePaths;

	%% Create A and b matrices

	max_card = KnowledgePath.find_maximum_cardinality_in_sequence();
	PhiI_dim = n_x * (T+1) + n_u*T + n_w*T*max_card + n_x*max_card;

	[ ibs ] = lcsas_in.internal_behavior_set( KnowledgePath );

	PwT = [];
	for word_idx = 1:L_end.cardinality()
		PwT = [ PwT ; ibs.projection(n_x*(T+1) + n_u*(T) + n_w*(T)*(word_idx-1) + [1:n_w*(T)]) ];
	end

	varargout{1} = PwT;

end

function [ lcsas_out , KnowledgePaths , Px0 , Pu , fb_type ] = find_hypothesis_generating_disturbances_input_processing(varargin)
	%Description:
	%	Manages the inputs appropriately.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%
	%% Input Checking %%
	%%%%%%%%%%%%%%%%%%%%

	if nargin < 2
		error(['The function find_hypothesis_generating_disturbances() requires at least two input arguments. Received ' num2str(nargin) '.'])
	end

	% Assigning required values

	lcsas_out = varargin{1};
	KnowledgePaths = varargin{2};
	% Px0 = varargin{4};
	% Pu = varargin{5};

	if ~isa(lcsas_out,'LCSAS')
		error(['The first input is not of type LCSAS. Received type ' class(lcsas_out) '.' ] )
	end

	if isempty(lcsas_out.X0)
		error(['get_feasible_combinations_of_beliefs() requires nonempty choice of X0 in lcsas object.'])
	end

	if ~isa(lcsas_out.X0,'Polyhedron')
		error(['get_feasible_combinations_of_beliefs() requires X0 to be a Polyhedron.'])
	end

	Px0 = lcsas_out.X0;

	if isempty(lcsas_out.U)
		error(['get_feasible_combinations_of_beliefs() requires nonempty choice of U in lcsas object.'])
	end

	if ~isa(lcsas_out.U,'Polyhedron')
		error(['get_feasible_combinations_of_beliefs() requires U to be a Polyhedron.'])
	end

	Pu = lcsas_out.U;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Assign Default Values for All Other Things %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	fb_type = 'state'; %State feedback option.

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Process Any Extra Inputs %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if nargin > 5
		argidx = 6;
		while argidx <= nargin
			switch varargin{argidx}
				case 'Feedback Type'
					fb_type = varargin{argidx+1};
					argidx = argidx + 2;
				otherwise
					error(['Unexpected input to find_hypothesis_generating_disturbances: ' varargin{argidx} ])
			end
		end
	end
end
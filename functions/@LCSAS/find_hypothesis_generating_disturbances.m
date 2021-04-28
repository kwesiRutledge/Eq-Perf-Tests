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

	%$ Create A and b matrices

	PhiI_A  = []; PhiI_b  = [];
	PhiI_Ae = []; PhiI_be = [];

	max_card = KnowledgePaths.find_maximum_cardinality_in_sequence();

	PhiI_dim = n_x * (T+1) + n_u*T + n_w*T*max_card + n_x*max_card;

	for t = 1:T

		KnowledgeAt_t = KnowledgePaths(t,1);

		Phi_t = lcsas_in.consistent_set( t , KnowledgeAt_t , Pu , Px0 );

		WT_vector = zeros(max_card,1);
		for word_index = 1:KnowledgeAt_t.cardinality()
			temp_word = KnowledgeAt_t.words{word_index};
			[ tf ] =  
		end
		WT_vector()

		WT_factor = kron(  )

		A_Prefactor = [ ...
			eye( n_x * (t+1) ), zeros( n_x*(t+1) , PhiI_dim - n_x * (t+1) ) ;
			zeros( n_u*(t) , n_x*(T+1) ), eye(n_u*t) , zeros( n_u*(t) , PhiI_dim - n_x*(T+1) - n_u*t ) ;
			zeros( n_w*(t)*card_L_end , n_x*(T+1) + n_u*T ) , kron( eye(card_L_end) , [ eye(n_w*t) , zeros(n_w*t,n_w*(T-t)) ] ) , zeros(n_w*t*card_L_end, PhiI_dim - n_x*(T+1) - n_u*T - n_w*T*card_L_end) ;
			zeros( n_x*card_L_end, PhiI_dim - n_x*card_L_end), eye(n_x*card_L_end) ...
			];

		PhiI_A = [ 	PhiI_A ;
					Phi_t.A * A_Prefactor ];

	end

	[ ~ , PhiI ] = lcsas_in.consistent_set( T , L_end , Pu , Px0 );

	PwT = [];
	for word_idx = 1:L_end.cardinality()
		PwT = [ PwT ; PhiI.projection(n_x*(T+1) + n_u*(T) + n_w*T*(word_idx-1) + [1:n_w*(T)]) ];
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
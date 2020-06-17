function [ initial_nodes ] = get_initial_beliefnodes( varargin )
	%Desctiption:
	%	This function uses the definition of the current belief nodes to
	%
	%Usage:
	%	[ initial_nodes ] = bg.get_initial_beliefnodes()
	%	[ initial_nodes ] = bg.get_initial_beliefnodes(P_x0)
	%	[ initial_nodes ] = bg.get_initial_beliefnodes(P_x0, 'verbosity' , 0)

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%
	
	bg = varargin{1};
	arg_idx = 2;

	%Check to see if the initial state set is embedded in the belief graph's LCSAS member.
	lcsas_in = bg.lcsas;
	if isprop(lcsas_in,'X0')

		%Check to see if LCSAS object has a polyhedron value in X0
		x0_undefined = false;
		if ~isa(lcsas_in.X0,'Polyhedron')
			x0_undefined = true;
		end

		if x0_undefined & (nargin == 1)
			error('X0 was not defined in LCSAS and needs to be provided to get_initial_beliefnodes().')
		end

		if x0_undefined & ~isa(varargin{2},'Polyhedron')
			error('X0 was not defined by second argument to get_initial_beliefnodes().')
		end

		%Officially Define the value of X0
		if x0_undefined
			X0 = varargin{2};
			arg_idx = arg_idx + 1;
		else
			X0 = lcsas_in.X0;

			if isa(varargin{2},'Polyhedron')
				disp('Ignoring input X0.')
				arg_idx = arg_idx + 1;
			end

		end
	end

	disp(['arg_idx = ',num2str(arg_idx)])

	while arg_idx <= nargin
		switch varargin{arg_idx}
			case 'verbosity'
				verbosity = varargin{arg_idx+1};
				arg_idx = arg_idx + 2;
			case 'OverrideProjFlag'
				proj_flag = varargin{arg_idx+1};
				arg_idx = arg_idx+2;
			otherwise
				error(['Unexpected input to get_initial_beliefnodes: ' varargin{arg_idx} ])
		end
	end

	if ~exist('verbosity')
		verbosity = 0;
	end

	if ~exist('proj_flag')
		proj_flag = bg.UsedProjection;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	L = bg.ModeLanguage;
	L_powerset = L.powerset();

	% Define the Powerset of word_idx values.
	word_idcs = [1:L.cardinality()];
	word_idx_powerset = {};
	for k = 1:L.cardinality()
		%Consider all combinations of k elements from L
		temp_combs = nchoosek([1:L.cardinality()],k);
		for elt_idx = 1:size(temp_combs,1)
			%Add each element of temp_combs to word_idx_powerset
			word_idx_powerset{length(word_idx_powerset)+1} = temp_combs(elt_idx,:);
		end
	end

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%% Collect all Y Sets
	[initial_eb_sets,initial_ib_sets] = lcsas_in.get_consistency_sets_for_language(0,L,Polyhedron(),X0, ...
																					'use_proj', proj_flag, ...
																					'ConsistencySetVersion', bg.ConsistencySetVersion );


	if proj_flag
		eb_sets = initial_eb_sets;

		for powerset_idx = (L.cardinality()+1):length(word_idx_powerset)
			powerset_elt = word_idx_powerset{powerset_idx};
			powerset_elt_L = L_powerset(powerset_idx);

			%For each word, create the output set of the system.
			temp_sigma = powerset_elt_L.words{1};
			[ ~ , L_idx_of_sigma ] = L.contains(temp_sigma);

			%Compute Sets
			eb_sets = [eb_sets, eb_sets(L_idx_of_sigma)];
			for sigma_idx = 2:powerset_elt_L.cardinality()
				%For each word, create the output set of the system.
				temp_sigma = powerset_elt_L.words{sigma_idx};
				[ ~ , L_idx_of_sigma ] = L.contains(temp_sigma);

				eb_sets(end) = eb_sets(end).intersect( eb_sets(L_idx_of_sigma) );
			end

		end
	end

	ib_sets = bg.get_all_consistent_internal_behavior_sets( initial_ib_sets );

	%% Find which Y_Sets contain others
	if proj_flag
		containment_matrix = false(length(word_idx_powerset));
		for word_idx = 1:length(word_idx_powerset)
			%word at word_idx always contains itself.
			containment_matrix(word_idx,word_idx) = true;
		end

		for x_idx = 1:size(containment_matrix,1)
			potential_y_idcs = [1:size(containment_matrix,2)];
			potential_y_idcs = potential_y_idcs( potential_y_idcs ~= x_idx );
			for y_idx = potential_y_idcs
				%Consider the temporary combination
				%disp(['x_idx = ' num2str(x_idx) ', y_idx = ' num2str(y_idx) ])

				% Observe if Y_Set of X is contained by the Y_Set of Y
				ObservationSetX = eb_sets(x_idx);
				ObservationSetY = eb_sets(y_idx);
				
				containment_matrix(x_idx,y_idx) = (ObservationSetX <= ObservationSetY);
			end
		end
	else
		containment_matrix = bg.internal_behavior_sets2containment_mat( ib_sets , 'verbosity' , verbosity );
	end

	% containment_matrix

	%% Create Nodes Based on containment_matrix and whether set is nonempty %%

	if proj_flag
		[ ~ , empty_set_flags ] = bg.find_empty_observation_polyhedra( eb_sets );
	else
		[ ~ , empty_set_flags ] = bg.find_empty_observation_polyhedra( ib_sets );
	end

	[ ~ , observation_set_is_observable ] = bg.containment_mat2observable_combos( containment_matrix );

	valid_powerset_idcs = [1:length(L_powerset)];
	valid_powerset_idcs = valid_powerset_idcs( (~empty_set_flags) & observation_set_is_observable );

	%Create Belief Nodes
	initial_nodes = [];
	for bn_idx = 1:length(valid_powerset_idcs)
		temp_powerset_idx = valid_powerset_idcs(bn_idx);
		temp_comb = word_idx_powerset{temp_powerset_idx};

		%Create Language for temp_comb
		temp_lang = L_powerset(valid_powerset_idcs(bn_idx));

		%Create Consistency Set for temp_comb
		if proj_flag
			temp_cs = eb_sets(valid_powerset_idcs(bn_idx));
			initial_nodes = [initial_nodes, BeliefNode( temp_lang , 0 , temp_cs )];
		else
			temp_ib_set = ib_sets(valid_powerset_idcs(bn_idx));
			initial_nodes = [initial_nodes, BeliefNode( temp_lang , 0 , ...
														'FullTrajectorySet', temp_ib_set )];
		end
	end

end
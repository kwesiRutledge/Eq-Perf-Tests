function [ initial_nodes ] = get_initial_beliefnodes( varargin )
	%Desctiption:
	%	This function uses the definition of the current belief nodes to
	%
	%Usage:
	%	[ initial_nodes ] = bg.get_initial_beliefnodes()
	%	[ initial_nodes ] = bg.get_initial_beliefnodes('X0',X0)
	%	[ initial_nodes ] = bg.get_initial_beliefnodes('X0',X0, 'verbosity' , 0)
	%
	%Todo:
	%	- Speed up the exhaustive search for all Y0's.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%
	
	bg = varargin{1};
	arg_idx = 2;

	while arg_idx <= nargin
		switch varargin{arg_idx}
			case 'X0'
				X0 = varargin{arg_idx+1};
				arg_idx = arg_idx + 2;
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

	lcsas_in = bg.lcsas;

	if isempty(lcsas_in.X0) && ~exist('X0')
		error('X0 was not defined in LCSAS and needs to be provided to get_initial_beliefnodes().')
	end

	%%%%%%%%%%%%%%%%%%%%
	%% Default Values %%
	%%%%%%%%%%%%%%%%%%%%

	if ~exist('verbosity')
		verbosity = 0;
	end

	if ~exist('proj_flag')
		proj_flag = bg.UsedProjection;
	end

	if ~exist('X0')
		X0 = lcsas_in.X0;
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

	initial_nodes = [];

	switch bg.FeedbackMethod
	case 'state'
		% If we are using state feedback and there is only one X0 set, then there should be a single output BeliefNode.
		initial_nodes = BeliefNode( L , 0 , 'ConsistencySet' , X0 );
	case 'output'
		Y0_list = [];
		for word_index = 1:L.cardinality()
			word_i = L.words{word_index};
			temp_first_symbol = word_i(1);

			Y0_list = [ Y0_list ; lcsas_in.Dyn( temp_first_symbol ).C * X0 + lcsas_in.Dyn( temp_first_symbol ).P_v ];
		end

		% Determine Which Combinations of Words Intersect.
		nonempty_combinations = [];
		Y0_combinations = [];
		for combination_index = 1:length(word_idx_powerset)
			temp_combination = word_idx_powerset{combination_index};

			%Intersect all Elements In The 
			temp_Y0 = Y0_list(temp_combination(1));
			for comb_vector_index = 1:length(temp_combination)
				temp_Y0 = temp_Y0.intersect( Y0_list(temp_combination(comb_vector_index)) );
			end

			%Save Results
			Y0_combinations = [ Y0_combinations ; temp_Y0 ];
			nonempty_combinations = [ nonempty_combinations ; ~temp_Y0.isEmptySet ];
		end

		% Determine which combinations Contain Others
		last_nonempty_comb_index = find(nonempty_combinations,1,'last');
		nonempty_combinations_indices = find(nonempty_combinations);
		
		% Initialize Loop with an Index to Start with and the valid combinations.
		subset_checker_index = last_nonempty_comb_index;
		viable_combinations = nonempty_combinations;

        %Ignore all combinations that have a lower combination index
		%than the current one.
		% By searching in this order, we can ignore combinations that have
		% higher cardinality than the current one.
        
		while subset_checker_index > 1

			% Check to see if any of the nonempty combinations are contained by the "larger"
			% ones.

			for combination_index = 1:subset_checker_index-1

				if ~viable_combinations(combination_index)
					%IF this combination was empty or if this combination was previously
					%shown to be unnecessary, then continue.
					continue;
                end

				if (Y0_combinations(subset_checker_index) <= Y0_combinations(combination_index)) && (Y0_combinations(subset_checker_index) >= Y0_combinations(combination_index))
					%If the two combinations are equal, then
					%ignore the one of smaller cardinality (combination_index)
					viable_combinations(combination_index) = false;
				end

			end

			%Get the last viable combination in the list that is "before" the current index
			subset_checker_index = find(find( viable_combinations ) < subset_checker_index,1,'last');

		end

		%Finally create the nodes
		for combination_index = 1:length(viable_combinations)
			
			if viable_combinations(combination_index)
				temp_combination = word_idx_powerset{combination_index};
				initial_nodes = [ initial_nodes , BeliefNode( Language(L.words{ temp_combination }) , 0 , 'ConsistencySet' , Y0_combinations(combination_index) ) ];
			end
        end
        1;

	otherwise
		error(['Unexpected FeedbackMethod for BeliefGraph: ' bg.FeedbackMethod ])
	end

end
function [ possibleSubsetsOfPaths , possible_choices , choices_as_binary_flags ] = get_feasible_combinations_of_beliefs( lcsas_in , KnowledgePaths )
	%Description:
	%	This function attempts to identify all possible combinations of
	%	knowledge sequences where the knowledge sequences come from LK{end}.
	%
	%Inputs:
	%	lcsas_in: The LCSAS object for which these beliefs will be with respect to.
	%
	%Usage:
	%	[ possibleSubsetsOfPaths , possible_choices , choices_as_binary_flags ] = lcsas_in.get_feasible_combinations_of_beliefs( KnowledgePaths )

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing/Checking %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if isempty(lcsas_in.X0)
		error(['get_feasible_combinations_of_beliefs() requires nonempty choice of X0 in lcsas object.'])
	end

	if ~isa(lcsas_in.X0,'Polyhedron')
		error(['get_feasible_combinations_of_beliefs() requires X0 to be a Polyhedron.'])
	end

	if isempty(lcsas_in.U)
		error(['get_feasible_combinations_of_beliefs() requires nonempty choice of U in lcsas object.'])
	end

	if ~isa(lcsas_in.U,'Polyhedron')
		error(['get_feasible_combinations_of_beliefs() requires U to be a Polyhedron.'])
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ TimeHorizon , numPaths ] = size(KnowledgePaths);

	Px0 = lcsas_in.X0;
	Pu  = lcsas_in.U;
    
	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%% Construct Disturbance Sets Associated with Each Element of LK
	LK_w = lcsas_in.find_hypothesis_generating_disturbances( KnowledgePaths );
	for LK_sequence_index = 1:numPaths
		LK_w{LK_sequence_index}.outerApprox;
	end

	%% Identify If Any Of The Disturbance Sets are Empty (=>Consistency Sets are also empty)
	[ LK_T_indices_trimmed , empty_LK_w_indices ] = RemovePathsWithEmptyConsistencySets( LK_w );

	%% Identify All Possible Combinations of Matching Behaviors That Might Exist	
	possible_choices = {};
	for cardinality = 1:length(LK_T_indices_trimmed)
		possible_choices{cardinality} = nchoosek(LK_T_indices_trimmed,cardinality);
	end

	% Remove All Choices that:
	% - Do Not Respect/Consider All Modes
	% - Do Not Have a Covering Set of Disturbances

	possible_choices = RemoveCombinationsThatDoNotConsiderAllModes( possible_choices , lcsas_in , KnowledgePaths );
	possible_choices = lcsas_in.RemoveCombinationsThatDoNotHaveCoveringLKw( possible_choices , KnowledgePaths , LK_w );

	% Transform possible_choices to a list of binary vectors
	choices_as_binary_flags = {};
	all_false_vec = false(numPaths,1);
	for cardinality = 1:length(possible_choices)

		choices_with_cardinality = possible_choices{cardinality};
		
		for choice_idx = 1:size(possible_choices{cardinality},1)
			temp_choice = choices_with_cardinality(choice_idx,:);
			choices_as_binary_flags{end+1} = all_false_vec;
			choices_as_binary_flags{end}(temp_choice) = true;
		end

	end

	% Translate this choice cell array into a cell array of sequences.
	possibleSubsetsOfPaths = {};
	for temp_length = 1:length(LK_T_indices_trimmed)
		pathSubsetsAtTL = {};
		for choice_idx = 1:size(possible_choices{temp_length},1)
			temp_choice = possible_choices{temp_length}(choice_idx,:);
			pathSubsetsAtTL{choice_idx} = KnowledgePaths(:,temp_choice);
			possibleSubsetsOfPaths{end+1} = KnowledgePaths(:,temp_choice);
		end
	end

end

function [ nonempty_LK_w_indices, empty_LK_w_indices ] = RemovePathsWithEmptyConsistencySets( LK_w )

	%% Variables
	numPaths = length(LK_w);

	%% Algorithm

	LK_w_set_is_empty = [];
	for LK_sequence_index = 1:numPaths
		LK_w_set_is_empty(LK_sequence_index) = LK_w{LK_sequence_index}(1).isEmptySet;
	end
	empty_LK_w_indices = find(LK_w_set_is_empty);

	% Results

	nonempty_LK_w_indices = setdiff( [1:numPaths] , empty_LK_w_indices );

end

function [ possible_choices_out ] = RemoveCombinationsThatDoNotConsiderAllModes( possible_choices , lcsas_in , KnowledgePaths )
	%Description:
	%	
	%Example:
	%	For a lcsas with L = {1,2,3}, it should be the case that an allowable combination has "leaves" that contain 1 and 2 and 3.

	%% Variables

	L = lcsas_in.L;

	%% Algorithm

	possible_choices_out = {};
	for cardinality = 1:length(possible_choices)

		unfiltered_choices_with_cardinality = possible_choices{cardinality};
		
		possible_choices_out{cardinality} = [];
		for choice_idx = 1:size(possible_choices{cardinality},1)
			temp_choice_indices = unfiltered_choices_with_cardinality(choice_idx,:);
			temp_choice = KnowledgePaths(:,temp_choice_indices);

			% Only Add Combinations to the final list (possible_choices_out) if the combinations contain all words
			% from L.

			last_langs = temp_choice(end,:);
			if length(last_langs) == 1
				if L == last_langs
					possible_choices_out{cardinality} = [ possible_choices_out{cardinality} ; temp_choice_indices ];
				end
			else
				ll_union = last_langs(1);
				ll_union = ll_union.union([last_langs(2:end)]);
				if L == ll_union
					possible_choices_out{cardinality} = [ possible_choices_out{cardinality} ; temp_choice_indices ];
				end
			end

		end

	end


end
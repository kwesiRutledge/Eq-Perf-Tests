function [ possibleSubsetsOfPaths , possible_choices , choices_as_binary_flags ] = get_feasible_combinations_of_beliefs( varargin )
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

	[ lcsas_in , KnowledgePaths , ibs_sets , gfcob_settings ] = ip_get_feasible_combinations_of_beliefs( varargin{:} );

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ TimeHorizon , numPaths ] = size(KnowledgePaths);

	Px0 = lcsas_in.X0;
	Pu  = lcsas_in.U;

	L = lcsas_in.L;
    
	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	if gfcob_settings.verbosity > 0
		disp('Starting get_feasible_combinations_of_beliefs()...')
	end

	%% Construct Disturbance Sets Associated with Each Element of LK
	if isempty(ibs_sets)
		% LK_w = lcsas_in.find_hypothesis_generating_disturbances( KnowledgePaths );
		% for LK_sequence_index = 1:numPaths
		% 	LK_w{LK_sequence_index}.outerApprox;
		% end
		for knowledge_path_index = 1:numPaths
			temp_path = KnowledgePaths(:,knowledge_path_index);
			temp_ibs = InternalBehaviorSet(lcsas_in,temp_path);
			ibs_sets(knowledge_path_index) = temp_ibs;
		end

		if gfcob_settings.verbosity > 0
			disp('- Created All Internal Behavior Sets in get_feasible_combinations_of_beliefs().')
		end

	end

	%% Identify If Any Of The Disturbance Sets are Empty (=>Consistency Sets are also empty)
	empty_flags = ibs_sets.IsEmpty();
	LK_T_indices_trimmed = setdiff( [1:length(ibs_sets)] , find(empty_flags) );

	if gfcob_settings.verbosity > 0
		disp('- Removed all paths with empty internal behavior sets (avoids projection).')
	end


	%% Identify All Possible Combinations of Matching Behaviors That Might Exist	
	possible_choices = {};
	for cardinality = 1:length(LK_T_indices_trimmed)
		possible_choices{cardinality} = nchoosek(LK_T_indices_trimmed,cardinality);
	end

	% Remove All Choices that:
	% - Do Not Respect/Consider All Modes
	% - Do Not Have a Covering Set of Disturbances

	possible_choices = RemoveCombinationsThatDoNotConsiderAllModes( possible_choices , lcsas_in , KnowledgePaths );
	if gfcob_settings.verbosity > 0
		disp('- Removed all combinations that did not consider all modes of the system lcsas_in.')
	end

	possible_choices = lcsas_in.RemoveCombinationsThatDoNotHaveCoveringBehaviorSets( possible_choices , KnowledgePaths , ibs_sets );
	if gfcob_settings.verbosity > 0
		disp('- Removed all combinations that did not consider all modes of the system lcsas_in.')
	end

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

function [ lcsas_in , KnowledgePaths , InternalBehaviorSets , gfcob_settings ] = ip_get_feasible_combinations_of_beliefs( varargin )
	%Description:
	%	Processes the inputs to the function get_feasible_combinations_of_beliefs().
	
	%% Extract the lcsas and KnowledgePaths variables

	lcsas_in = varargin{1};
	KnowledgePaths = varargin{2};

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

	%% Get Extra Variables
	
	InternalBehaviorSets = InternalBehaviorSet.empty;
	gfcob_settings = struct('verbosity',0);

	arg_index = 3;
	while arg_index <= nargin
		switch varargin{arg_index}
			case 'InternalBehaviorSets'
				InternalBehaviorSets = varargin{arg_index+1};
				arg_index = arg_index + 2;
			case 'verbosity'
				gfcob_settings.verbosity = varargin{arg_index+1};
				arg_index = arg_index + 2;
			otherwise
				error(['Unexpected input to get_feasible_combinations_of_beliefs: ' varargin{arg_index} ])
		end
	end



end

function [ nonempty_ibs_indices, empty_ibs_indices ] = RemovePathsWithEmptyConsistencySets( InternalBehaviorSets )
	%Description:
	%	This function identifies which of the InternalBehaviorSets are empty.
	%	If a given path's internal behavior set is empty, then we will mark that set appropriately.

	%% Variables
	numPaths = length(InternalBehaviorSets);

	%% Algorithm

	ibs_is_empty = [];
	for LK_sequence_index = 1:numPaths
		ibs_is_empty(LK_sequence_index) = InternalBehaviorSets{LK_sequence_index}(1).isEmptySet;
	end
	empty_ibs_indices = find(ibs_is_empty);

	% Results

	nonempty_ibs_indices = setdiff( [1:numPaths] , empty_ibs_indices );

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
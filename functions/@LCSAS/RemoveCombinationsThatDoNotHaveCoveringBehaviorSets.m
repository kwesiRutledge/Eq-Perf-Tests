function [ combinations_out ] = RemoveCombinationsThatDoNotHaveCoveringBehaviorSets( lcsas_in , combinations_in , KnowledgePaths , InternalBehaviorSets )
	%Description:
	%	Identifies the combinations in combinations_in which do or do not properly cover the disturbance set as they should.

	%% Input Processing
	if ~isa('InternalBehaviorSet',InternalBehaviorSets)
		error(['The input InternalBehaviorSets must be of class InternalBehaviorSet; instead it is of type ' class(InternalBehaviorSets)])
	end

	%% Variables

	L 			= lcsas_in.L;
	TimeHorizon = size(KnowledgePaths,1);

	%% Algorithm

	% Create the Disturbance Sequences Expected of Each Word
	WT = {};
	for word_index = 1:L.cardinality()
		temp_word = L.words{word_index};
		WT{word_index} = 1;
		for time_index = 0:TimeHorizon-1
			temp_symbol = temp_word(time_index+1);
			WT{word_index} = WT{word_index} * lcsas_in.Dyn(temp_symbol).P_w;
		end
	end

	% Verify that Each Choice "Covers" Each WT
	combinations_out = {};
	for cardinality = 1:length(combinations_in)

		unfiltered_choices_with_cardinality = combinations_in{cardinality};
		
		combinations_out{cardinality} = [];
		for choice_idx = 1:size(combinations_in{cardinality},1)
			temp_choice_indices = unfiltered_choices_with_cardinality(choice_idx,:);
			temp_choice = KnowledgePaths(:,temp_choice_indices);

			% Find All Active Paths Containing word index 1, then 2, then 3 ...
			[ path_subset_list , index_subset_list ] = temp_choice.find_paths_with_ending_from( L ); % Result stored in path_subset_list{1}, path_subset_list{2}, ...

			% % Find the Disturbance Sets associated with each element of temp_choice.
			% for subset_elt_index = 1:size(temp_choice,2)
			% 	temp_path = temp_choice(:,subset_elt_index);
			% 	LK_w{subset_elt_index} = lcsas0.find_hypothesis_generating_disturbances( temp_path(end) );
			% end

			% Only Add Combinations to the final list (combinations_out) if the combinations contain all words
			% from L.

			all_words_covered = true;

			% Check each word in L.
			for word_index = 1 : L.cardinality()
				
				paths_in_subset_containing_word_i = path_subset_list{word_index};
				indices_in_subset_containing_word_i = index_subset_list{word_index};

				% Check for covering

				% path_subset_polyunion = PolyUnion( [ LK_w{ temp_choice_indices(indices_in_subset_containing_word_i) } ] );
				
				% all_words_covered = all_words_covered && p_subseteq_pu( WT{word_index} , path_subset_polyunion );

			end

			% If all words are covered, then add the choice.
			if all_words_covered
				combinations_out{cardinality} = [ combinations_out{cardinality} ; temp_choice_indices ];
			end

		end

	end


end

% function tf = InternalBehaviorSetDisturbancesCoverPolyhedron( t , InternalBehaviorSets , TargetPolytope1 )
% 	%Description:
% 	%	Identifies if there exists an internal behavior set in InternalBehaviorSets such that the disturbances for those sets
% 	%	include all points from the target polytope, TargetPolytope1.

% 	%% Constants

% 	R_T = [ zeros( n_w*T ,  ) ]

% 	%% Algorithm



% end
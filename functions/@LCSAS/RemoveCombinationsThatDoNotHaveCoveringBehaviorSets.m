function [ combinations_out ] = RemoveCombinationsThatDoNotHaveCoveringBehaviorSets( lcsas_in , combinations_in , KnowledgePaths , InternalBehaviorSets )
	%Description:
	%	Identifies the combinations in combinations_in which do or do not properly cover the disturbance set as they should.

	%% Input Processing
	if ~isa(InternalBehaviorSets,'InternalBehaviorSet')
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

			% If all words are covered, then add the choice.
			if temp_choice.EvaluateIfAllWordBehaviorsAreCovered( lcsas_in )
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
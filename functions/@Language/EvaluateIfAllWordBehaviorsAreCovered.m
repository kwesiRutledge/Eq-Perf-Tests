function [tf] = EvaluateIfAllWordBehaviorsAreCovered( ChoiceOfPaths , System )
	%Description:
	%	Determines if the set of paths ChoiceOfPaths covers the behaviors that can be ellicited by
	%	each word from System.
	%
	%Usage:
	%	tf = Paths.EvaluateIfAllWordBehaviorsAreCovered( System )

	%% Input Processing %%

	if ~isa(System,'LCSAS')
		error(['System input to PathSubsetCoversBehaviorsFromWord() is of class ' class(System) ' but expected class LCSAS.' ])
	end

	%% Constants

	L = System.L;

	%% Algorithm

	% Find All Active Paths Containing word index 1, then 2, then 3 ...
	[ path_subset_list , index_subset_list ] = ChoiceOfPaths.find_paths_with_ending_from( L ); % Result stored in path_subset_list{1}, path_subset_list{2}, ...

	tf = true;

	% Check each word in L.
	for word_index = 1 : L.cardinality()

		temp_word = L.words{word_index};
		
		paths_in_subset_containing_word_i = path_subset_list{word_index};
		indices_in_subset_containing_word_i = index_subset_list{word_index};

		% Check for covering
		tf = tf && paths_in_subset_containing_word_i.PathSubsetCoversBehaviorsFromWord( temp_word , System );

	end

end
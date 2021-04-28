function [ paths_subset_containing_word , index_subset_containing_word ] = find_paths_with_ending_from( obj , L )
	%Description:
	%	For each word w in L, finds the subset of paths in obj (aka KnowledgePaths) which contains w.
	%
	%Usage:
	%	[ path_subset_list , index_subset_list ] = find_paths_with_ending_from( obj , L )

	%% Variables

	KnowledgePaths = obj;

	paths_subset_containing_word = {};
	index_subset_containing_word = {};

	%% Algorithms

	for word_index = 1:L.cardinality()

		word_i = L.words{word_index};
		[ paths_subset_containing_word{word_index} , index_subset_containing_word{word_index} ] = KnowledgePaths.FindPathsWhoseEndingContains(word_i);

	end


end
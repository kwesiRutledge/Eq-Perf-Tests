function [ PathsOut , PathIndices ] = FindPathsWhoseEndingContains(obj,target_word)
	%Description:
	%	This function finds which paths lead to a final knowledge set (i.e. Language) that contain the
	%	target word. Each column of obj (aka KnowledgePaths) is a sequence of Language objects and we
	%	need to test whether or not the final Language object contains target_word.
	%	KnowledgePaths (aka obj) is meant to be a matrix where each column represents a "path"
	%	in a belief/knowledge tree. Each element of the matrix is a Language object.
	%
	%Usage:
	%	[ P_out , P_indices ] = KnowledgePaths.FindPathsWhoseEndingContains(target_word)
	

	%% Variables

	KnowledgePaths = obj;

	[ T , numPaths ] = size(KnowledgePaths);

	%% Algorithm

	PathsOut = [];
	PathIndices = [];

	for path_index = 1:numPaths
		path_i = KnowledgePaths(:,path_index);

		final_knowledge = path_i(end);

		if final_knowledge.contains(target_word)
			PathsOut = [PathsOut,path_i];
			PathIndices = [PathIndices,path_index];
		end

	end

end
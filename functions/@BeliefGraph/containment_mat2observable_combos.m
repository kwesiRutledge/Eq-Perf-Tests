function [ observable_set_idcs , obsv_set_is_observable ] = containment_mat2observable_combos( beliefgraph_in , contain_matrix )
	%Description:
	%	

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n_polys = size(contain_matrix,1);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%% Observe which elements of the powerset are "empty"
	% i.e. each element word_idx_powerset{i} of the powerset corresponds to an intersection
	%      of each element

	obsv_set_is_contained = false(n_polys,1);
	observable_set_idcs = [];

	for x_idx = 1:n_polys
		%This BeliefNode or Observation Set is not observed if it is contained
		%by another node of LARGER cardinality.
		for y_idx = (x_idx+1):n_polys %powerset_idx2 = powerset_idx+1:length(word_idx_powerset)
			if contain_matrix(x_idx,y_idx)
				obsv_set_is_contained(x_idx) = true;
			end
		end

		if obsv_set_is_contained(x_idx)
			observable_set_idcs = [observable_set_idcs ; x_idx];
		end

	end

	obsv_set_is_observable = ~obsv_set_is_contained;

end
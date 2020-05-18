function [ empty_set_idcs , empty_set_flags ] = find_empty_observation_polyhedra( beliefgraph_in , observation_polys_in )
	%find_empty_observation_polyhedra()
	%Description:
	%
	%Notes:
	%	find_empty_observation_polyhedra( beliefgraph_in , observation_polys_in ) is equivalent to using
	%	find_empty_observation_polyhedra( beliefgraph_in , observation_polys_in , 'Use BG.ModeLanguage' )
	%
	%Usage:
	%	[empty_set_flags , empty_set_idcs] = find_empty_observation_polyhedra( beliefgraph_in , observation_polys_in , 'L_in' , L_in )
	%	[empty_set_flags , empty_set_idcs] = find_empty_observation_polyhedra( beliefgraph_in , observation_polys_in )
	%	[empty_set_flags , empty_set_idcs] = find_empty_observation_polyhedra( beliefgraph_in , observation_polys_in , 'Use BG.ModeLanguage' )

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	% beliefgraph_in = varargin{1};
	% observation_polys_in = varargin{2};

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	empty_set_flags = false(length(observation_polys_in),1);
	empty_set_idcs = [];

	for set_idx = 1:length(observation_polys_in)
		temp_obsv_set = observation_polys_in(set_idx);

		empty_set_flags(set_idx) = temp_obsv_set.isEmptySet;
		if temp_obsv_set.isEmptySet
			empty_set_idcs = [empty_set_idcs ; set_idx ];
		end
	end


end
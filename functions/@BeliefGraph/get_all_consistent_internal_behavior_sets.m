function [ extended_internal_behavior_sets ] = get_all_consistent_internal_behavior_sets( belief_graph_in , initial_internal_behavior_sets )
	% get_consistent_internal_behavior_sets.m
	%Description:
	%	Computes the consistent internal behavior sets for combinations of the sets in initial_internal_behavior_sets.
	%	Another interpretation of this is as follows:
	%		- Let iibs = initial_internal_behavior_sets
	%		- iibs(1) is the internal behavior set for word 1, iibs(2) is the internal behavior set for word 2, etc. 	
	%		- This script computes the sets for combinations such as the internal behavior set consistent with word 1 AND 2
	%		  this can be done by creating a higher dimensional polyhedron.
	%
	%Usage:
	%	extended_internal_behavior_sets = bg.get_all_consistent_internal_behavior_sets( initial_internal_behavior_sets )

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dummy_words = {};
	for word_idx = 1:length(initial_internal_behavior_sets)
		dummy_words{word_idx} = word_idx;
	end

	dummy_L = Language(dummy_words);
	[~, dummy_L_index_powerset] = dummy_L.powerset();

	dim = initial_internal_behavior_sets(1).Dim;

	fb_method = belief_graph_in.FeedbackMethod;

	%Find the value of time horizon t
	lcsas_in = belief_graph_in.lcsas;
	n_x = size(lcsas_in.Dyn(1).A,1);
	n_u = size(lcsas_in.Dyn(1).B,2);
	n_w = size(lcsas_in.Dyn(1).B_w,2);
	if strcmp(fb_method,'output')
		n_y = size(lcsas_in.Dyn(1).C,1);
		n_v = size(lcsas_in.Dyn(1).C_v,2);

		t = (dim-n_y-n_v-2*n_x)/(n_y+n_u+n_w+n_v+n_x);
		
	elseif strcmp(fb_method,'state')
		t = (dim - 2*n_x)/(n_x+n_u+n_w);
	else
		error('Unrecognized feedback method.')
	end
		


	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Pipe the initial_internal_behavior_sets into extended_internal_behavior_sets
	extended_internal_behavior_sets = initial_internal_behavior_sets;

	%Begin combining sets to make the extended internal_behavior_sets
	for powerset_idx = (length(initial_internal_behavior_sets)+1):length(dummy_L_index_powerset)
		temp_combination = dummy_L_index_powerset{powerset_idx};

		temp_polytope = initial_internal_behavior_sets(temp_combination(1));

		for tc_idx = 2:length(temp_combination)
			temp_polytope = temp_polytope*initial_internal_behavior_sets(temp_combination(tc_idx));
		end

		%Add an Equality constraint between each of the combined polytopes.
		temp_comb_length = length(temp_combination);
		temp_Ae = [];
		for tc_idx = 1:(length(temp_combination)-1)
			if strcmp(fb_method,'output')
				eb_size = n_y*(t+1)+n_u*t; %Size of the external behavior vector
				temp_Ae = [	temp_Ae;
						[zeros(eb_size,(tc_idx-1)*dim),[eye(eb_size),zeros(eb_size,dim-eb_size)],-[eye(eb_size),zeros(eb_size,dim-eb_size)],zeros(eb_size, (temp_comb_length - tc_idx - 1)*dim)] ];
			else
				error('This function is not fully implemented for state feedback yet.')
			end
		end

		temp_be = zeros(size(temp_Ae,1),1);

		extended_internal_behavior_sets(powerset_idx) = temp_polytope.intersect( Polyhedron('Ae',temp_Ae,'be',temp_be) );

	end


end
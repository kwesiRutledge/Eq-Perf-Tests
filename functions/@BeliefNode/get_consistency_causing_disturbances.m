function [disturb_polytope_arr] = get_consistency_causing_disturbances(bn,lcsas_in)
	%Description:
	%	This function takes the FullTrajectorySet and projects it onto the axes that define the
	%	disturbances. Such a projected polytope represents the disturbances which can cause the
	%	hypothesis (subL) of this node.
	%
	%Usage:
	%	disturb_polytope = bn.get_consistency_causing_disturbances()

	%% Input Processing
	if ~bn.phi_set_present
		error('This function cannot be called when the FullTrajectorySet is not defined for this BeliefNode.')
	end

	if ~isa(lcsas_in,'LCSAS')
		error(['The input system lcsas_in should be an LCSAS object. Instead it is of class ' class(lcsas_in) '.' ] );
	end

	%% Constants

	n_x = size(lcsas_in.Dyn(1).A,2);
	n_u = size(lcsas_in.Dyn(1).B,2);
	n_w = size(lcsas_in.Dyn(1).B_w,2);
	n_y = size(lcsas_in.Dyn(1).C,1);
	n_v = size(lcsas_in.Dyn(1).C_v,2);

	t = bn.t;

	%% Algorithm

	temp_fts = bn.FullTrajectorySet;

	if bn.subL.cardinality() == 1
		%Return only a single disturbance polytope.
		disturb_polytope_arr = temp_fts.projection( n_y*(t+1) + n_u*t + ...
				[1:(n_w+n_v)*t+n_x]);

	elseif bn.subL.cardinality() > 1

		ib_size = n_y*(t+1)+n_u*t + (n_w+n_v*t+n_x+n_x*(t+1));

		%Return the polytopes for each word which can lead to the measurements.
		disturb_polytope_arr = [];

		for word_idx = 1:bn.subL.cardinality()
			disturb_polytope_arr = [ disturb_polytope_arr , ...
				temp_fts.projection( ib_size*(word_idx-1) + n_y*(t+1) + n_u*t + [1:(n_w+n_v)*t+n_x]) ];
		end
end
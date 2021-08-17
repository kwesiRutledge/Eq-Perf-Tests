function [ explaining_node, explaining_node_index ] = NodeWhichExplainsExternalBehavior( BG , external_behavior_in )
	%Description:
	%
	%
	%Inputs:
	%	external_behavior_in: a n_y*(t+1)+n_u*t or n_x*(t+1)+n_u*t vector containing the external behavior received so far.
	%
	%Usage:
	%	node_index = BG.NodeWhichExplainsExternalBehavior( external_behavior_in )

	%% Constants

	FeedbackMethod = BG.FeedbackMethod;
	System = BG.lcsas;

	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	eb_dim = length(external_behavior_in);

	%% Algorithm

	switch FeedbackMethod
	case 'state'
	
		t = (eb_dim - n_x)/(n_x+n_u); %Get Time

	case 'output'

		t = (eb_dim - n_y)/(n_y+n_u); %Get time

	otherwise
		error(['Unexpected FeedbackMethod value (''' FeedbackMethod ''') received by NodeWhichExplainsExternalBehavior().'])
	end

	N_t = BG.get_all_nodes_at_time(t);

	cand_BNs = []; cand_BN_indices = [];
	for node_index = 1:length(N_t)
		temp_bn = N_t(node_index);
		
		%Perform set membership check.
		if ~isempty(temp_bn.c_set)
			consistency_set = temp_bn.c_set;
			obsvd_traj_is_contained = consistency_set.contains(external_behavior_in);
		elseif ~isempty(temp_bn.FullTrajectorySet)
			%Use Optimization to Verify Consistency Set Containment
			external_behavior_set = temp_bn.FullTrajectorySet;
			obsvd_traj_is_contained = external_behavior_set.containsExternalBehavior( external_behavior_in );
			%obsvd_traj_is_contained = contr.external_behavior_explained_by( obsvd_traj , external_behavior_set , 'verbosity' , verbosity );
		end

		if obsvd_traj_is_contained
			cand_BNs = [cand_BNs,temp_bn];
			cand_BN_indices = [cand_BN_indices,node_index];
		end
	end

	% Select the Node with Highest Cardinality
	explaining_node = cand_BNs.GetNodeWithHighestCardinality();
	explaining_node_index = BG.find_node_idx(explaining_node);

end
function gain_idx = prefix_detection( varargin )
	%Description:
	%	Performs detection of the prefix of the discrete mode signal when given:
	%	- Observed Mode Sequence, or
	%	- Measured Output only.
	%	
	%Usage:
	%	gain_idx = contr.prefix_detection( observed_w )
	%	gain_idx = contr.prefix_detection( y_mat )
	%
	%Inputs:
	%	- y_mat:A n_y x ow_len matrix that represents the inputs observed from time
	%			t0 in the first column to time t0+ow_len-1 in the final column.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	contr = varargin{1};

	%This controller uses a prefix of beliefs and thus its mode is not observed.
	mode_type = 'unobserved';
	y_mat = varargin{2};

	argidx = 3;
	while argidx <= nargin
		switch varargin{argidx}
			case 'verbosity'
				verbosity = varargin{argidx+1};
				argidx = argidx + 2;
			otherwise
				error('Unrecognized input to prefix_detection.')
		end
	end

	if ~exist('verbosity')
		verbosity = 0;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	BG = contr.BG;
	ConsistencySetAvailable = BG.UsedProjection; %A Consistency Set exists only if the Belief Graph was made using projection.

	u_hist = contr.u_hist;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	y_vec = reshape(y_mat,prod(size(y_mat)),1);
	u_vec = reshape(contr.u_hist,prod(size(contr.u_hist)),1);

	%Get the current time and find the nodes that correspond to this time.
	t = size(y_mat,2) - 1;

	N_t = contr.BG.get_all_nodes_at_time(t);

	%Perform Set inclusion on each of the nodes at time t and report the results.
	cand_BNs = [];
	for node_idx = 1:length(N_t)
		temp_bn = N_t(node_idx);
		
		%Perform set membership check.
		obsvd_traj = [y_vec;u_vec];

		if ConsistencySetAvailable
			consistency_set = temp_bn.c_set;
			obsvd_traj_is_contained = consistency_set.contains(obsvd_traj);
		else
			%Use Optimization to Verify Consistency Set Containment
			external_behavior_set = temp_bn.FullTrajectorySet;
			obsvd_traj_is_contained = external_behavior_set.containsExternalBehavior( obsvd_traj );
			%obsvd_traj_is_contained = contr.external_behavior_explained_by( obsvd_traj , external_behavior_set , 'verbosity' , verbosity );
		end

		if obsvd_traj_is_contained
			
			if verbosity > 0
				disp(['Belief Node ' num2str(BG.find_node_idx(temp_bn)) ' contains the trajectory at time t = ' num2str(t) '.'])
			end

			cand_BNs = [cand_BNs,temp_bn];
		end
	end

	cardinality_of_cand_BNs = [];
	for bn_idx = 1:length(cand_BNs)
		cardinality_of_cand_BNs(bn_idx) = cand_BNs(bn_idx).subL.cardinality();
	end

	%Search through the BNs that are consistent with this and find the one with the largest cardinality of
	%the sublanguage.
	[~,detected_bn_idx] = max(cardinality_of_cand_BNs);
	detected_bn = cand_BNs(detected_bn_idx);

	contr.b_hist = [contr.b_hist,BG.find_node_idx(detected_bn)]; 

	[gain_idx,~] = contr.BG.BeliefLanguage.find_a_word_with_pref(contr.b_hist);
	%error('Currently this type of system is not supported.');

end
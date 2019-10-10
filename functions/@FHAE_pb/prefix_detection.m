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

	%Identify if this is a controller that uses prefix of beliefs or prefixes of mode
	%signals in its feedback.

	if (~isempty(contr.L)) && (isempty(contr.BG))
		mode_type = 'observed';
	elseif isempty(contr.L) && (~isempty(contr.BG))
		mode_type = 'unobserved';
	else
		error('L or BG may be defined in the controller but not both. Only one of them may be.')
	end

	%Compute other inputs based on this.
	switch mode_type
		case 'observed'
			observed_w = varargin{2};
		case 'unobserved'
			y_mat = varargin{2};
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	y_vec = reshape(y_mat,prod(size(y_mat)),1);
	u_vec = reshape(contr.u_hist,prod(size(contr.u_hist)),1);

	BG = contr.BG;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	switch mode_type
		case 'observed'
			L = contr.L;
			gain_idx = L.find_a_word_with_pref(observed_w);
						
		case 'unobserved'
			u_hist = contr.u_hist;

			%Get the current time and find the nodes that correspond to this time.
			t = size(y_mat,2) - 1;

			N_t = contr.BG.get_all_nodes_at_time(t);

			%Perform Set inclusion on each of the nodes at time t and report the results.
			cand_BNs = [];
			cand_BNs_subL_card = [];
			for node_idx = 1:length(N_t)
				temp_bn = N_t(node_idx);
				%Perform set membership check.
				obsvd_traj = [y_vec;u_vec]
				temp_bn.c_set.contains(obsvd_traj)
				if temp_bn.c_set.contains(obsvd_traj)
					disp(['Belief Node ' num2str(BG.find_node_idx(temp_bn)) ' contains the trajectory at time t = ' num2str(t) '.'])
					cand_BNs = [cand_BNs,temp_bn];
					cand_BNs_subL_card = length(temp_bn.subL.words);
				end
			end

			%Search through the BNs that are consistent with this and find the one with the largest cardinality of
			%the sublanguage.
			[~,detected_bn_idx] = max(cand_BNs_subL_card);
			detected_bn = cand_BNs(detected_bn_idx);

			contr.b_hist = [contr.b_hist,contr.BG.find_node_idx(detected_bn)]; 

			[gain_idx,~] = contr.BG.BeliefLanguage.find_a_word_with_pref(contr.b_hist);
			%error('Currently this type of system is not supported.');
	end
	


end
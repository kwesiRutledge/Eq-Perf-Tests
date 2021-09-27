function [ tf , norm_matrix_diff , vector_diff ] = check_reachability_condition( cbc , ReachabilityDualVariables , X_Target , use_ol_matrices )
	%Description:
	%	


	%% Input Processing %%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%
	
	[ ibs_array , H_array , h_array , H_T_array , h_T_array ] = create_containment_target_matrices1( cbc , X_Target , use_ol_matrices );

	eps0 = 10^(-3);

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	norm_matrix_diff = {};
	vector_diff = {};

	tf = true;

	for sequence_index = 1:length(ibs_array)
		% Get Dual Variables
		Lambda_si = ReachabilityDualVariables{sequence_index};

		% Get Polytope Variables
		H_si = H_array{sequence_index}; h_si = h_array{sequence_index};
		H_T_si = H_T_array{sequence_index}; h_T_si = h_T_array{sequence_index};

		norm_matrix_diff{sequence_index} = norm(Lambda_si * H_si - H_T_si,inf);
		vector_diff{sequence_index} = Lambda_si * h_si - h_T_si;

		tf = tf && ( norm_matrix_diff{sequence_index} < eps0 ) && (max(vector_diff{end}) < 0);

	end

end

function [ ibs_array , H_array , h_array , H_T_array , h_T_array ] = create_containment_target_matrices1( cbc , X_Target , use_ol_matrices )
	%Description:
	%	

	%% Constants %%
	System = cbc.System;
	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	KnowledgeSequences = cbc.KnowledgeSequences;
	[ T , num_sequences ] = size(KnowledgeSequences);

	x0 = System.X0.V';

	%% Create IBS Arrays %%

	ibs_array = [];
	H_array = {}; h_array = {};
	H_T_array = {}; h_T_array = {};

	for sequence_index = 1:num_sequences
		% Get Sequence Index
		temp_knowledge_sequence = KnowledgeSequences(:,sequence_index);
		ibs = InternalBehaviorSet(	System,temp_knowledge_sequence, ...
									'OpenLoopOrClosedLoop','Closed',cbc.K_set{sequence_index},cbc.k_set{sequence_index});
		ibs_array = [ ibs_array ; ibs ];

		%Use first word in Last Language (of Knowledge Sequence) to Help Create H
		last_language = temp_knowledge_sequence(end);
		word1 = last_language.words{1};
		lastSymbol1 = word1(end);

		W_final = System.Dyn(lastSymbol1).P_w;

		% Create the Set Which Will Represent The Possible Sequence of Disturbances
		if use_ol_matrices
			ibs_ol = InternalBehaviorSet(System,temp_knowledge_sequence); %Get Open Loop Version of ibs

			A = ibs_ol.A; b = ibs_ol.b; 
			Ae = ibs_ol.Ae; be = ibs_ol.be; %Get Polytope Matrices

		else
			A = ibs.A; b = ibs.b; Ae = ibs.Ae; be = ibs.be; %Get Polytope Matrices
		end

		H = [ A ; Ae ; -Ae];
		H = [ H , zeros(size(H,1),n_w) ; zeros(size(W_final.A,1),ibs.Dim) , W_final.A ];

		h = [ b ; be ; -be ; W_final.b];

		if ~isempty(W_final.Ae)
			H = [ H ; zeros(size(W_final.Ae,1),ibs.Dim) , W_final.Ae ; zeros(size(W_final.Ae,1),ibs.Dim) , -W_final.Ae ];
			h = [ h ; W_final.be ; -W_final.be ];
		end

		%Add these to the array for output
		H_array{sequence_index} = H; h_array{sequence_index} = h;

		% Create Containment Based on Full Length Disturbances from H
		SelectX0Matrices = ibs.SelectX0();

		selectW1 = ibs.SelectW();
		selectX0 = SelectX0Matrices{1};

		selectWAll = [ 	selectW1, zeros(n_w*(T-1),n_w) ; 
					[ zeros(n_w,ibs.Dim) , eye(n_w) ] ];

		% Create Target Matrices
		[ S_w , S_u , S_C , S_x0 , S_k , S_Bw , S_Cv ] = System.get_mpc_matrices('word',word1); %Get MPC Matrices
		G = [ S_w*S_Bw + S_u*cbc.K_set{sequence_index} ];

		selectFinalState = [ zeros(n_x,n_x*T) , eye(n_x) ];

		H_T = X_Target.A * selectFinalState * G * selectWAll;
		h_T = X_Target.b - X_Target.A * selectFinalState * ( S_x0*x0 + S_w * S_k + S_u * cbc.k_set{sequence_index});

		H_T_array{sequence_index} = H_T; h_T_array{sequence_index} = h_T;

	end	

end
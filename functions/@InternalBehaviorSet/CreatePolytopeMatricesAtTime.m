function [A,b,Ae,be] = CreatePolytopeMatricesAtTime(ibs,t,L_t,ibs_settings)
	%Description:
	%	Creates the set of behaviors that are consistent with the hypothesis L_t at time t.
	%	This set does not consider hypotheses at previous times and is thus a sometimes over approximation
	%	of what can be expected.
	%
	%Usage:
	%	[A,b,Ae,be] = ibs.CreatePolytopeMatricesAtTime(t,L_t)
	%
	%Notes:
	%	It's possible that this can be simplified. We are enforcing a lot of constraints that don't need to be enforced.
	%	i.e. the set is enforcing constraints also on times t-1,t-2, etc. but it does not need to.
	%
	%	Each trajectory contained in the internal behavior set starts from t=0. So, the internal behavior_set_at_t=1
	%	will contain 
	%
	%			[ x0 ]
	%		x = [ x1 ] , u = [ u0 ] , w = [ w0 ]
	%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%
	System = ibs.System;
	KnowledgeSequence = ibs.KnowledgeSequence;

	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% X0 Stuff
	P_x0_L = 1;
	for word_index = 1:KnowledgeSequence(1).cardinality()
		%Initial Condition Set for Each Word in L
		P_x0_L = P_x0_L * System.X0;
	end

	% Create large disturbance set from Cartesian products
	P_wT = 1;
	for word_idx = 1:L_t.cardinality()
		for symb_idx = 1:t
			P_wT = P_wT * System.Dyn( L_t.words{word_idx}(symb_idx) ).P_w;
		end
		%P_x0_L = P_x0_L * System.X0;
    end
    
    % Created Disturbance Sets
    P_uT = 1;
    for t_idx = 1:t
        P_uT = P_uT * System.U;
    end

    % Create mpc matrices for each word in the language L
	Hc = {}; Sc = {}; Jc = {}; fc = {}; Cc = {}; Bwc = {}; Cvc = {};
	for word_ind = 1:length(L_t.words)
		[Hc{word_ind},Sc{word_ind},Cc{word_ind},Jc{word_ind},fc{word_ind},Bwc{word_ind},Cvc{word_ind}] = System.get_mpc_matrices('word',L_t.words{word_ind}(1:t));
	end

	H_block = []; S_block = []; J_block = []; f_block = [];
	I_blockx = []; I_blockx2 = [];  I_blocky = [];
	C_block = []; Cv_block = [];
	for word_ind = 1:length(L_t.words)
		H_block(end+[1:size(Hc{word_ind},1)],end+[1:size(Bwc{word_ind},2)]) = Hc{word_ind}*Bwc{word_ind};
		S_block(end+[1:size(Sc{word_ind},1)],[1:size(Sc{word_ind},2)]) = Sc{word_ind};
		%J_block(end+[1:size(Jc{word_ind},1)],[1:size(Jc{word_ind},2)]) = Jc{word_ind};
		f_block(end+[1:size(Hc{word_ind}*fc{word_ind},1)],1) = Hc{word_ind}*fc{word_ind};

		I_blockx(end+[1:n_x*(t+1)],[1:n_x*(t+1)]) = eye(n_x*(t+1));
		I_blocky(end+[1:n_y*(t+1)],[1:n_y*(t+1)]) = eye(n_y*(t+1));

	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Constructing the Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	switch ibs_settings.fb_type
	case 'state'
	    
	    %Create disturbance Polyhedron
		P_eta = P_uT * P_wT * System.X0;

		% Create the
		for word_ind = 1:L_t.cardinality()
			J_block(end+[1:size(Jc{word_ind},1)],[1:size(Jc{word_ind},2)]) = Jc{word_ind};
		end

		%Create the set of feasible (x,u,w,x0) tuples
		A = [zeros(size(P_eta.A,1),n_x*(t+1)),P_eta.A];
		b = P_eta.b;

		Ae = [-I_blockx, S_block, H_block, J_block];
		be = -f_block;

		% ib_at_t = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_x*(t+1)),P_eta.A],'b',P_eta.b, ...
		% 						'Ae',[-I_blockx, S_block, H_block, J_block],'be',-f_block );

		% Add Equality Constraints Between Words in Internal Behavior Sets
		% for word_ind = 2:L_t.cardinality()
		% 	Ltc = L_t.cardinality();
		%
		% 	%matching_x_block = [ eye(n_x*(t+1)) , zeros( n_x*(t+1) , (-2+word_ind)*n_x*(t+1) ) , -eye(n_x*(t+1)) , zeros(n_x*(t+1) , (Ltc-word_ind)*n_x*(t+1) ) ];
		%   matching_x0_block = [ eye(n_x) , zeros(n_x,(-2+word_ind)*n_x ) , -eye(n_x) , zeros(n_x , (Ltc-word_ind)*n_x ) ];
		%
		% 	Ae = [	Ae ;
		% 			zeros(n_x,n_x*(t+1)) , zeros(n_x,size(S_block,2)) , zeros(n_x,size(H_block,2)) , matching_x0_block  ];

		% 	be = [ be ; zeros(n_x,1) ];
		% end

	case 'output'

		%error('This part of CreatePolytopeMatricesAtTime() has not been tested.')

		%Also introduce the measurement disturbance into the equation
		P_vT = 1; 
		for word_idx = 1:L_t.cardinality()
			for symb_idx = 1:t+1
				P_vT = P_vT * System.Dyn( L_t.words{word_idx}(symb_idx) ).P_v;
			end
    	end

    	% Introduce Sets
    	for word_ind = 1:L_t.cardinality()
    		J_block(end+[1:size(Jc{word_ind},1)],end+[1:size(Jc{word_ind},2)]) = Jc{word_ind};

    		C_block(end+[1:n_y*(t+1)],end+[1:n_x*(t+1)]) = [Cc{word_ind} ; zeros(n_y,n_x*t), System.Dyn( L_t.words{word_ind}(t+1) ).C ];
			Cv_block(end+[1:n_y*(t+1)],end+[1:n_v*(t+1)]) = [Cvc{word_ind},zeros(size(Cvc{word_ind},1),n_v);zeros(n_y,size(Cvc{word_ind},2)), System.Dyn( L_t.words{word_ind}(t+1) ).C_v ];
			I_blockx2(end+[1:n_x*(t+1)],end+[1:n_x*(t+1)]) = eye(n_x*(t+1));
		end

    	P_eta = P_uT * P_wT * P_vT * P_x0_L;

    	%Create the set of feasible (x,u,w,x0) tuples
    	A = [zeros(size(P_eta.A,1),n_y*(t+1)),P_eta.A,zeros(size(P_eta.A,1),L_t.cardinality()*n_x*(t+1))];
    	b = P_eta.b;

    	Ae = [	zeros(size(S_block,1),size(I_blocky,2)),S_block, H_block, zeros(size(S_block,1),size(Cv_block,2)), J_block, -I_blockx2; ...
    			I_blocky, zeros(size(I_blocky,1),size(S_block,2)+size(H_block,2)), -Cv_block , zeros(size(I_blocky,1),size(J_block,2)) , -C_block ];
    	be = [-f_block;zeros(size(I_blocky,1),1)];

    	% ib_at_t = Polyhedron(	'A',[zeros(size(P_eta.A,1),n_y*(t+1)),P_eta.A,zeros(size(P_eta.A,1),length(L.words)*n_x*(t+1))],'b',P_eta.b, ...
    	% 						'Ae',[zeros(size(S_block,1),size(I_blocky,2)),S_block, H_block, zeros(size(S_block,1),size(Cv_block,2)), J_block, -I_blockx2; ...
    	% 							  I_blocky, zeros(size(I_blocky,1),size(S_block,2)+size(H_block,2)), -Cv_block , zeros(size(I_blocky,1),size(J_block,2)) , -C_block ], ...
    	% 						'be', [-f_block;zeros(size(I_blocky,1),1)] );
    otherwise
    	error(['Unexpected fb_type :' fb_type ])

    end

end
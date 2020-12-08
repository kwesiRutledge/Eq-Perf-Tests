function [results] = observer_comparison84( varargin )
	%observer_comparison84.m
	%Description:
	%	Attempts to use PenBMI to solve the bilinear optimization
    %   problem that is created by when using our finite horizon
	%   control methods.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison84';

	disp(['Beginning ' test_name '.' ])
	disp(' ')

	disp('Defining Constants.')

    %% Problem 1: Simple One-Dimensional Problem

    % Create Some Constants

    T = 4;

    D = Polyhedron('lb',-0.25,'ub',0.25);
    D_T = 1;
    for t=0:T-1
        D_T = D_T * D;
    end

    X0 = Polyhedron('lb',-0.3,'ub',0.3);

    P_eta = D_T * X0;
    n_P_eta = size(P_eta.A,1);

    %Create Target Set
    XT = Polyhedron('lb',1,'ub',2);
    n_XT = size(XT.A,1);

    %Create R_T
    RT = [zeros(1,4),1];

    %Create K matrix
    %syms K00 K10 K11 K20 K21 K22 K30 K31 K32 K33
    K00 = sdpvar(1,1,'full');
    K10 = sdpvar(1,1,'full');
    K11 = sdpvar(1,1,'full');
    K20 = sdpvar(1,1,'full');
    K21 = sdpvar(1,1,'full');
    K22 = sdpvar(1,1,'full');
    K30 = sdpvar(1,1,'full');
    K31 = sdpvar(1,1,'full');
    K32 = sdpvar(1,1,'full');
    K33 = sdpvar(1,1,'full');
    
    K = [ K00 , zeros(1,4);
			K10, K11, zeros(1,3);
			K20, K21, K22, zeros(1,2);
			K30, K31, K32, K33, 0]
    
    k = sdpvar(4,1);

    constraints = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Constraint Matrices for word 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Create S_w1, J_1 and S_u1 matrices
    S_w1 = [	zeros(1,4);
            	1,zeros(1,3);
                ones(1,2), zeros(1,2);
                ones(1,3), zeros(1,1);
                ones(1,4)];

    S_u1 = S_w1;
    J_1 = ones(5,1);

    %Create G_xw_1 and G_xx0_1
    % You can create them symbolically but I don't yet understand how to
    % translate symbolic representation to YALMIP.
    
%     Q1 = (eye(5) + S_u1*K)^(-1)
    Q1 = ...
        [ ...
            [                                                                                                                                                                                                                              1,                                                                     0,                   0,    0, 0];
            [                                                                                                                                                                                                                           -K00,                                                                     1,                   0,    0, 0];
            [                                                                                                                                                                                                            K00*K11 - K10 - K00,                                                                  -K11,                   1,    0, 0];
            [                                                                                                                                                          K00*K11 - K10 - K20 - K00 + K00*K21 + K00*K22 + K10*K22 - K00*K11*K22,                                                   K11*K22 - K21 - K11,                -K22,    1, 0];
            [K00*K11 - K10 - K20 - K30 - K00 + K00*K21 + K00*K22 + K00*K31 + K00*K32 + K10*K22 + K00*K33 + K10*K32 + K10*K33 + K20*K33 - K00*K11*K22 - K00*K11*K32 - K00*K11*K33 - K00*K21*K33 - K00*K22*K33 - K10*K22*K33 + K00*K11*K22*K33, K11*K22 - K21 - K31 - K11 + K11*K32 + K11*K33 + K21*K33 - K11*K22*K33, K22*K33 - K32 - K22, -K33, 1]
        ];


    G_xw_1 = S_w1+S_u1*K*Q1*S_w1;
    G_xx0_1 = (eye(5) + S_u1*K*Q1)*J_1

%     G_xw_1 = ...
%         [ ...
%             [                                                                                                                0,                                   0,         0, 0];
%             [                                                                                                                1,                                   0,         0, 0];
%             [                                                                                                        K11 + 1/2,                                   1,         0, 0];
%             [                                                                              K11/2 + K21 + K22/2 - K11*K22 + 1/4,                           K22 + 1/2,         1, 0];
%             [K11/4 + K21/2 + K22/4 + K31 + K32/2 + K33/4 - K33*(K11/2 + K21 - K11*K22) - (K22*K33)/2 - K11*(K22/2 + K32) + 1/8, K22/2 + K32 + K33/2 - K22*K33 + 1/4, K33 + 1/2, 1]
%         ];

    G_1 = [ G_xw_1, G_xx0_1 ]; 

    x_tilde = (eye(5) + S_u1*K*Q1)*S_u1*k;

    %Create Dual Variable matrix
    Pi1 = sdpvar(n_XT,n_P_eta,'full');
    constraints = constraints + [ Pi1 >= 0 ];
    
    %Write Constraints
    constraints = constraints + [ Pi1*P_eta.A == XT.A * RT*G_1 ]
    constraints = constraints + [ Pi1*P_eta.b <= XT.b - XT.A*RT*x_tilde ]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Constraint Matrices for word 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Create S_w1, J_1 and S_u1 matrices
    S_w2 = [zeros(1,4);
                    1,zeros(1,3);
                    0.5,1, zeros(1,2);
                    0.25,0.5,1, zeros(1,1);
                    0.125,0.25,0.5,1];

    S_u2 = S_w2;
    J_2 = [1;0.5;0.25;0.125;0.0625];

    %Create G_xw_1 and G_xx0_1
    
    % Q2 = (eye(5) + S_u2*K)^(-1)
    Q2 = [ ...
            [                                                                                                                                                                                                                                                                            1,                                                                                 0,                     0,    0, 0];
            [                                                                                                                                                                                                                                                                         -K00,                                                                                 1,                     0,    0, 0];
            [                                                                                                                                                                                                                                                        K00*K11 - K10 - K00/2,                                                                              -K11,                     1,    0, 0];
            [                                                                                                                                                                                            (K00*K11)/2 - K10/2 - K20 - K00/4 + K00*K21 + (K00*K22)/2 + K10*K22 - K00*K11*K22,                                                             K11*K22 - K21 - K11/2,                  -K22,    1, 0];
            [(K00*K11)/4 - K10/4 - K20/2 - K30 - K00/8 + (K00*K21)/2 + (K00*K22)/4 + K00*K31 + (K00*K32)/2 + (K10*K22)/2 + (K00*K33)/4 + K10*K32 + (K10*K33)/2 + K20*K33 - (K00*K11*K22)/2 - K00*K11*K32 - (K00*K11*K33)/2 - K00*K21*K33 - (K00*K22*K33)/2 - K10*K22*K33 + K00*K11*K22*K33, (K11*K22)/2 - K21/2 - K31 - K11/4 + K11*K32 + (K11*K33)/2 + K21*K33 - K11*K22*K33, K22*K33 - K32 - K22/2, -K33, 1]
        ;]

    G_xw_2 = S_w2+S_u2*K*Q2*S_w2
    G_xx0_2 = (eye(5) + S_u2*K*Q2)*J_2

    G_2 = [ G_xw_2, G_xx0_2 ]; 

    x_tilde = (eye(5) + S_u2*K*Q2)*S_u2*k;

    %Create Dual Variable matrix
    Pi2 = sdpvar(n_XT,n_P_eta,'full');
    constraints = constraints + [ Pi2 >= 0 ];
    
    %Write Constraints
    constraints = constraints + [Pi2*P_eta.A == XT.A * RT*G_2];
    constraints = constraints + [Pi2*P_eta.b <= XT.b - XT.A*RT*x_tilde];
    
    % Optimize !!
   	ops = sdpsettings('verbose',1,'solver','penlab');
    opt_out1 = optimize(constraints, ...
						[] , ...
						ops)
    
    results.Experiment1.opt_out = opt_out1;
    results.Experiment1.K = value(K);
    results.Experiment1.k = value(k);
    
    return
    clear all
                    
    %% Problem 2: Bouncing Ball Type of Dynamics
    
	disp('- Defining Constants.')

	m = 1;
	l = 1;
	g = 10;
	k = 100;
	d = .1;
	h = .01;

	%% Create Dynamics Matrices
	%
	%	x-dot = { A1 x + B1 u, 		if (x,u) \in D1
	%			{ A2 x + B2 u + c2,	if (x,u) \in D2
	%

	A1 = [ 	0 , 1 ;
			g/l , 0 ];

	B1 = [ 0 ; 1/(m*l^2) ];

	A2 = [ 	0 , 1 ;
			(g/l) - (k/m), 0 ];

	B2 = B1;

	c2 = [0; (k*d)/(m*l)];

	%% Create the Domains of each of these functions
	n = 2;
	eta_x = 1.5;
	X_template = Polyhedron('lb',-eta_x*ones(1,n),'ub',eta_x*ones(1,n));
	eta_u = 4;
	Pu = Polyhedron('lb',-eta_u*ones(1,n),'ub',eta_u*ones(1,n));

	D1 = X_template.intersect( Polyhedron('A',[1,0],'b',d/l) );
	D1 = D1 * Pu;

	D2 = X_template.intersect( Polyhedron('A',-[1,0],'b',-d/l) );
	D2 = D2 * Pu;


	eta_w = 0.4;
	Pw_template = Polyhedron('lb',-eta_w*ones(1,n),'ub',eta_w*ones(1,n));
	Pw1 = Pw_template + [eta_w;0];
	Pw2 = Pw_template + [0;eta_w];

	eta_v = 0.2;
	Pv = Polyhedron('lb',-eta_v*ones(1,n),'ub',eta_v*ones(1,n));

	ad1 = Aff_Dyn(A1,B1,zeros(n,1),eye(n),Pw1,Pv);
	ad2 = Aff_Dyn(A2,B2,c2,eye(n),Pw2,Pv);

	sys = LCSAS( [ad1,ad2] , Language([1,1,1],[2,2,2]) );
	T = 3;

	eta_u = 4;
	Pu = Polyhedron('lb',-eta_u*ones(1,m),'ub',eta_u*ones(1,m));

	eta_x0 = 0.1;
	Px0 = Polyhedron('lb',-eta_x0*ones(1,n),'ub',eta_x0*ones(1,n));

	results.Parameters.sys = sys;

	%Define Target
	eta_target = 1.5;
	P_target = Polyhedron('lb',-eta_target*ones(1,n),'ub',eta_target*ones(1,n));
	P_target = [2;1] + P_target;

	results.Parameters.P_target = P_target;

	% Plotting Information
	verbosity = 1;

	flags.ConstructBeliefGraph = true;

	results.Parameters.verbosity = verbosity;
	results.Parameters.flags = flags;

	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Display Belief Graph %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%

	if flags.ConstructBeliefGraph
		%% Get the Belief Graph for this Simple 2D System
		BG = BeliefGraph(sys,Pu,Px0,'verbosity',verbosity);
		results.Experiment0.BeliefGraph = BG;

		disp('- Created BeliefGraph object.')

		figure;
		BG.plot()
	end

	%%%%%%%%%%%%%%%%%
	%% T=1 Problen %%
	%%%%%%%%%%%%%%%%%

	disp('- Beginning test 1.')
	disp('  + Creating a T=1 optimization to find the first F.')
	
	%% Create Constants
	P_target1 = (1/3)*[2;1] + 2*Polyhedron('lb',-eta_target*ones(1,n),'ub',eta_target*ones(1,n));
	results.Experiment1.P_target = P_target1;

	ConstraintGenerator = constr_gen(0);

	ops = sdpsettings('verbose',verbosity);

	word1 = sys.L.words{1};
	P_eta1 = sys.Dyn(word1(1)).P_w * sys.Dyn(word1(1)).P_v * Px0;

	A_w1 = sys.Dyn(word1(1)).A;
	B_w1 = sys.Dyn(word1(1)).B;
	C_w1 = sys.Dyn(word1(1)).C;
	f_w1 = sys.Dyn(word1(1)).f;
	Cv_w1 = sys.Dyn(word1(1)).C_v;

	Pv_w1 = sys.Dyn(word1(1)).P_v;

	m_w1 = size(B_w1,1);

	word2 = sys.L.words{2};
	P_eta2 = sys.Dyn(word2(1)).P_w * sys.Dyn(word2(1)).P_v * Px0;

	A_w2 = sys.Dyn(word1(2)).A;
	B_w2 = sys.Dyn(word1(2)).B;
	C_w2 = sys.Dyn(word1(2)).C;
	f_w2 = sys.Dyn(word1(2)).f;
	Cv_w2 = sys.Dyn(word1(2)).C_v;	

	%% Create Optimization Variables
	F = sdpvar(m,n,'full');
	f = sdpvar(m,1,'full');

	%% Create Constraints
    
    dual_var_idx = 0;
    
	for blang_idx = 1:length(BG.BeliefLanguage.words)
		bword_ut = BG.BeliefLanguage.words{blang_idx};

		final_lang = BG.N(bword_ut(end)).subL;
		T_i = final_lang.find_longest_length();

		%Create Q,r Gains for Each Belief Word
		Q{blang_idx} = sdpvar(n_u*T_i,n_y*T_i,'full');
		r{blang_idx} = sdpvar(n_u*T_i,1,'full');

		%Create reachability constraints for each word that could be occurring
		%(all words compatible with this belief sequence)
		for word_idx = 1:length(final_lang.words)
			possible_word = final_lang.words{word_idx};

            dual_var_idx = dual_var_idx + 1;
            [ Pi1{dual_var_idx} , Piu{dual_var_idx} , temp_constrs ] = cg.get_robust_reachability_constraints(	in_lcsas, possible_word, ...
                                                                                    P_x0, Q{blang_idx},r{blang_idx}, ...
                                                                                    'P_des', P_target ,'P_u' , P_u);
            constraints = constraints + temp_constrs;

            [ Pi2{dual_var_idx} , Piu2{dual_var_idx} , temp_constrs ] = cg.get_robust_invariance_constraints(	in_lcsas, possible_word , ...
                                                                                    P_x0,Q{blang_idx},r{blang_idx}, ...
                                                                                    'eta_des', M2 , 'P_u' , P_u );

            constraints = constraints + temp_constrs;

		end

		disp(['Created the constraints for ' num2str(blang_idx) ' belief words.' ])

    end
    
    %Add Block Lower Diagonal Constraints
	l_diag_constrs = cg.get_causal_constr_on_extd_gains(in_lcsas,Q);
	constraints = constraints + l_diag_constrs;

	%Insert prefix constraint
	pref_constrs = cg.create_prefix_constr_on_gains( in_lcsas , Q , r );
	constraints = constraints + pref_constrs;
    
    results.constraints = constraints;

	%% Optimize Using PenBMI
    objective = [] ; %This is a feasibility problem; no objective.
    
	opt_out = optimize(	inclusion_constr1 + inclusion_constr2 + [ Px0.A*x <= Px0.b , Pv_w1.A * v <= Pv_w1.b ], ...
						objective , ...
						ops);

	results.Experiment1.OptimizationData1 = opt_out;
	results.Experiment1.Objective1 = value(objective); %This tells us the maximum norm of the input right?
	results.Experiment1.F1 = value(F);
	results.Experiment1.f1 = value(f);

	opt_out = optimize(	inclusion_constr1 + inclusion_constr2 + input_constr1 + input_constr2, ...
						[] , ...
						ops);

	results.Experiment1.OptimizationData2 = opt_out;
	results.Experiment1.F2 = value(F);
	results.Experiment1.f2 = value(f);

	
        
end
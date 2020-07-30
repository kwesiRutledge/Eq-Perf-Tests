function [results] = observer_comparison82( varargin )
	%observer_comparison82.m
	%Description:
	%	Playing around with the pendulum + elastic wall example from Tobia Marcucci's paper.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison82';

	disp(['Beginning ' test_name '.' ])
	disp('Testing Pendulum with Elastic Walls systems.')
	disp(' ')

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

	flags.ConstructBeliefGraph = false;

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
	[Lambda1,inclusion_constr1] = ConstraintGenerator.get_H_polyt_inclusion_constr( ...
		P_eta1.A , P_eta1.b , ...
		P_target1.A * [ eye(n) , B_w1*F , A_w1 + B_w1*F*C_w1 ] , ...
		P_target1.b - P_target1.A * B_w1 * f );

	[Lambda2,inclusion_constr2] = ConstraintGenerator.get_H_polyt_inclusion_constr( ...
		P_eta2.A , P_eta2.b , ...
		P_target1.A * [ eye(n) , B_w2*F , A_w2 + B_w2*F*C_w2 ] , ...
		P_target1.b - P_target1.A * B_w2 * f );

	[LambdaU1,input_constr1] = ConstraintGenerator.get_H_polyt_inclusion_constr( ...
		P_eta1.A , P_eta1.b , ...
		Pu.A*[ zeros(m,n) , F*Cv_w1 , F*C_w1 ] , ...
		Pu.b - Pu.A*f);

	[LambdaU2,input_constr2] = ConstraintGenerator.get_H_polyt_inclusion_constr( ...
		P_eta2.A , P_eta2.b , ...
		Pu.A*[ zeros(m,n) , F*Cv_w2 , F*C_w2 ] , ...
		Pu.b - Pu.A*f);

	% No need to constrain the F matrix to lower block-diagonal

	%% Optimize for Word 1 and/or Word 2 simultaneously
	x = sdpvar(n,1,'full');
	v = sdpvar(size(Cv_w1,2),1,'full');

	objective = -norm( F*C_w1*x + F*Cv_w1*v + f , inf);

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

	%%%%%%%%%%%%%%%%%
	%% T=2 Problem %%
	%%%%%%%%%%%%%%%%%

	disp('- Beginning test 2.')
	disp('  + Creating a T=2 optimization to find the first F.')

	%% Create Constants
	P_target2 = (2/3)*[2;1] + 2*Polyhedron('lb',-eta_target*ones(1,n),'ub',eta_target*ones(1,n));
	results.Experiment2.P_target = P_target2;

	ConstraintGenerator = constr_gen(0);

	ops = sdpsettings('verbose',verbosity);

	for word_idx = 1:sys.L.cardinality()
		word_i = sys.L.words{word_idx};
		word_i = word_i([1:T]);
		[ H{word_idx}, S{word_idx}, C_bar{word_idx}, J{word_idx}, f_bar{word_idx} ] = sys.get_mpc_matrices('word',word_i);

		P_wT{word_idx} = 1;
		P_vT{word_idx} = 1;
		for symbol_idx = 1:length(word_i)
			symb_i = word_i(symbol_idx);

			P_wT{word_idx} = P_wT{word_idx} * sys.Dyn(symb_i).P_w;
			P_vT{word_idx} = P_vT{word_idx} * sys.Dyn(symb_i).P_v;
		end
		P_eta{word_idx} = P_wT{word_idx} * P_vT{word_idx} * Px0;
	end

	n_x = size(sys.Dyn(1).A,2);
	n_u = size(sys.Dyn(1).B,2);

	T = 2;

	select_m = @(t,T_r) [zeros(n_x,t*n_x), eye(n_x), zeros(n_x,(T_r-t)*n_x) ];

	%% Create Optimization Variables
	for word_idx = 1:sys.L.cardinality()
		Q{word_idx} = sdpvar(n_u*T,n_x*(T+1),'full');
		r{word_idx} = sdpvar(n_u*T,1,'full');
	end

	%% Constraints

	%Create Target Constraints
	target_constraints = [];
	for word_idx = 1:sys.L.cardinality()
		word_i = sys.L.words{word_idx};
		[ Pi{word_idx} , Piu{word_idx} , constraints_i ] = ConstraintGenerator.get_robust_reachability_constraints( ...
			sys,word_i,Px0,Q{word_idx},r{word_idx},'P_des',P_target2, 'P_u' , Pu )

	end

	%Constrain Q's upper left block element to be defined by F.

end
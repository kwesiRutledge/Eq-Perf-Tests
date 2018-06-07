function [results] = observer_comparison29(varargin)
%observer_comparison29.m	
%	Description:
%		The objective of this test is to test new functionality that will be used
%		in the synthesis of estimators/filters when the constraint is an automaton instead of a language.
%		1. Synthesis when the provided language has different length words.
%		2. Synthesis according to the 'Free' recovery problem.

	%%%%%%%%%%%%%%%%%%%
	%% Manage Inputs %%
	%%%%%%%%%%%%%%%%%%%

	if nargin == 0
		verbosity	= 0;
	elseif nargin == 1
		verbosity	= varargin{1};
	else
		error('Unacceptable number of arguments.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	L = [];
	L = [	1,0,1,1,1,1;
			1,1,0,1,1,1;
			1,1,1,0,1,1;
			1,1,1,1,0,1];

	L2 = {};
	for ind = 1:size(L,1)
		L2{ind} = L(ind,:);
	end

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	%Create Aff_Dyn object with the data from acc_e
	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));

	acc_ad = Aff_Dyn(acc_e.A,acc_e.B,zeros(size(acc_e.A,1),1), acc_e.C, acc_e.d , acc_e.m, acc_e.E , eye(size(acc.C,1)) );

	n = size(acc_ad.A,1);
	m = size(acc_ad.B,2);
	p = size(acc_ad.C,1);
	wd = size(acc_ad.B_w,2);
	vd = size(acc_ad.C_v,2);

	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	M1 = 1;
	T = size(L,2);
	ops = sdpsettings('verbose',verbosity);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Language Consisting of Different Length Words Tests %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('1. Replace old synthesis method (where array is a constraint),')
	disp('     with one where a cell array is the constraint.')
	disp('     Problem Details:')
	disp('		- Minimize M2')

	%Perform Optimization USING MATRIX constraint

	% Optimization Variables
	% ++++++++++++++++++++++
	w		= sdpvar(wd*T,1,'full');
	ad.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');
	alpha_l 	= sdpvar(T+1,1,'full');

	for pattern_ind = 1 : size(L,1)
		v{pattern_ind} = sdpvar(vd*T,1,'full');

		% Feedback Variables
		Q{pattern_ind} = sdpvar(m*T,p*T,'full');
		r{pattern_ind} = sdpvar(m*T,1,'full');

		% Dual Variables
		Pi_1{pattern_ind} = sdpvar(2*n*T,2*(wd+vd)*T+2*n,'full');
		Pi_2{pattern_ind} = sdpvar(2*n,2*(wd+vd)*T+2*n,'full');
	end

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];

	for pattern_ind = 1 : size(L,1)
		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(acc_ad,T,'missing',find(L(pattern_ind,:) == 0)-1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
		end

		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ acc_ad.eta_w * ones(2*wd*T,1) ; acc_ad.eta_v * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*kron(eye(T),acc_ad.B)*r{pattern_ind} ];
		noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ acc_ad.eta_w * ones(2*wd*T,1) ; acc_ad.eta_v * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*kron(eye(T),acc_ad.B)*r{pattern_ind} ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T
			pre_xi = [ pre_xi ; acc_ad.A^i];
		end

		G = [ 	(eye(n*(T+1))+H0*Q{pattern_ind}*Cm0)*S0*B_w_big ...
				H0*Q{pattern_ind}*C_v_big ...
				(eye(n*(T+1))+H0*Q{pattern_ind}*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,vd*T+n) ;
									zeros(2*vd*T,wd*T) [ eye(vd*T) ; -eye(vd*T) ] zeros(2*vd*T,n) ;
									zeros(2*n,(vd+wd)*T) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = dual_equal_constrs + [Pi_1{pattern_ind} * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2{pattern_ind} * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

		%Lower Diagonal Constraint
		for bl_row_num = 1 : T-1
			l_diag_constr = l_diag_constr + [ Q{pattern_ind}(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
																[bl_row_num*p+1:end] ) == 0 ];
		end

		%Awd joint constraints for all 
		for patt_i = pattern_ind+1:size(L,1)
			%Match
			p1 = L(pattern_ind,:);
			p2 = L(patt_i,:);
			% Awd xor
			p_overlap = bitxor(p1,p2);
			%Bit at which things end up being different
			ind_identical = find(p_overlap,1) - 1;
			%Awd constraints
			shared_Q_constrs = shared_Q_constrs +  [Q{pattern_ind}( [1:ind_identical*m] , [1:ind_identical*p] ) == Q{patt_i}( [1:ind_identical*m] , [1:ind_identical*p] )];
			shared_r_constrs = shared_r_constrs +  [r{pattern_ind}( [1:ind_identical*m] ) == r{patt_i}([1:ind_identical*m]) ];
		end

	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim0 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr+shared_Q_constrs+shared_r_constrs, ...
			alpha_2, ...
			ops);

	opt_out = optim0;
	if opt_out.problem ~= 0
		contr = [];
	else
		% Save Feedback Matrices
		% ++++++++++++++++++++++
		Q_set = {}; r_set = {};
		F_set = {}; u0_set = {};
		for pattern_ind = 1 : size(L,1)
			%Get Parameters
			[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(acc_ad,T,'missing',find(L(pattern_ind,:) == 0)-1);

			Q_set{pattern_ind} = value(Q{pattern_ind});
			r_set{pattern_ind} = value(r{pattern_ind});
			F_set{pattern_ind} = value( (inv(value(eye(size(S0,2)) + kron(eye(T),acc_ad.B)*Q{pattern_ind}*Cm0*S0)) ) *kron(eye(T),acc_ad.B)* Q{pattern_ind});
			u0_set{pattern_ind} = value( inv(value(eye(size(S0,2)) + kron(eye(T),acc_ad.B)*Q{pattern_ind}*Cm0*S0)) *kron(eye(T),acc_ad.B)* r{pattern_ind} );
		
			%Fix up F and u0 to avoid NaN
			F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
			% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;
		end
	end

	opt_out.M2 = value(alpha_2);
	contr = FHAE_pb(L,F_set,u0_set);

	% Save Results of the Standard Matrix Tests
	% +++++++++++++++++++++++++++++++++++++++++

	results.exp1.matrix_opt_out = opt_out;
	results.exp1.matrix_contr = contr;
	results.exp1.L = L;
	results.exp1.L2 = L2;
	results.exp1.M1 = M1;

	%Perform Optimization USING CELL constraint

	% Clear old variables
	% +++++++++++++++++++

	clear w Pi_1 Pi_2 alpha_2 alpha_l Q r L T

	% Create a new one
	% ++++++++++++++++

	T_max = -1;
	for ind = 1:length(L2)
		if length(L2{ind}) > T_max
			T_max = length(L2{ind});
		end
	end

	% Optimization Variables
	% ++++++++++++++++++++++
	w			= sdpvar(wd*T_max,1,'full');
	acc_ad.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');
	alpha_l 	= sdpvar(T_max+1,1,'full');

	for pattern_ind = 1 : length(L2)
		T_i = length(L2{pattern_ind});
		v{pattern_ind} = sdpvar(vd*T_i,1,'full');

		% Feedback Variables
		Q{pattern_ind} = sdpvar(m*T_i,p*T_i,'full');
		r{pattern_ind} = sdpvar(m*T_i,1,'full');

		% Dual Variables
		Pi_1{pattern_ind} = sdpvar(2*n*T_i,2*(wd+vd)*T_i+2*n,'full');
		Pi_2{pattern_ind} = sdpvar(2*n,2*(wd+vd)*T_i+2*n,'full');
	end

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];

	for pattern_ind = 1 : length(L2)
		% Creating Constraints
		% ++++++++++++++++++++
		T_i = length(L2{pattern_ind});

		[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(acc_ad,T_i,'missing',find(L2{pattern_ind} == 0)-1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T_i
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
		end

		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ acc_ad.eta_w * ones(2*wd*T_i,1) ; acc_ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T_i,1) - [eye(n*T_i);-eye(n*T_i)]*sel_influenced_states*S0*kron(eye(T_i),acc_ad.B)*r{pattern_ind} ];
		noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ acc_ad.eta_w * ones(2*wd*T_i,1) ; acc_ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T_i,T_i)*S0*kron(eye(T_i),acc_ad.B)*r{pattern_ind} ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T_i
			pre_xi = [ pre_xi ; acc_ad.A^i];
		end

		G = [ 	(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*S0*B_w_big ...
				H0*Q{pattern_ind}*C_v_big ...
				(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T_i) ; -eye(wd*T_i) ] zeros(2*wd*T_i,vd*T_i+n) ;
									zeros(2*vd*T_i,wd*T_i) [ eye(vd*T_i) ; -eye(vd*T_i) ] zeros(2*vd*T_i,n) ;
									zeros(2*n,(vd+wd)*T_i) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = dual_equal_constrs + [Pi_1{pattern_ind} * bounded_disturb_matrix == [eye(n*T_i); -eye(n*T_i)]*sel_influenced_states*G ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2{pattern_ind} * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T_i,T_i)*G];

		%Lower Diagonal Constraint
		for bl_row_num = 1 : T_i-1
			l_diag_constr = l_diag_constr + [ Q{pattern_ind}(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
																[bl_row_num*p+1:end] ) == 0 ];
		end

		%Add joint constraints for all 
		for patt_i = pattern_ind+1:length(L2)
			%Match
			p1 = L2{pattern_ind};
			p2 = L2{patt_i};
			%Truncate one if necessary.
			if length(p1) < length(p2)
				p2 = p2(1:length(p1));
			elseif length(p1) > length(p2)
				p1 = p1(1:length(p2));
			end
			% Add xor
			p_overlap = bitxor(p1,p2);
			%Bit at which things end up being different
			ind_identical = find(p_overlap,1) - 1;
			%Awd constraints
			shared_Q_constrs = shared_Q_constrs +  [Q{pattern_ind}( [1:ind_identical*m] , [1:ind_identical*p] ) == Q{patt_i}( [1:ind_identical*m] , [1:ind_identical*p] )];
			shared_r_constrs = shared_r_constrs +  [r{pattern_ind}( [1:ind_identical*m] ) == r{patt_i}([1:ind_identical*m]) ];
		end

	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim1 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr+shared_Q_constrs+shared_r_constrs, ...
			alpha_2, ...
			ops);

	opt_out = optim1;
	if opt_out.problem ~= 0
		contr = [];
	else
		% Save Feedback Matrices
		% ++++++++++++++++++++++
		Q_set = {}; r_set = {};
		F_set = {}; u0_set = {};
		for pattern_ind = 1 : length(L2)
			T_i = length(L2{pattern_ind});
			%Get Parameters
			[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(acc_ad,T_i,'missing',find(L2{pattern_ind} == 0)-1);

			Q_set{pattern_ind} = value(Q{pattern_ind});
			r_set{pattern_ind} = value(r{pattern_ind});
			F_set{pattern_ind} = value( (inv(value(eye(size(S0,2)) + kron(eye(T_i),acc_ad.B)*Q{pattern_ind}*Cm0*S0)) ) *kron(eye(T_i),acc_ad.B)* Q{pattern_ind});
			u0_set{pattern_ind} = value( inv(value(eye(size(S0,2)) + kron(eye(T_i),acc_ad.B)*Q{pattern_ind}*Cm0*S0)) *kron(eye(T_i),acc_ad.B)* r{pattern_ind} );
		
			%Fix up F and u0 to avoid NaN
			F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
			% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;
		end
	end

	opt_out.M2 = value(alpha_2);
	contr = FHAE_pb(L2,F_set,u0_set);

	% Save Results of the Cell Matrix Tests
	% +++++++++++++++++++++++++++++++++++++

	results.exp2.T_max = T_max;
	results.exp2.cell_opt_out = opt_out;
	results.exp2.cell_contr = contr;

	%% Compare entries
	disp('Comparing the entries of test 1 and test 2''s Prefix-Based Controllers:')
	
	matrix_contr = results.exp1.matrix_contr;
	cell_contr = results.exp2.cell_contr;

	for prefix_num = 1 : length(results.exp2.cell_contr.F_set)
		%Compare the F's associated with this pattern index
		if sum(sum( matrix_contr.F_set{prefix_num} == cell_contr.F_set{prefix_num} )) == prod(size(matrix_contr.F_set{prefix_num}))
			disp(['Prefix #' num2str(prefix_num) ' F: Match.' ])
		else
			disp(['Prefix #' num2str(prefix_num) ' F: DO NOT MATCH.'])
		end

		if sum( matrix_contr.u0_set{prefix_num} == cell_contr.u0_set{prefix_num} ) == prod(size(matrix_contr.u0_set{prefix_num}))
			disp(['Prefix #' num2str(prefix_num) ' u0: Match.' ])
		else
			disp(['Prefix #' num2str(prefix_num) ' u0: DO NOT MATCH' ])
		end

	end

	disp('Identical Prefix-Based Feedback is created.')

	%%%%%%%%%%%%
	%% Test 3 %%
	%%%%%%%%%%%%

	disp('======================================================')
	disp('3. Modify the FHAE_pb class AND the synthesis function')
	disp('   eq_rec_design_pb() to reflect this new change.')
	disp(' ')

	disp('FHAE_pb has been accurately updated. See observer_comparison28().')


	L3 = L2;
	L3{length(L3)+1} = [1,0,1,0,1];

	[ opt_data, contr ] = eq_rec_design_pb( acc_ad , 'Feasible Set' , M1 , 3 , 6 );

	results.exp3.L = L3;
	results.exp3.opt_data_feas = opt_data;
	results.exp3.contr_feas = contr;

	[ opt_data , contr ] = free_rec_design_pb( acc_ad , 'Feasible Set' , M1 , 3 , M1 , 6 );
	results.exp3.opt_data_free_rec = opt_data;
	results.exp3.contr_free_rec = contr;

	[ opt_data , contr ] = free_rec_design_pb( acc_ad , 'Min_M3' , M1 , 3 , 6 );
	results.exp3.opt_data_minM3 = opt_data;
	results.exp3.contr_minM3 = contr;

	disp('Introduced the free recovery problem synthesis function.')

	%%%%%%%%%%%%
	%% Test 4 %%
	%%%%%%%%%%%%

	disp('===================================================================================================')
	disp('4. Testing how the new free recovery problem handles strange languages with different length words.')
	disp(' ')

	[ opt_data , contr ] = free_rec_design_pb( acc_ad , 'Min_M3' , M1 , 3 , L3 );
	results.exp4.L = L3;
	results.exp4.opt_data = opt_data;
	results.exp4.contr = contr;

	%%%%%%%%%%%%
	%% Test 5 %%
	%%%%%%%%%%%%

	disp('======================================================')
	disp('5. Testing the Synthesis Based on Observable Automaton')
	disp(' ')

	X2   = [0:5]';
	X0_2 = [1];
	Y2   = [0;1];
	H2 = [ 	0,0 ;
			1,1 ;
			2,0 ;
			3,1 ;
			4,0 ;
			5,0 ];
	Delta2 = [	0,1;
			 	1,0;
			 	1,2;
			 	2,3;
			 	3,4;
			 	4,5;
			 	5,1 ];
    
	fsm3 = FSM0(X2,X0_2,Y2,H2,Delta2);

	disp('Automata created.')
	disp(' ')

	oa3 = OBSV_AUT(fsm3);

	disp('Created Observable Automaton?')
	disp(' ')

	[L,S] = oa3.find_all_2R_paths()

	results.exp5.L = L;
	results.exp5.S = S;
	results.exp5.oa = oa3;
	results.exp5.fsm = fsm3;

	%%%%%%%%%%%%
	%% Test 6 %%
	%%%%%%%%%%%%

	disp('===========================================================================================')
	disp('6. Testing Synthesis When the state of the Observable Automaton Guides Constraint Selection')
	disp(' ')

	clear G

	%Experiment Parameters
	L = {results.exp5.L{[2 3]} };
	s_0 = S{2}(1);
	M2 = 3;
	M3 = sdpvar(1,1,'full');

	obj_fcn = M3;
	obj_constrs = [M3 >= 0];

	% Optimization Variables
	% ++++++++++++++++++++++
	ad.x0 	= sdpvar(n,1,'full');

	max_T_i = -1;
	for pattern_ind = 1 : length(L)
		T_i = length(L{pattern_ind});
		% Feedback Variables
		Q{pattern_ind} = sdpvar(m*T_i,p*T_i,'full');
		r{pattern_ind} = sdpvar(m*T_i,1,'full');

		% Dual Variables
		Pi_1{pattern_ind} = sdpvar(2*n*T_i,2*(wd+vd)*T_i+2*n,'full');
		Pi_2{pattern_ind} = sdpvar(2*n,2*(wd+vd)*T_i+2*n,'full');

		%Find the maximum T_i
		if T_i > max_T_i
			max_T_i = T_i;
		end
	end
	w	= sdpvar(wd*max_T_i,1,'full');

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];

	for pattern_ind = 1 : length(L)
		T_i = length(L{pattern_ind});
		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(acc_ad,T_i,'missing',find(L{pattern_ind} == 0)-1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T_i
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
		end

		%There are a pair of constraints that enforce that bounds M1 is obeyed (or bound M2, M3). One will be placed in noise_constrs the other will be in dual equal constrs.
		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ acc_ad.eta_w * ones(2*wd*T_i,1) ; acc_ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M2 * ones(2*n*T_i,1) - [eye(n*T_i);-eye(n*T_i)]*sel_influenced_states*H0*r{pattern_ind} ];
		if S{pattern_ind} == s_0
			noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ acc_ad.eta_w * ones(2*wd*T_i,1) ; acc_ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M3 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T_i,T_i)*H0*r{pattern_ind} ];
		else
			noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ acc_ad.eta_w * ones(2*wd*T_i,1) ; acc_ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M3 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T_i,T_i)*H0*r{pattern_ind} ];
		end
			
		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T_i
			pre_xi = [ pre_xi ; acc_ad.A^i];
		end

		G{pattern_ind} = [ 	(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*S0*B_w_big ...
							H0*Q{pattern_ind}*C_v_big ...
							(eye(n*(T_i+1))+H0*Q{pattern_ind}*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T_i) ; -eye(wd*T_i) ] zeros(2*wd*T_i,vd*T_i+n) ;
									zeros(2*vd*T_i,wd*T_i) [ eye(vd*T_i) ; -eye(vd*T_i) ] zeros(2*vd*T_i,n) ;
									zeros(2*n,(vd+wd)*T_i) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = dual_equal_constrs + [Pi_1{pattern_ind} * bounded_disturb_matrix == [eye(n*T_i); -eye(n*T_i)]*sel_influenced_states*G{pattern_ind} ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2{pattern_ind} * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T_i,T_i)*G{pattern_ind}];

		%Lower Diagonal Constraint
		for bl_row_num = 1 : T_i-1
			l_diag_constr = l_diag_constr + [ Q{pattern_ind}(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
																[bl_row_num*p+1:end] ) == 0 ];
		end

		%Awd joint constraints for all 
		for patt_i = pattern_ind+1:length(L)
			%Match
			p1 = L{pattern_ind};
			p2 = L{patt_i};
			%Truncate one if necessary.
			if length(p1) < length(p2)
				p2 = p2(1:length(p1));
			elseif length(p1) > length(p2)
				p1 = p1(1:length(p2));
			end
			% Add xor
			p_overlap = bitxor(p1,p2);
			%Bit at which things end up being different
			ind_identical = find(p_overlap,1) - 1;
			%Awd constraints
			shared_Q_constrs = shared_Q_constrs +  [Q{pattern_ind}( [1:ind_identical*m] , [1:ind_identical*p] ) == Q{patt_i}( [1:ind_identical*m] , [1:ind_identical*p] )];
			shared_r_constrs = shared_r_constrs +  [r{pattern_ind}( [1:ind_identical*m] ) == r{patt_i}([1:ind_identical*m]) ];
		end

	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim0 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr+shared_Q_constrs+shared_r_constrs+obj_constrs, ...
			obj_fcn, ...
			ops);

	opt_out = optim0;
	if opt_out.problem ~= 0
		contr = [];
	else
		% Save Feedback Matrices
		% ++++++++++++++++++++++
		Q_set = {}; r_set = {};
		F_set = {}; u0_set = {};
		for pattern_ind = 1 : length(L)
			T_i = length(L{pattern_ind});
			%Get Parameters
			[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(acc_ad,T_i,'missing',find(L{pattern_ind} == 0)-1);

			Q_set{pattern_ind} = value(Q{pattern_ind});
			r_set{pattern_ind} = value(r{pattern_ind});
			F_set{pattern_ind} = value( (inv(value(eye(size(S0,2)) + Q{pattern_ind}*Cm0*H0)) ) * Q{pattern_ind});
			u0_set{pattern_ind} = value( inv(value(eye(size(S0,2)) + Q{pattern_ind}*Cm0*H0)) * r{pattern_ind} );
			
			%Fix up F and u0 to avoid NaN
			F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
			% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

			%Create Function Outputs
			opt_out.Q_set = Q_set;
			opt_out.r_set = r_set;

			contr = FHAE_pb(L,F_set,u0_set);

		end
	end

	% Save Results
	% ++++++++++++
	results.exp6.opt_out = opt_out;
	results.exp6.contr = contr;

end
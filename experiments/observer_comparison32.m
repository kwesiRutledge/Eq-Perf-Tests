function [results] = observer_comparison32(varargin)
%observer_comparison32.m	
%	Description:
%		The objective of this test is to verify whether or not 2 nonlinear change of variables will be useful in defining the 'modularity'
%		concept that I've been wrestling with. We will solve a problem (an arbitrary one) where the language is the following one:
%		
%		L = { 1,0,1,1,0,1,1 }
%
%		which will correspond to the following sequence in state space
%
%		q = { q_R , q_1, q_2 , q_3 , q_R , q_4 , q_5 , q_R }
%		
%		What I would like to test, is if we can synthesize a feedback matrix F with the following form
%		
%		F = [F_1 0  ;
%			 0   F_2 ];
%
%		Which is of a block diagonal structure instead of the block lower diagonal structure that was previously used.

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

	L = [1,0,1,1,0,1,1];
	M1 = 1;

	%Using ACC System
	load('data/system_examples/acc_p.mat');
	
	%Create Aff_Dyn object with the data from acc_e
	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));

	acc_ad = Aff_Dyn(acc_e.A,acc_e.B,zeros(size(acc_e.A,1),1), acc_e.C, acc_e.d , acc_e.m, acc_e.E , eye(size(acc.C,1)) );

	disp('#32')
	disp('- Using ACC System')

	L_str = '- L = [';
	for ind_L = 1:length(L)
		L_str = [ L_str num2str(L(ind_L)) ];
		if ind_L ~= length(L)
			L_str = [ L_str ',' ];
		end
	end
	L_str = [ L_str ']'];

	disp(L_str)
	disp(['- M1 = ' num2str(M1)]);

	results.params.L = L;
	results.params.M1 = M1;
	results.params.sys = acc_ad;

	n = size(acc_ad.A,1);
	m = size(acc_ad.B,2);
	p = size(acc_ad.C,1);
	wd = size(acc_ad.B_w,2);
	vd = size(acc_ad.C_v,2);

	%Selection vector function
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%Unit Cube in MPT3
	unit_cube = Polyhedron('A',[ eye(n) ; -eye(n) ],'b',ones(2*n,1));

	results.constants.unit_cube = unit_cube;

	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Optimization Set Up %%
	%%%%%%%%%%%%%%%%%%%%%%%%%

	ops = sdpsettings('verbose',verbosity);

	T_tot = length(L);
	T(1) = 4;
	T(2) = T_tot - T(1);

	w		= sdpvar(wd*T_tot,1,'full');
	acc_ad.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');
	alpha_l 	= sdpvar(T_tot+1,1,'full');

	for seq_ind = 1 : 2
		v{seq_ind} = sdpvar(vd*T_tot,1,'full');

		% Feedback Variables
		Q{seq_ind} = sdpvar(m*T(seq_ind),p*T(seq_ind),'full');
		r{seq_ind} = sdpvar(m*T(seq_ind),1,'full');

		% Dual Variables
		Pi_1{seq_ind} = sdpvar(2*n*T(seq_ind),2*(wd+vd)*T(seq_ind)+2*n,'full');
		Pi_2{seq_ind} = sdpvar(2*n,2*(wd+vd)*T(seq_ind)+2*n,'full');
	end

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% First Half Optimization (F_1) %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% T(1) = 4;
	seq_num = 1;

	[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(acc_ad,T(1),'missing',find(L(1,[1:T(1)]) == 0)-1);

	positive_constr = positive_constr + [ Pi_1{seq_num} >= 0, Pi_2{seq_num} >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T(seq_num)
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T(seq_num)) ];
	end

	noise_constrs = noise_constrs + [ Pi_1{seq_num} * [ acc_ad.eta_w * ones(2*wd*T(seq_num),1) ; acc_ad.eta_v * ones(2*p*T(seq_num),1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T(seq_num),1) - [eye(n*T(seq_num));-eye(n*T(seq_num))]*sel_influenced_states*S0*kron(eye(T(seq_num)),acc_ad.B)*r{seq_num} ];
	noise_constrs = noise_constrs + [ Pi_2{seq_num} * [ acc_ad.eta_w * ones(2*wd*T(seq_num),1) ; acc_ad.eta_v * ones(2*p*T(seq_num),1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T(seq_num),T(seq_num))*S0*kron(eye(T(seq_num)),acc_ad.B)*r{seq_num} ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T(seq_num)
		pre_xi = [ pre_xi ; acc_ad.A^i];
	end

	G = [ 	(eye(n*(T(seq_num)+1))+H0*Q{seq_num}*Cm0)*S0*B_w_big ...
			H0*Q{seq_num}*C_v_big ...
			(eye(n*(T(seq_num)+1))+H0*Q{seq_num}*Cm0)*pre_xi ];

	G1_saved = G;

	bounded_disturb_matrix = [ [ eye(wd*T(seq_num)) ; -eye(wd*T(seq_num)) ] zeros(2*wd*T(seq_num),vd*T(seq_num)+n) ;
								zeros(2*vd*T(seq_num),wd*T(seq_num)) [ eye(vd*T(seq_num)) ; -eye(vd*T(seq_num)) ] zeros(2*vd*T(seq_num),n) ;
								zeros(2*n,(vd+wd)*T(seq_num)) [ eye(n) ; -eye(n) ] ];

	dual_equal_constrs = dual_equal_constrs + [Pi_1{seq_num} * bounded_disturb_matrix == [eye(n*T(seq_num)); -eye(n*T(seq_num))]*sel_influenced_states*G ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2{seq_num} * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T(seq_num),T(seq_num))*G];

	%Lower Diagonal Constraint
	for bl_row_num = 1 : T(seq_num)-1
		l_diag_constr = l_diag_constr + [ Q{seq_num}(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
															[bl_row_num*p+1:end] ) == 0 ];
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% First Half Optimization (F_2) %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear S0 H0 Cm0 xi0m B_w_big C_v_big G

	seq_num = 2;

	[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(acc_ad,T(seq_num),'missing',find(L(1,[T(1):end]) == 0)-1);

	positive_constr = positive_constr + [ Pi_1{seq_num} >= 0, Pi_2{seq_num} >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T(seq_num)
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T(seq_num)) ];
	end

	noise_constrs = noise_constrs + [ Pi_1{seq_num} * [ acc_ad.eta_w * ones(2*wd*T(seq_num),1) ; acc_ad.eta_v * ones(2*p*T(seq_num),1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T(seq_num),1) - [eye(n*T(seq_num));-eye(n*T(seq_num))]*sel_influenced_states*S0*kron(eye(T(seq_num)),acc_ad.B)*r{seq_num} ];
	noise_constrs = noise_constrs + [ Pi_2{seq_num} * [ acc_ad.eta_w * ones(2*wd*T(seq_num),1) ; acc_ad.eta_v * ones(2*p*T(seq_num),1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T(seq_num),T(seq_num))*S0*kron(eye(T(seq_num)),acc_ad.B)*r{seq_num} ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T(seq_num)
		pre_xi = [ pre_xi ; acc_ad.A^i];
	end

	G = [ 	(eye(n*(T(seq_num)+1))+H0*Q{seq_num}*Cm0)*S0*B_w_big ...
			H0*Q{seq_num}*C_v_big ...
			(eye(n*(T(seq_num)+1))+H0*Q{seq_num}*Cm0)*pre_xi ];

	%Use MPT to find where this controller will start.

	xi_t0_poly = unit_cube*M1;
	w_traj_poly = kron( ones(T(1),1) ,[eye(wd);-eye(wd)]) * Polyhedron('A',[eye(wd);-eye(wd)],'b',[ones(2*wd,1)]) * acc_ad.eta_w ;
	v_traj_poly = kron( ones(T(1),1) ,[eye(vd);-eye(vd)]) * Polyhedron('A',[eye(vd);-eye(vd)],'b',[ones(2*vd,1)]) * acc_ad.eta_v ;

	xi_t1_poly = select_m(T(1),T(1)) * ( 	G1_saved(:,[1:size(acc_ad.B_w,2)*T(1)]) * w_traj_poly + ...
											G1_saved(:,size(acc_ad.B_w,2)*T(1)+[1:size(acc_ad.C_v,2)*T(1)]) * v_traj_poly + ...
											G1_saved(:,size(acc_ad.B_w,2)*T(1)+size(acc_ad.C_v,2)*T(1)+[1:n]) * xi_t0_poly )


	bounded_disturb_matrix = [ [ eye(wd*T(seq_num)) ; -eye(wd*T(seq_num)) ] zeros(2*wd*T(seq_num),vd*T(seq_num)+n) ;
								zeros(2*vd*T(seq_num),wd*T(seq_num)) [ eye(vd*T(seq_num)) ; -eye(vd*T(seq_num)) ] zeros(2*vd*T(seq_num),n) ;
								zeros(2*n,(vd+wd)*T(seq_num)) [ eye(n) ; -eye(n) ] ];

	dual_equal_constrs = dual_equal_constrs + [Pi_1{seq_num} * bounded_disturb_matrix == [eye(n*T(seq_num)); -eye(n*T(seq_num))]*sel_influenced_states*G ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2{seq_num} * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T(seq_num),T(seq_num))*G];

	%Lower Diagonal Constraint
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q{seq_num}(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
															[bl_row_num*p+1:end] ) == 0 ];
	end
end
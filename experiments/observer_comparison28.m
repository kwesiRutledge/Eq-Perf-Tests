function [ results ] = observer_comparison28( varargin )
%	observer_comparison28.m
%		Description:
%			The objective of this experiment is to design a finite horizon, affine
%			estimator (fhae) that is robust against the possibility of 1 observation
%			missing in the entire sequence of length T.
%			We will compare: Worst Case Language method v. My proposal
%
%			What makes this experiment different from 13 and 14,25,27 is that:
%				- It uses the recent optimization presented by Prof. Yong.
%				- We make explicit use of the "E" matrix in the ACC example.
%				- Everything is done through functions.
%
%			This Equalized Recovery Problem:
%				Let M1 T
%				Let ||xi(0)||<= M1. What is the minimum value for M1 such that,
%				- ||xi(t)||<= M2 \forall t \in [1,T-1]
%				- ||xi(T)||<= M1 
%				OR
%				- ||xi(t)||<= M1 \forall t
%
%		Inputs:
%			verbosity - 
%			T -			Time Horizon
%			M1 - 		


	%%%%%%%%%%%%%%%%%%%
	%% Manage Inputs %%
	%%%%%%%%%%%%%%%%%%%

	if nargin == 0
		verbosity	= 0;
	elseif nargin == 1
		verbosity	= varargin{1};
	elseif nargin == 2
		verbosity	= varargin{1};
		T 			= varargin{2};
	elseif nargin == 3
		verbosity	= varargin{1};
		T	= varargin{2};
		M1 = varargin{3};
	else
		error('Unacceptable number of arguments.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%
	T_s = 0.5;

	%Default Values for observer_comparison
	if nargin < 2
		T = 6;
		M1 = 1;
	end

	L = [];
	L = [	1,0,1,1,1,1;
			1,1,0,1,1,1;
			1,1,1,0,1,1;
			1,1,1,1,0,1];

	L_star = ones(1,6);
	for word_in_L = 1:size(L,1)
		L_star = bitand(L_star,L(word_in_L,:));
	end

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	%Create Aff_Dyn object with the data from acc_e
	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));

	acc_ad = Aff_Dyn(acc_e.A,acc_e.B,zeros(size(acc_e.A,1),1), acc_e.C, acc_e.d , acc_e.m, acc_e.E , eye(size(acc.C,1)) );
	[ oc27_opt1 , oc27_contr1 ] = eq_rec_design_t( acc_ad , 'Min_M2' , M1 , T , L );

	disp('Testing a feasibility problem that is known to be feasible.')
	[ oc27_feas_opt1 , oc27_feas_contr1 ] = eq_rec_design_t( acc_ad , 'Feasible Set' , M1 , oc27_contr1.opt_obj, T , L );

	disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	disp('Testing a feasibility problem for the prefix-based feedback law version.')
	[ opt_data2 , new_contr ] = eq_rec_design_pb( acc_ad , 'Feasible Set' , M1 , oc27_contr1.opt_obj, L );

	new_contr.apply_control([1,1,0],ones(3*size(acc.C,1),1))

	disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	disp('Simulate one round with prefix-based control.')
	pb_sim1 = new_contr.simulate_1run( acc_ad , M1 );

	figure;
	plot(pb_sim1')

	disp('Figure is plotted.')

	disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	num_runs = 100;
	disp(['Simulate ' num2str(num_runs) ' rounds with prefix-based control.'])
	[ pb_sim2 , pb_sim2_norm ] = new_contr.simulate_n_runs( acc_ad , M1 , num_runs );

	pb_sim2_mod = [];
	for ind = 1:length(pb_sim2)
		pb_sim2_mod(:,:,ind) = pb_sim2{ind};
	end

	pb_sim2_norm_mod = [];
	for ind = 1:length(pb_sim2_norm)
		pb_sim2_norm_mod(:,:,ind) = pb_sim2_norm{ind};
	end

	figure;
	plot([0:T],reshape(pb_sim2_norm_mod,T+1,num_runs));
	axis([0 T 0 oc27_contr1.opt_obj+0.2 ])

	title('State Based Observer with Feasible Set of (M_1,M_2,T)')


	%Test Minimization of M2 with prefix based observers.
	disp('+++++++++++++++++++++++++++++++++++++++++++++++++++')
	disp('Test Minimization of M2 with prefix based observer.')
	[ pb_opt3 , pb_contr3 ] = eq_rec_design_pb( acc_ad , 'Min_M2' , M1 , L );
	[ pb3_sim , pb3_sim_norm ] = pb_contr3.simulate_n_runs( acc_ad , M1 , num_runs );

	pb3_sim_norm_mod = [];
	for ind = 1:length(pb_sim2_norm)
		pb3_sim_norm_mod(:,:,ind) = pb3_sim_norm{ind};
	end

	figure;
	bar_heights = [ M1 pb_opt3.M2*ones(1,T-1) M1];

	subplot(1,2,1)
	hold on;
	bar([0:T],bar_heights,'w')
	plot([0:T],reshape(pb3_sim_norm_mod,T+1,num_runs));
	
	axis([0-0.5 T+0.5 0 oc27_contr1.opt_obj+0.2 ])
	xlabel('Time Index')
	ylabel('$||x(t)-\hat{x}(t)||_{\infty}$','Interpreter','latex')
	title('Prefix-Based Observer with Minimized M_2')
	disp('Plotted Prefix-based solution.')

	%Time Based Simulation
	[ pb_opt4 , pb_contr4 ] = eq_rec_design_pb(acc_ad , 'Min_M2' , M1 , L_star );
	[ pb4_sim , pb4_sim_norm ] = pb_contr4.simulate_n_runs( acc_ad , M1 , num_runs );

	pb4_sim_norm_mod = [];
	for ind = 1:length(pb_sim2_norm)
		pb4_sim_norm_mod(:,:,ind) = pb4_sim_norm{ind};
	end

	bar_heights = [ M1 pb_opt4.M2*ones(1,T-1) M1];

	subplot(1,2,2);
	hold on;
	bar([0:T],bar_heights,'w')
	plot([0:T],reshape(pb4_sim_norm_mod,T+1,num_runs));
	
	axis([0-0.5 T+0.5 0 oc27_contr1.opt_obj+0.2 ])
	xlabel('Time Index')
	ylabel('$||x(t)-\hat{x}(t)||_{\infty}$','Interpreter','latex')
	title('(Pseudo) Time-Based Observer with Minimized M_2')
	disp('Plotted time-based solution.')

	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc_ad;
	results.L = L;
	results.L_star = L_star;

	results.oc27_opt1 = oc27_opt1;
	results.oc27_contr1 = oc27_contr1;
	
	results.oc27_feas_opt1 = oc27_feas_opt1;
	results.oc27_feas_contr1 = oc27_feas_contr1;

	results.pb_exp1_opt = opt_data2;
	results.pb_exp1_contr = new_contr;

	results.pb_sim1 = pb_sim1;
	results.pb_sim2.x = pb_sim2_mod;
	results.pb_sim2.x_norms = pb_sim2_norm_mod;

	results.pb3.opt = pb_opt3;
	results.pb3.contr = pb_contr3;
	results.pb3.x = pb3_sim;
	results.pb3.x_norm = pb3_sim_norm_mod;

	results.pb4.opt = pb_opt4;
	results.pb4.contr = pb_contr4;
	results.pb4.x = pb4_sim;
	results.pb4.x_norm = pb4_sim_norm_mod;

end


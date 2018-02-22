function [ results ] = observer_comparison21( varargin )
%	observer_comparison21.m
%		Description:
%			The objective of this experiment is to test the new functionality in the
%			create fhae function, which should help make the data pattern customizable.
%
%			What makes this experiment different from 13,14,16 and 18 is that:
%				- It uses the recent optimization presented by Prof. Yong.
%				- We make explicit use of the "E" matrix in the ACC example.
%				- A reduced order observer is incorporated to try to achieve Equalized Performance if necessary.
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
%			verbosity - Desired Verbosity (amount of status updates/warnings that you receive)
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

	%Default Values for observer_comparison
	if nargin < 2
		T = 6;
		M1 = 1;
	end

	%Define our example system
	%Using ACC System
	load('data/system_examples/acc_p.mat');

	%Create Reduced Order Observer's Error Dynamics
	acc_roo.A = acc.A(3,3);
	acc_roo.B = eye(1);
	acc_roo.C = acc.A([1:2],3);
	acc_roo.E = [acc.E(3,1) acc.A(3,[1:2]) zeros(1,3) ];
	% acc_roo.G1 = [acc.E([1:2],1) acc.A([1:2],[1:2])];
	% acc_roo.G2 = [acc.E([1:2],1) zeros(2) ];
	acc_roo.G = [acc.E([1:2],1) acc.A([1:2],[1:2]) zeros(2,1) eye(2) ];
	acc_roo.m = acc.m;
	acc_roo.d = acc.d; 

	%Create Error System
	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));

	test_sys = acc_e;

	%Create dimensions
	n = size(test_sys.A,1);
	p = size(test_sys.C,1);

	m = size(test_sys.B,2);
	dd = size(test_sys.E,2); %Delta Dimension

	b_dim = size(acc_roo.G,2);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Synthesize Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Announce New Experiments
	disp('========================================================================================')
	disp('Experiment 1: Synthesis of default fhae function vs. with data pattern explicitly chosen')
	disp(' ')
	disp(['The standard function is believed to assume the following data pattern for n = ' num2str(n) ' is:' ])

	exp1.pattern1 = [];
	for i = 2 : T-2
		temp_row = ones(1,T);
		temp_row(i) = 0;
		exp1.pattern1 = [ exp1.pattern1 ; temp_row ];
	end

	disp(exp1.pattern1)

	exp1.pattern2 = ones(1,T);
	disp(exp1.pattern2)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Standard function's design %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[contr1_fcn,optim1_fcn] = create_fhae_w_eq_recovery('min_M2',M1,T,test_sys,'verbosity',2)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Custom Pattern Function's Design %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[contr2_fcn,optim2_fcn] = create_fhae_w_eq_recovery('min_M2',M1,T,test_sys,'verbosity',2,'pattern',exp1.pattern2)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot Results for First Experiment Test %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	num_rollouts = 1000;
	% contr1.F = contr1_fcn.F;
	% contr1.u0 = contr1.u0;



	test_sys.x0 = zeros(n,1);

	[exp1.xi_t,  exp1.xi_mag_t] 	= apply_controller_to_rollouts(test_sys,contr1_fcn,T,num_rollouts,M1);
	[exp1.xi2_t, exp1.xi2_mag_t ] = apply_controller_to_rollouts(test_sys,contr2_fcn,T,num_rollouts, M1);

	rollouts_per_plot = 1000;

	figure;
	subplot(1,2,1)
	hold on;
	bar([0:T],[ M1 optim1_fcn.opt_obj*ones(1,T-1) M1 ],'w')
	for i = 1:rollouts_per_plot
		plot([0:T],exp1.xi_mag_t(:,i))
	end

	xlabel('Time')
	ylabel('\infty Norm of the Estimation Error')
	legend('Guarantees')
	title('Estimator''s Error when ALL Data is available')
	% xt = get(gca, 'XTick');
	% set(gca, 'FontSize', 16)

	subplot(1,2,2)
	hold on;
	bar([0:T],[ M1 optim2_fcn.opt_obj*ones(1,T-1) M1 ],'w')
	for i = 1:rollouts_per_plot
		plot([0:T],exp1.xi2_mag_t(:,i))
	end

	xlabel('Time')
	ylabel('\infty Norm of the Estimation Error')
	legend('Guarantees')
	title('Estimator''s Error when ALL Data is available')
	% xt = get(gca, 'XTick');
	% set(gca, 'FontSize', 16)


	disp('Feedback Gain (for Nominal function):')
	disp(contr1_fcn.F)

	disp('Feedback Gain (for given missing data pattern):')
	disp(contr2_fcn.F)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Experiment 2: Testing the effect of changing the pattern in function %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	disp('Experiment 2: Testing the handling of different pattern inputs')
	disp(' ')

	exp2.pattern1 = [1 0 0 0 1 1];

	[exp2.contr1,exp2.optim1] = create_fhae_w_eq_recovery('min_M2',M1,T,test_sys,'verbosity',2,'pattern',exp1.pattern1);

	[exp2.contr2,exp2.optim2] = create_fhae_w_eq_recovery('min_M2',M1,T,test_sys,'verbosity',2,'pattern',exp2.pattern1)

	[exp2.xi_t,  exp2.xi_mag_t] = apply_controller_to_rollouts(test_sys,exp2.contr2,T,num_rollouts,M1,'missing',[1]);

	figure;
	hold on;
	bar([0:T],[ M1 exp2.optim2.opt_obj*ones(1,T-1) M1 ],'w')
	for i = 1:rollouts_per_plot
		plot([0:T],exp2.xi_mag_t(:,i))
	end

	xlabel('Time')
	ylabel('\infty Norm of the Estimation Error')
	legend('Guarantees')
	title('Estimator''s Error when ALL Data is available')


	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = test_sys;
	results.experim_params.T = T;
	results.experim_params.M1 = M1;

	results.implicit.optimization = optim1_fcn;
	results.implicit.controller = contr1_fcn;

	results.explicit.optimization = optim2_fcn;
	results.explicit.controller   = contr2_fcn;

	results.e2 = exp2;

end
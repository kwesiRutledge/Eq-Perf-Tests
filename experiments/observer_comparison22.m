function [ results ] = observer_comparison22( varargin )
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
	disp('Experiment 1: Generalizability of training')
	disp(' ')
	disp('Given: Pattern 101111.')
	disp('Testing: PAttern 1011.')
	disp(['The standard function is believed to assume the following data pattern for n = ' num2str(n) ' is:' ])

	pattern1 = [1 0 1 1 1 1];
	num_rollouts = 10^4;
	rollouts_per_plot = 100;

	[contr2,optim2] = create_fhae_w_eq_recovery('min_M2',M1,T,test_sys,'verbosity',2,'pattern',pattern1)

	test_sys.x0 = zeros(n,1);
	[xi_t,  xi_mag_t] = apply_controller_to_rollouts(test_sys,contr2,T,num_rollouts,M1,'missing',[1]);

	figure;
	hold on;
	bar([0:T],[ M1 optim2.opt_obj*ones(1,T-1) M1 ],'w')
	for i = 1:rollouts_per_plot
		plot([0:T],xi_mag_t(:,i))
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

	results.optimization = optim2;
	results.controller   = contr2;

end
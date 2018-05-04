function [ results ] = observer_comparison28( varargin )
%	observer_comparison27.m
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
			1,1,1,1,0,1 ];

	%Using ACC System
	load('data/system_examples/acc_p.mat');

	%Create AFF_Dyn object with the data from acc_e
	acc_e = acc;
	acc_e.B = eye(size(acc.A,1));

	acc_ad = Aff_Dyn(acc_e.A,acc_e.B,zeros(size(acc_e.A,1),1), acc_e.C, acc_e.d , acc_e.m, acc_e.E , eye(size(acc.C,1)) );
	[ oc27_opt1 , oc27_contr1 ] = eq_rec_design_t( acc_ad , 'Min_M2' , M1 , T , L );

	disp('Testing a feasibility problem that is known to be feasible.')
	[ oc27_feas_opt1 , oc27_feas_contr1 ] = eq_rec_design_t( acc_ad , 'Feasible Set' , M1 , oc27_contr1.opt_obj, T , L );

	disp('Testing a feasibility problem for the prefix-based feedback law version.')
	[ opt_data2 , new_contr ] = eq_rec_design_pb( acc_ad , 'Feasible Set' , M1 , oc27_contr1.opt_obj, T , L );

	new_contr.apply_control([1,1,0],ones(3*size(acc.C,1),1))

	%%%%%%%%%%%%%%%%%%%%
	%% Saving Results %%
	%%%%%%%%%%%%%%%%%%%%

	results.sys = acc_ad;
	results.L = L;

	results.oc27_opt1 = oc27_opt1;
	results.oc27_contr1 = oc27_contr1;
	
	results.oc27_feas_opt1 = oc27_feas_opt1;
	results.oc27_feas_contr1 = oc27_feas_contr1;

	results.pb_exp1_opt = opt_data2;
	results.pb_exp1_contr = new_contr;

end


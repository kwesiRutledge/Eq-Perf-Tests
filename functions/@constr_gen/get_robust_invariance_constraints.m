function [ varargout ] = get_robust_invariance_constraints(varargin)
	%Description:
	%	Gets the variables necessary to define a finite horizon reachability problem.
	%	- Q-Parameterized versions of the feedback gains (F,f) [Notation used is the NAHS submission variety.]
	%	- Dual variable that satisfy polytope inclusion of the state
	%	- Dual variable that is used to define satisfaction of the input constraint.
	%
	%Usage:
	%	[ Pi2 , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,P_des,Q,r)
	%	[ Pi2 , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,eta_x0,Q,r)
	%
	%	[ Pi1 , Piu , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,P_des,Q,r, 'P_u' , P_u)
	%	[ Pi1 , Piu , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,P_des,Q,r, 'P_u' , P_u , 'u_des' , u_d)
	%
	%	[ Pi1 , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,P_des,gain1,gain2,'param_type',param_flag)
	%	[ Pi1 , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,P_des,Q,r,'param_type','Q')
	%	[ Pi1 , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,P_des,F,f,'param_type','F1',L,xhat0)
	%
	%Inputs:
	%	lcsas - A language constrained switched affine system (LCSAS) defined in an LCSAS object.
	%	word - A word that is feasible/allowed given the definition of the lcsas.
	%	param_flag 	- This flag describes the type of feedback parameterization that is being used.
	%				  (1) 'Q' indicates that output feedback is being used and that the system is using the nonlinear mapping between Q and F.
	%				  (2) 'F1' indicates that disturbance feedback is being used with a very simple estimator in addition.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 6
		error('Not enough input arguments!')
	end

	cg = varargin{1};
	lcsas = varargin{2};
	word = varargin{3};
	P_x0 = varargin{4};
	P_des = varargin{5};

	if ~(isa(P_x0,'Polyhedron') && isa(P_des,'Polyhedron'))
		error('P_x0 and P_des must be given as a Polyhedron.')
	end

	varargin_idx = 8;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case 'P_u'
				P_u = varargin{varargin_idx+1};
				varargin_idx = varargin_idx + 2;
			case 'u_des'
				u_d = varargin{varargin_idx+1};
				if size(u_d,1) ~= length(word)*size(lcsas(1).B,2)
					error('u_d appears to be improperly sized. Please give the full trajectory of desired inputs.')
				end
				varargin_idx = varargin_idx + 2;
			case 'param_type'
				param_flag = varargin{varargin_idx+1};
				switch param_flag
					case 'Q'
						varargin_idx = varargin_idx + 2;
					case 'F1'
						L = varargin{varargin_idx+2};
						xhat0 = varargin{varargin_idx+3};
						varargin_idx = varargin_idx + 4;
					otherwise
						error(['Unexpected parameterization flag: ' param_flag])
				end
			otherwise
				error('Unexpected extra input.')
		end
	end

	if ~exist('u_d')
		u_d = zeros(length(word)*size(lcsas.Dyn(1).B,2),1);
	end

	if ~exist('param_flag')
		param_flag = 'Q';
	end

	% Create Gain Variables
	switch param_flag
	case 'Q'
		Q = varargin{6};
		r = varargin{7};
	case 'F1'
		F = varargin{6};
		f = varargin{7};
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);
	p = size(lcsas.Dyn(1).C,1);
	wd = size(lcsas.Dyn(1).B_w,2);
	vd = size(lcsas.Dyn(1).C_v,2);

	T = length(word);

	q_x0 = size(P_x0.A,1);
	q_des = size(P_des.A,1);

	select_m = @(t,T_r) [zeros(n,t*n), eye(n), zeros(n,(T_r-t)*n) ];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define the Dual Variables for Robust Reachability %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	constraints = [];

	P_wT = 1; P_vT = 1; P_desT = P_des;
	for symb_idx = 1:length(word)
		P_wT = P_wT * lcsas.Dyn( word(symb_idx) ).P_w;
		P_vT = P_vT * lcsas.Dyn( word(symb_idx) ).P_v;
		P_desT = P_desT * P_des;
	end

	P_eta = P_wT * P_vT * P_x0;

	%Get Special Matrices
	[H0,S0,Cm0,J0,k_bar,B_w_big,C_v_big] = lcsas.get_mpc_matrices('word',word);

	switch param_flag
		case 'Q'

			G = [ 	(eye(n*(T+1))+S0*Q*Cm0)*H0*B_w_big ...
					S0*Q*C_v_big ...
					(eye(n*(T+1))+S0*Q*Cm0)*J0 ];

			[Pi2,temp_constrs] = cg.get_H_polyt_inclusion_constr( P_eta.A, P_eta.b , P_desT.A*G, P_des.b-P_desT.A*select_m(T,T)*(S0*r+(eye(n*(T+1))+S0*Q*Cm0)*H0*k_bar) );

			constraints = constraints + temp_constrs;

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Define the Dual Variables for Input Constraint %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			if exist('P_u')

				P_uT = 1;
				for symb_idx = 1:length(word)
					P_uT = P_uT * P_u;
				end

				%Create the G matrix for the input, H_u
				G_u = [ Q*Cm0*H0*B_w_big, Q*C_v_big, Q*Cm0*J0 ];

				[Piu, temp_constrs] = cg.get_H_polyt_inclusion_constr( P_eta.A, P_eta.b , P_uT.A*G_u, P_uT.b- P_uT.A*(r+Q*Cm0*H0*k_bar+u_d) );

				constraints = constraints + temp_constrs;

			end

		


	end

	%%%%%%%%%%%%%%%%%%%%%%%
	%% Output Processing %%
	%%%%%%%%%%%%%%%%%%%%%%%

	varargout_idx = 1;

	if exist('Pi2')
		varargout{varargout_idx} = Pi2;
		varargout_idx = varargout_idx + 1;
	end

	if exist('Piu')
		varargout{varargout_idx} = Piu;
		varargout_idx = varargout_idx + 1;
	end

	varargout{varargout_idx} = constraints;

end
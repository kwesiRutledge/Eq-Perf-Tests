function varargout = get_robust_reachability_constraints(varargin)
	%Description:
	%	Gets the variables necessary to define a finite horizon reachability problem.
	%	- Q-Parameterized versions of the feedback gains (F,f) [Notation used is the NAHS submission variety.]
	%	- Dual variable that satisfy polytope inclusion of the state
	%	- Dual variable that is used to define satisfaction of the input constraint.
	%
	%Allowed Flags:
	%	- 'P_des'
	%	- 'P_u'
	%	- 'u_des'
	%	- 'selection_variable'
	%	- 'param_type'
	%
	%Usage:
	%	[ Pi1 , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,Q,r)
	%	[ Pi1 , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,Q,r,'P_des',P_des)
	%	[ Pi1 , Piu , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,Q,r,'P_des',P_des, 'P_u' , P_u)
	%	[ Pi1 , Piu , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,Q,r,'eta_des',M3, 'P_u' , P_u)
	%	[ Pi1 , Piu , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,Q,r,'P_des',P_des, 'P_u' , P_u , 'u_des' , u_d)
	%	[ Pi1 , Piu , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,Q,r,'P_des',P_des, 'P_u' , P_u , 'selection_variable' , selection_variable)
	%
	%	[ Pi1 , Piu , constraints , bv1 ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,Q,r,'P_des',P_des, 'P_u' , P_u , 'u_des' , u_d)
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

	[ cg , lcsas , word , P_x0 , P_u , P_des , u_des , L , xhat0 , selection_variable , LinGain , OffsetGain , settings_struct ] = grrc_input_processing(varargin{:});

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n , m , p , wd , vd ] = lcsas.Dimensions();

	T = length(word);

	q_x0 = size(P_x0.A,1);

	select_m = @(t,T_r) [zeros(n,t*n), eye(n), zeros(n,(T_r-t)*n) ];

	switch settings_struct.GainFormat
		case 'Q'
			Q = LinGain;
			r = OffsetGain;
		case 'F1'
			F = LinGain;
			f = OffsetGain;
		otherwise
			error(['Unexpected GainFormat value: ' settings_struct.GainFormat ])
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Use Helper Function If Possible %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% If the target set, P_target, is a Polyhedron, then use the simple helper function get_robust_reachability_constraints_polytope.m
	switch settings_struct.PDesiredClass
		case 'Polyhedron'
			if length(P_des) == 1
				[ Pi1 , Piu , constraints ] = get_robust_reachability_constraints_polytope( ...
					cg , lcsas , word , LinGain , OffsetGain , selection_variable , settings_struct ...
					); %See below
				varargout = { Pi1 , Piu , constraints };
				return
			else
				error(['There were multiple polyhedra given as P_des!'])
			end
		case 'double'
			[ Pi1 , Piu , constraints ] = get_robust_reachability_constraints_infnormbound( ...
				cg , lcsas , word , LinGain , OffsetGain , selection_variable , settings_struct ...
				);
			varargout = { Pi1 , Piu , constraints };
			return
		otherwise
			body
	end
	

	%%%%%%%%%%%%%%%%%%%%%%%
	%% Computing Outputs %%
	%%%%%%%%%%%%%%%%%%%%%%%

	varargout{1} = Pi1;

	argout_idx = 2;

	if ~exist('Piu')
		varargout{argout_idx} = Piu;
		argout_idx = argout_idx + 1;
	end

	varargout{argout_idx} = constraints;

end

function [ cg , lcsas , word , P_x0 , P_u , P_des_info , u_des , L , xhat0 , selection_variable , LinGain , OffsetGain , settings_struct ] = grrc_input_processing(varargin)
	%Description:
	%	Processes the inputs received by get_robust_reachability_constraints()

	% Initialize Settings Struct 

	settings_struct = [];
	settings_struct.IncorporateBinaryVariable = false;
	settings_struct.EtaDesGiven = false;

	% Check the number of arguments

	if nargin < 6
		error('Not enough input arguments!')
	end

	cg = varargin{1};

	if isa(varargin{2},'LCSAS')
		lcsas = varargin{2};
	elseif isa(varargin{2},'Aff_Dyn')
		error('This function currently does not support Aff_Dyn objects. Maybe in the future?')
	else
		error('Unexpected input for the second input. Expected an LCSAS or Aff_Dyn type of object.')
	end

	word = varargin{3};
	P_x0 = varargin{4};
	lcsas.X0 = P_x0;

	% Check if the polyhedron P_x0 is of the right type.
	if ~isa(P_x0,'Polyhedron')
		error('P_x0 and P_des must be given as a Polyhedron.')
	end

	varargin_idx = 7;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
			case 'P_u'
				P_u = varargin{varargin_idx+1};
				varargin_idx = varargin_idx + 2;
			case 'P_des'
				P_des = varargin{varargin_idx+1};
				settings_struct.PDesiredClass = class(P_des);
				P_des_info = P_des;
				varargin_idx = varargin_idx + 2;
			case 'eta_des'
				eta_des = varargin{varargin_idx+1};
				settings_struct.PDesiredClass = class(eta_des);
				P_des_info = eta_des;
				varargin_idx = varargin_idx + 2;
			case 'u_des'
				u_des = varargin{varargin_idx+1};
				if size(u_d,1) ~= length(word)*size(lcsas(1).B,2)
					error('u_d appears to be improperly sized. Please give the full trajectory of desired inputs.')
				end
				varargin_idx = varargin_idx + 2;
			case 'param_type'
				settings_struct.GainFormat = varargin{varargin_idx+1};
				switch settings_struct.GainFormat
					case 'Q'
						varargin_idx = varargin_idx + 2;
					case 'F1'
						L = varargin{varargin_idx + 2};
						xhat0 = varargin{varargin_idx + 3};
						varargin_idx = varargin_idx + 4;
					otherwise
						error(['Unexpected parameterization flag: ' param_flag])
				end
			case 'selection_variable'
				selection_variable = varargin{varargin_idx + 1};
				settings_struct.IncorporateBinaryVariable = true;
				varargin_idx = varargin_idx + 2;
			otherwise
				error('Unexpected extra input.')
		end
	end

	if ~exist('u_des')
		u_des = zeros(length(word)*lcsas.Dim_u(),1);
	end

	if ~isfield(settings_struct,'GainFormat')
		settings_struct.GainFormat = 'Q';
	end

	if ~(exist('P_des') || exist('eta_des'))
		P_des = P_x0;
		settings_struct.PDesiredClass = class(P_des);
		P_des_info = P_des;
    end
    
    if ~exist('L')
        L = Language();
    end

    if ~exist('xhat0')
    	xhat0 = zeros(lcsas.Dim_x(),1);
    end

    if ~exist('selection_variable')
    	selection_variable = 0;
    end

	% Create Gain Variables

	switch settings_struct.GainFormat
	case 'Q'
		LinGain = varargin{5};		%Q
		OffsetGain = varargin{6};	%r
	case 'F1'
		LinGain = varargin{5}; 		%F
		OffsetGain = varargin{6}; 	%f
	end

end

function [ Pi1 , Piu , constraints ] = get_robust_reachability_constraints_polytope( cg , lcsas , word , LinGain , OffsetGain , selection_variable , settings_struct )
	%Description:
	%	Gets the variables necessary to define a finite horizon reachability problem.
	%	- Q-Parameterized versions of the feedback gains (F,f) [Notation used is the NAHS submission variety.]
	%	- Dual variable that satisfy polytope inclusion of the state
	%	- Dual variable that is used to define satisfaction of the input constraint.
	%
	%Usage:
	%	[ Pi1 , constraints ] = cg.get_robust_reachability_constraints(lcsas,word,P_x0,Q,r)
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

	if isempty(lcsas.X0)
		error(['Must include complete LCSAS object with well-defined X0 in call to get_robust_reachability_constraints_polytope.'])
	end

	P_x0 = lcsas.X0;

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n , m , p , wd , vd ] = lcsas.Dimensions();

	T = length(word);

	q_x0 = size(P_x0.A,1);

	select_m = @(t,T_r) [zeros(n,t*n), eye(n), zeros(n,(T_r-t)*n) ];

	switch settings_struct.GainFormat
		case 'Q'
			Q = LinGain;
			r = OffsetGain;
		case 'F1'
			F = LinGain;
			f = OffsetGain;
		otherwise
			error(['Unexpected GainFormat value: ' settings_struct.GainFormat ])
	end

	if ~(exist('P_des') || exist('eta_des'))
		P_des = P_x0;
	end

	if ~exist('M')
		M = 10^3;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Define the Dual Variables for Robust Reachability %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	constraints = [];
	Piu = [];

	P_wT = 1; P_vT = 1;
	for symb_idx = 1:length(word)
		P_wT = P_wT * lcsas.Dyn( word(symb_idx) ).P_w;
		P_vT = P_vT * lcsas.Dyn( word(symb_idx) ).P_v;
	end

	P_eta = P_wT * P_vT * P_x0;

	%Get Special Matrices
	[H0,S0,Cm0,J0,k_bar,B_w_big,C_v_big] = lcsas.get_mpc_matrices('word',word);

	switch settings_struct.GainFormat
		case 'Q'
			
			G = [ 	(eye(n*(T+1))+S0*Q*Cm0)*H0*B_w_big ...
					S0*Q*C_v_big ...
					(eye(n*(T+1))+S0*Q*Cm0)*J0 ];

			[Pi1,temp_constrs] = cg.get_H_polyt_inclusion_constr( ...	
										P_eta.A, P_eta.b , ...
										P_des.A*select_m(T,T)*G, ...
										P_des.b + selection_variable*M*ones(size(P_des.b)) - P_des.A*select_m(T,T)*(S0*r+(eye(n*(T+1))+S0*Q*Cm0)*H0*k_bar) );

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

		case 'F1'

			if settings_struct.IncorporateBinaryVariable
				error('F1 flag was not tested for use in this function.')
			end

			L_big = kron(eye(T),L);
			temp_IpLC = (eye(size(H0,1))+ L_big*Cm0)^(-1);

			G = [ 	(S0*F*temp_IpLC+eye(m*T))*H0*B_w_big, ...
					(S0*F*temp_IpLC*(-L_big)+S0*F)*C_v_big, ...
					(S0*F*temp_IpLC+eye(m*T))*J0 ];

			[Pi1,temp_constrs] = cg.get_H_polyt_inclusion_constr( 	P_eta.A, ...
																	P_eta.b , ...
																	P_des.A*select_m(T,T)*G, ...
																	P_des.b-P_des.A*select_m(T,T)*(-S0*F*Cm0*temp_IpLC*J0*xhat0 + S0*f + H0*k_bar ) );

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
				%temp_IpLC = (eye(size(H0,1))+ L_big*Cm0)^(-1);
				G_u = [ F*Cm0*temp_IpLC*B_w_big, F*Cm0*temp_IpLC*(-L_big)*C_v_big, F*Cm0*temp_IpLC*J0 ];

				[Piu, temp_constrs] = cg.get_H_polyt_inclusion_constr( 	P_eta.A, P_eta.b , ...
																		P_uT.A*G_u, P_uT.b- P_uT.A*(f-F*Cm0*temp_IpLC*J0*xhat0) );

				constraints = constraints + temp_constrs;

			end


	end

end

function [ Pi1 , Piu , constraints ] = get_robust_reachability_constraints_infnormbound( cg , lcsas , word , LinGain , OffsetGain , selection_variable , settings_struct )
	%Description:
	%	Handles the case when P_des is not given as a 

	if isempty(lcsas.X0)
		error(['Must include complete LCSAS object with well-defined X0 in call to get_robust_reachability_constraints_polytope.'])
	end

	P_x0 = lcsas.X0;

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	[ n , m , p , wd , vd ] = lcsas.Dimensions();

	T = length(word);

	q_x0 = size(P_x0.A,1);

	select_m = @(t,T_r) [zeros(n,t*n), eye(n), zeros(n,(T_r-t)*n) ];

	switch settings_struct.GainFormat
		case 'Q'
			Q = LinGain;
			r = OffsetGain;
		case 'F1'
			F = LinGain;
			f = OffsetGain;
		otherwise
			error(['Unexpected GainFormat value: ' settings_struct.GainFormat ])
	end

	if ~(exist('P_des') || exist('eta_des'))
		P_des = P_x0;
	end

	if ~exist('M')
		M = 10^3;
	end

	% Construct COnstraints

	constraints = [];
	Piu = [];

	P_wT = 1; P_vT = 1;
	for symb_idx = 1:length(word)
		P_wT = P_wT * lcsas.Dyn( word(symb_idx) ).P_w;
		P_vT = P_vT * lcsas.Dyn( word(symb_idx) ).P_v;
	end

	P_eta = P_wT * P_vT * P_x0;

	%Get Special Matrices
	[H0,S0,Cm0,J0,k_bar,B_w_big,C_v_big] = lcsas.get_mpc_matrices('word',word);

	switch param_flag
		case 'Q'
			if exist('P_des')
				G = [ 	(eye(n*(T+1))+S0*Q*Cm0)*H0*B_w_big ...
						S0*Q*C_v_big ...
						(eye(n*(T+1))+S0*Q*Cm0)*J0 ];

				if settings_struct.IncorporateBinaryVariable
					[Pi1,temp_constrs] = cg.get_H_polyt_inclusion_constr( P_eta.A, P_eta.b , P_des.A*select_m(T,T)*G, P_des.b-P_des.A*select_m(T,T)*(S0*r+(eye(n*(T+1))+S0*Q*Cm0)*H0*k_bar) );
				else
					[Pi1,temp_constrs] = cg.get_H_polyt_inclusion_constr( 	P_eta.A, P_eta.b + selection_variable*M*ones(size(P_eta.b)) , ...
																			P_des.A*select_m(T,T)*G, P_des.b-P_des.A*select_m(T,T)*(S0*r+(eye(n*(T+1))+S0*Q*Cm0)*H0*k_bar) );
				end

				constraints = constraints + temp_constrs;
			elseif exist('eta_des')

				G = [ 	(eye(n*(T+1))+S0*Q*Cm0)*H0*B_w_big ...
						S0*Q*C_v_big ...
						(eye(n*(T+1))+S0*Q*Cm0)*J0 ];

				I_t = [eye(n);-eye(n)];

				if settings_struct.IncorporateBinaryVariable
					[Pi1,temp_constrs] = cg.get_H_polyt_inclusion_constr( P_eta.A, P_eta.b , I_t*select_m(T,T)*G, eta_des*ones(2*n,1)-I_t*select_m(T,T)*(S0*r+(eye(n*(T+1))+S0*Q*Cm0)*H0*k_bar) );
				else
					[Pi1,temp_constrs] = cg.get_H_polyt_inclusion_constr( 	P_eta.A, P_eta.b, ...
																			I_t*select_m(T,T)*G, ...
																			eta_des*ones(2*n,1)-I_t*select_m(T,T)*(S0*r+(eye(n*(T+1))+S0*Q*Cm0)*H0*k_bar) + selection_variable*M*ones(2*n,1) );
				end

				constraints = constraints + temp_constrs;
			else
				error('Unrecognized input combination.')
			end

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

				if settings_struct.IncorporateBinaryVariable
					[Piu, temp_constrs] = cg.get_H_polyt_inclusion_constr( P_eta.A, P_eta.b , P_uT.A*G_u, P_uT.b- P_uT.A*(r+Q*Cm0*H0*k_bar+u_d) );
				else
					[Piu, temp_constrs] = cg.get_H_polyt_inclusion_constr( 	P_eta.A, P_eta.b, ...
																			P_uT.A*G_u, ...
																			P_uT.b- P_uT.A*(r+Q*Cm0*H0*k_bar+u_d)  + selection_variable*M*ones(2*n,1) );
				end

				constraints = constraints + temp_constrs;

			end

		case 'F1'

			if exist('selection_variable')
				error('Didn''t expect for selection_variable to be defined for F1')
			end


			L_big = kron(eye(T),L);
			temp_IpLC = (eye(size(H0,1))+ L_big*Cm0)^(-1);

			G = [ 	(S0*F*temp_IpLC+eye(m*T))*H0*B_w_big, ...
					(S0*F*temp_IpLC*(-L_big)+S0*F)*C_v_big, ...
					(S0*F*temp_IpLC+eye(m*T))*J0 ];

			[Pi1,temp_constrs] = cg.get_H_polyt_inclusion_constr( 	P_eta.A, ...
																	P_eta.b , ...
																	P_des.A*select_m(T,T)*G, ...
																	P_des.b-P_des.A*select_m(T,T)*(-S0*F*Cm0*temp_IpLC*J0*xhat0 + S0*f + H0*k_bar ) );

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
				%temp_IpLC = (eye(size(H0,1))+ L_big*Cm0)^(-1);
				G_u = [ F*Cm0*temp_IpLC*B_w_big, F*Cm0*temp_IpLC*(-L_big)*C_v_big, F*Cm0*temp_IpLC*J0 ];

				[Piu, temp_constrs] = cg.get_H_polyt_inclusion_constr( 	P_eta.A, P_eta.b , ...
																		P_uT.A*G_u, P_uT.b- P_uT.A*(f-F*Cm0*temp_IpLC*J0*xhat0) );

				constraints = constraints + temp_constrs;

			end


	end

end
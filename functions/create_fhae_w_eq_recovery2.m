function [ controller , optim_info ] = create_fhae_w_eg_recovery2( varargin )
	% create_fhae_w_eq_recovery.m
	%	Description:
	%		Create a Finite-Horizon Affine Estimator that satisfies Equalized Recovery,		
	%		optionally, with missing data occuring at certain times.
	%
	%	Usage:
	%		[controller,optim_info] = create_fhae_w_eq_recovery( 'feasibility' , M1 , M2 , T , sys )
	%		[controller,optim_info] = create_fhae_w_eq_recovery( 'min_M2' , M1 , T , sys )
	%		[controller,optim_info] = create_fhae_w_eq_recovery( 'min_M2' , M1 , T , sys , 'verbosity', verbosity )
	%		[controller,optim_info] = create_fhae_w_eq_recovery( 'min_M2' , M1 , T , sys , 'pattern' , m_patt)
	%
	%	Inputs:
	%		input_str - 	The first argument to this function should be a string telling the function
	%						what kind of design problem this will be.
	%						If it is just feasibility, then the constants M1,M2 and T are assumed to be given.
	%						
	%		M1 - 			The desired bound for the estimation error at times t = 0 and t = T.
	%		M2 - 			The desired bound for the estimation error for times t in [0,T]. (Assume M1 <= M2)
	%		T - 			Time horizon that we consider for recovery.
	%		sys -
	%		param_name - 	Name of the parameter which we will now tune.
	%		
	%		pattern_array - Array representing the missing data patterns that we will be robust against.
	%						Array (n_patt x n(T+1)). Should contain only 0 and 1 values.
	%						For each row in the matrix, it represents a feasible data pattern.
	%						i.e. if the row contains, [1,0,1,1] that means the following observations are available:
	%						y(t),y(t+2), and y(t+3) while the following observation is missing y(t+1)

	%%%%%%%%%%%%%%%%%%%
	%% Manage Inputs %%
	%%%%%%%%%%%%%%%%%%%

	if nargin <= 3
		error('Need at least 4 inputs.')
	end

	if strcmp(varargin{1},'feasibility')
		if (nargin < 5)
			error(['Need more than 4 inputs. Received ' num2str(nargin) '.' ])
		end

		M1 = varargin{2};
		M2 = varargin{3};
		T  = varargin{4};

		sys = varargin{5};

	elseif strcmp(varargin{1},'min_M2')

		if (nargin < 4)
			error(['Need more than 3 inputs. Received ' num2str(nargin) '.' ] )
		end

		M1 = varargin{2};
		T = varargin{3};

		sys = varargin{4};

		for ind = 5:2:nargin-1
			if strcmp(varargin{ind},'verbosity')
				verbosity = varargin{ind+1};
			elseif strcmp(varargin{ind},'pattern')
				if size(varargin{ind+1},2) > T
					error('Pattern defines more than ''T'' data points. Please choose T, pattern matrices appropriatetly.')
				end
				m_patt = varargin{ind+1};
			else
				error(['Unrecognized param_name string. Input number ' num2str(ind) '.' ])
			end
		end

		if ~exist('m_patt')
			m_patt = ones(1,T);
		end

	end

	sys_args = {'A','B','C','E','F','x0','m','d','C_v'};
	args_exist = [];

	missing_strs = [''];

	for arg_num = 1 : length(sys_args)
		args_exist(arg_num) = isfield(sys,sys_args(arg_num));
		
		if ~args_exist(arg_num) & (arg_num-1 == sum(args_exist))
			missing_strs = [missing_strs sys_args{arg_num}];
		elseif ~args_exist(arg_num)
			missing_strs = [ missing_strs ',' sys_args{arg_num} ];
		end

	end

	%There are some required arguments
	if ~all( args_exist([1 2 3 7 8]) )
		error('Missing some of the necessary system arguments.')
	end

	if ~isempty(args_exist)
		warning(['Missing some system arguments: ' missing_strs ])
	end

	%Fix the missing arguments
	for arg_num = 1 : length(sys_args)
		if (strcmp(sys_args{arg_num},'C_v')) & (~args_exist(arg_num))
			sys.C_v = eye(size(sys.C,1));
		end

		if (strcmp(sys_args{arg_num},'E')) & (~args_exist(arg_num))
			sys.E = eye(size(sys.A,1));
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%
	%% Define Constants %%
	%%%%%%%%%%%%%%%%%%%%%%

	%Defining error systems
	sys.B = eye(size(sys.A,1));

	n = size(sys.A,1);
	m = size(sys.B,2);
	p = size(sys.C,1);
	wd = size(sys.E,2);
	vd = size(sys.C_v,2);

	%Select matrix
	select_m = @(t,T_r) [ zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	ops = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%%%%%
	%% Synthesis %%
	%%%%%%%%%%%%%%%

	if strcmp(varargin{1},'feasibility')
			
		% FEASIBILITY

		% Optimization Variables
		% ++++++++++++++++++++++
		w 		= sdpvar(wd*T,1,'full');
		v 		= sdpvar(vd*T,1,'full');
		sys.x0 	= sdpvar(n,1,'full');

		% Feedback Variables
		Q = sdpvar(m*T,p*T,'full');
		r = sdpvar(m*T,1,'full');

		% Dual Variables
		Pi_1 = sdpvar(2*n*T,2*(wd+p)*T+2*n,'full');
		Pi_2 = sdpvar(2*n,2*(wd+p)*T+2*n,'full');

		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,Bw0,Cv0] = create_skaf_n_boyd_matrices(sys,T);

		positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
		end

		noise_constrs = [ Pi_1 * [ sys.d * ones(2*wd*T,1) ; sys.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
		noise_constrs = noise_constrs + [ Pi_2 * [ sys.d * ones(2*wd*T,1) ; sys.m * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T
			pre_xi = [ pre_xi ; sys.A^i];
		end

		G = [ (eye(n*(T+1))+S0*Q*Cm0)*S0*Bw0 S0*Q (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,p*T+n) ;
									zeros(2*p*T,wd*T) [ eye(p*T) ; -eye(p*T) ] zeros(2*p*T,n) ;
									zeros(2*n,(p+wd)*T) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = [ Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

		%Lower Diagonal Constraint
		l_diag_constr = [];
		for bl_row_num = 1 : T-1
			l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(sys.B,2)+1:bl_row_num*size(sys.B,2)], ...
													[bl_row_num*size(sys.C,1)+1:end] ) == 0 ];
		end

		% OPTIMIZATION
		% ++++++++++++
		optim = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr, ...
				ops)

		if optim.problem ~= 0
			error(['The design problem was not completely solved.' optim1.info ])
		end

	elseif strcmp(varargin{1},'min_M2')
		
		% Optimization Variables
		% ++++++++++++++++++++++
		w 		= sdpvar(wd*T,1,'full');
		v 		= sdpvar(vd*T,1,'full');
		sys.x0 	= sdpvar(n,1,'full');

		alpha_2 = sdpvar(1,1,'full');

		alpha_l = sdpvar(T+1,1,'full');

		% Feedback Variables
		Q = sdpvar(m*T,p*T,'full');
		r = sdpvar(m*T,1,'full');

		% Dual Variables
		Pi_1 = sdpvar(2*n*T,2*(wd+vd)*T+2*n,'full');
		Pi_2 = sdpvar(2*n,2*(wd+vd)*T+2*n,'full');

		% Creating Constraints
		% ++++++++++++++++++++

		[S0,H0,Cm0,xi0m,Bw0,Cv0] = create_skaf_n_boyd_matrices(sys,T);

		positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
		end

		noise_constrs = [ Pi_1 * [ sys.d * ones(2*wd*T,1) ; sys.m * ones(2*vd*T,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
		noise_constrs = noise_constrs + [ Pi_2 * [ sys.d * ones(2*wd*T,1) ; sys.m * ones(2*vd*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T
			pre_xi = [ pre_xi ; sys.A^i];
		end

		G = [ (eye(n*(T+1))+S0*Q*Cm0)*S0*Bw0 S0*Q*Cv0 (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

		bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,p*T+n) ;
									zeros(2*p*T,wd*T) [ eye(p*T) ; -eye(p*T) ] zeros(2*p*T,n) ;
									zeros(2*n,(p+wd)*T) [ eye(n) ; -eye(n) ] ];

		dual_equal_constrs = [ Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

		%Lower Diagonal Constraint
		l_diag_constr = [];
		for bl_row_num = 1 : T-1
			l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*size(sys.B,2)+1:bl_row_num*size(sys.B,2)], ...
													[bl_row_num*size(sys.C,1)+1:end] ) == 0 ];
		end

		% Constraints due to missing data pattern
		for comb_num = 1:size(m_patt,1)

			missing_locs = find(m_patt(comb_num,:)==0)-1;

			%Calculate Big C Matrix
			[~,~,Cm,~] = create_skaf_n_boyd_matrices(sys,T,'missing',missing_locs);

			%Create Special selection matrices for selecting the proper variables
			R = [ zeros(n,n*T) eye(n) ];
				% mu_select = [];
				% for i = 1:T+1
				% 	%Takes out the row representing the observation at time t=i+1
				% 	if any(missing_loc+1+1 == i)
				% 		mu_select = [ mu_select ; zeros(wd+vd,(wd+vd)*(T+1)) ];
				% 	else
				% 		mu_select = [ mu_select ; [ zeros((wd+vd),vd*(i-1)) eye(wd+vd) zeros(wd+vd,p*(T-i)) ] ];
				% 	end
				% end

				% Pxd = [ Pxd ; (eye(n*(T+1))+S*Q*Cm)*S ];%E_bar*[ eye(b_dim*T) zeros(b_dim*T,b_dim) ];
				% Pxm = [ Pxm ; S*Q*mu_select ];
				% xi_tilde = [ xi_tilde ; (eye(n*(T+1)) + S*Q*Cm)*xi0m + S*r];

			pre_v = eye(vd*T);
			for t = [0:T-1]
				if any(missing_locs == t)
					pre_v([t*vd+1:(t+1)*vd],:) = 0; 
				end
			end
			

			G = [ (eye(n*(T+1))+S0*Q*Cm)*S0*Bw0 S0*Q*pre_v (eye(n*(T+1))+S0*Q*Cm)*pre_xi ];

			%Awd to the constraint set
			dual_equal_constrs = dual_equal_constrs + [Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G];
			dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

		end

		% OPTIMIZATION
		% ++++++++++++
		ops = sdpsettings('verbose',verbosity);
		optim = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr, ...
				alpha_2, ...
				ops)

		if optim.problem ~= 0
			error(['The design problem was not completely solved.' optim1.info ])
		end

		% Save some results
		optim.opt_obj = value(alpha_2);

	end

	%%%%%%%%%%%%%
	%% Outputs %%
	%%%%%%%%%%%%%

	Q1 = value(Q);
	r1 = value(r);
	F1 = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
	u0_1 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

	controller.F = F1;
	controller.u0 = u0_1;
	
	optim_info = optim;
	optim_info.Q = Q1;
	optim_info.r = r1;

end
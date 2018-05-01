function [ opt_out, fb ] = eq_rec_design_t( varargin )
%Description:
%	Searches for a feasible, time-based feedback that satisfies the parameters given in the
%	Equalized Recovery Problem:
%		(M1,M2,T)
%
%	There are 2 different problems that we consider:
%	'Min_M2' , 'Feasible Set'. Which describe the purpose of our optimization.
%
%Usage:
%	eq_rec_design_tf( ad , 'Feasible Set' , M1 , M2 , T )	
%	eq_rec_design_tf( ad , 'Min_M2' , M1 , T )
%	eq_rec_design_tf( ad , 'Feasible Set' , M1 , M2 , T , L )
%	eq_rec_design_tf( ad , 'Min_M2' , M1 , T , L)

if nargin < 2
	error('Not enough initial inputs given.')
end

%Manage Inputs
ad 		= varargin{1};
str_in 	= varargin{2}; 

% Constants
n = size(ad.A,1);
m = size(ad.B,2);
p = size(ad.C,1);
wd = size(ad.B_w,2);
vd = size(ad.C_v,2);

verbosity = 1;
ops = sdpsettings('verbose',verbosity);

%Select matrix
select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

switch str_in
case 'Feasible Set'
	%Collect Inputs
	if nargin < 4
		error('Not enough inputs.')
	end

	M1 = varargin{3};
	M2 = varargin{4};
	T  = varargin{5};

	if nargin > 5
		L = varargin{6};
	else
		L = ones(1,T);
	end

	% Create L_star
	L_star = ones(1,T);
	for sig_i = 1:size(L,1)
		L_star = bitand(L_star,L(sig_i,:));
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform Optimization %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Optimization Variables
	% ++++++++++++++++++++++
	w		= sdpvar(wd*T,1,'full');
	v		= sdpvar(vd*T,1,'full');
	ad.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');
	alpha_l 	= sdpvar(T+1,1,'full');

	% Feedback Variables
	Q = sdpvar(n*T,p*T,'full');
	r = sdpvar(n*T,1,'full');

	% Dual Variables
	Pi_1 = sdpvar(2*n*T,2*(wd+vd)*T+2*n,'full');
	Pi_2 = sdpvar(2*n,2*(wd+vd)*T+2*n,'full');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Solving for Worst Case Word while overconstraining Q %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[S0,H0,Cm0,xi0m,B_w_big] = create_skaf_n_boyd_matrices(ad,T,'missing',find(L_star==0)-1);

	positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
	end

	noise_constrs = [ Pi_1 * [ ad.eta_w * ones(2*wd*T,1) ; ad.eta_v * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
	noise_constrs = noise_constrs + [ Pi_2 * [ ad.eta_w * ones(2*wd*T,1) ; ad.eta_v * ones(2*p*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T
		pre_xi = [ pre_xi ; ad.A^i];
	end

	G = [ (eye(n*(T+1))+S0*Q*Cm0)*S0*B_w_big S0*Q (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

	bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,vd*T+n) ;
								zeros(2*vd*T,wd*T) [ eye(vd*T) ; -eye(vd*T) ] zeros(2*vd*T,n) ;
								zeros(2*n,(vd+wd)*T) [ eye(n) ; -eye(n) ] ];

	dual_equal_constrs = [ Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

	%Lower Diagonal Constraint
	l_diag_constr = [];
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
												[bl_row_num*p+1:end] ) == 0 ];
	end

	% Disturbance v=0
	v_constr = [];
	for v_ind = find(L_star==0)
		v_constr = v_constr + [ v([(v_ind-1)*p+1:v_ind*p])==0 ];
	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim0 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr + v_constr, ...
			[], ...
			ops);

	opt_out = optim0;
	if opt_out.problem ~= 0
		fb = [];
	else
		% Save Feedback Matrices
		% ++++++++++++++++++++++
		fb.Q  = value(Q);
		fb.r  = value(r);
		fb.F  = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
		fb.u0 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

		% fb.opt_obj = value(alpha_2);
	end

case 'Min_M2'
	%Collect Inputs
	if nargin < 4
		error('Not enough inputs.')
	end

	M1 = varargin{3};
	T  = varargin{4};

	if nargin > 4
		L = varargin{5};
	else
		L = ones(1,T);
	end

	% Create L_star
	L_star = ones(1,T);
	for sig_i = 1:size(L,1)
		L_star = bitand(L_star,L(sig_i,:));
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform Optimization %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Optimization Variables
	% ++++++++++++++++++++++
	w		= sdpvar(wd*T,1,'full');
	v		= sdpvar(vd*T,1,'full');
	ad.x0 	= sdpvar(n,1,'full');

	alpha_2 	= sdpvar(1,1,'full');
	alpha_l 	= sdpvar(T+1,1,'full');

	% Feedback Variables
	Q = sdpvar(n*T,p*T,'full');
	r = sdpvar(n*T,1,'full');

	% Dual Variables
	Pi_1 = sdpvar(2*n*T,2*(wd+vd)*T+2*n,'full');
	Pi_2 = sdpvar(2*n,2*(wd+vd)*T+2*n,'full');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Solving for Worst Case Word while overconstraining Q %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(ad,T,'missing',find(L_star==0)-1);

	positive_constr = [ Pi_1 >= 0, Pi_2 >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T) ];
	end

	noise_constrs = [ Pi_1 * [ ad.eta_w * ones(2*wd*T,1) ; ad.eta_v * ones(2*vd*T,1) ; M1 * ones(2*n,1) ] <= alpha_2 * ones(2*n*T,1) - [eye(n*T);-eye(n*T)]*sel_influenced_states*S0*r ];
	noise_constrs = noise_constrs + [ Pi_2 * [ ad.eta_w * ones(2*wd*T,1) ; ad.eta_v * ones(2*vd*T,1) ; M1 * ones(2*n,1) ] <= M1 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T,T)*S0*r ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T
		pre_xi = [ pre_xi ; ad.A^i];
	end

	G = [ (eye(n*(T+1))+H0*Q*Cm0)*S0*B_w_big H0*Q*C_v_big (eye(n*(T+1))+S0*Q*Cm0)*pre_xi ];

	bounded_disturb_matrix = [ [ eye(wd*T) ; -eye(wd*T) ] zeros(2*wd*T,vd*T+n) ;
								zeros(2*vd*T,wd*T) [ eye(vd*T) ; -eye(vd*T) ] zeros(2*vd*T,n) ;
								zeros(2*n,(vd+wd)*T) [ eye(n) ; -eye(n) ] ];

								size(Pi_1)
								size(bounded_disturb_matrix)
								size(G)

	dual_equal_constrs = [ Pi_1 * bounded_disturb_matrix == [eye(n*T); -eye(n*T)]*sel_influenced_states*G ];
	dual_equal_constrs = dual_equal_constrs + [Pi_2 * bounded_disturb_matrix == [eye(n);-eye(n)]*select_m(T,T)*G];

	%Lower Diagonal Constraint
	l_diag_constr = [];
	for bl_row_num = 1 : T-1
		l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*m+1:bl_row_num*m], ...
												[bl_row_num*p+1:end] ) == 0 ];
	end

	% Disturbance v=0
	v_constr = [];
	for v_ind = find(L_star==0)
		v_constr = v_constr + [ v([(v_ind-1)*p+1:v_ind*p])==0 ];
	end

	% OPTIMIZE
	% ++++++++

	% ops = sdpsettings('verbose',verbosity);
	optim0 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr + v_constr, ...
			alpha_2, ...
			ops);

	opt_out = optim0;
	if opt_out.problem ~= 0
		fb = [];
	else
		% Save Feedback Matrices
		% ++++++++++++++++++++++
		fb.Q  = value(Q);
		fb.r  = value(r);
		fb.F  = value( (inv(value(eye(size(Q,1)) + Q*Cm0*S0)) ) * Q);
		fb.u0 = value( inv(value(eye(size(Q,1)) + Q*Cm0*S0)) * r );

		fb.opt_obj = value(alpha_2);
	end

otherwise
	error(['Unrecognized String: ' str_in ] )

end


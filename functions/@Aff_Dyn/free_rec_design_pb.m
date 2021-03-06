function [ opt_out, contr ] = free_rec_design_pb( varargin )
%Description:
%	Searches for a feasible, prefix-based feedback that satisfies the parameters given in the
%	Equalized Recovery Problem:
%		(M1,M2,M3,T,L)
%
%	There are 3 different problems that we consider:
%		1. 'Min_M2'
%		2. 'Feasible Set'
%		3. 'Min_M3'
%	Which describe the purpose of our optimization.
%
%Usage:
%	[ opt_out, contr ] = free_rec_design_pb( ad , 'Feasible Set' , M1 , M2 , M3 , T )
%	[ opt_out, contr ] = free_rec_design_pb( ad , 'Feasible Set' , M1 , M2 , M3 , L )	
%	[ opt_out, contr ] = free_rec_design_pb( ad , 'Min_M2' , M1 , M3 , T )
%	[ opt_out, contr ] = free_rec_design_pb( ad , 'Min_M2' , M1 , M3 , L )
%	[ opt_out, contr ] = free_rec_design_pb( ad , 'Min_M3' , M1 , M2 , T )
%	[ opt_out, contr ] = free_rec_design_pb( ad , 'Min_M3' , M1 , M2 , L )
%	[ opt_out, contr ] = free_rec_design_pb( ad , 'Min_M3' , M1 , M2 , T , 'verbosity' , 0 )
%
%Notes:
%	Previously, a 'Max_M3' mode was tried, but the problem is unbounded. Typically any M3 larger than the minimum M3 will do for this problem.

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

%Select matrix
select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

nargin_evaltd = -1;

%Constraints on Objective
obj_constrs = [];
switch str_in
case 'Feasible Set'
	%Collect Inputs
	if nargin < 6
		error('Not enough inputs.')
	end

	M1 = varargin{3};
	M2 = varargin{4};
	M3 = varargin{5};

	if iscell(varargin{6})
		L = varargin{6};
	elseif isscalar(varargin{6})
		L = {ones(1,varargin{6})};
	elseif isnumeric(varargin{6})
		L_in = varargin{6};
		L = {};
		for ind = 1:size(L_in,1)
			L{ind} = L_in(ind,:);
		end
	else
		error(['Unrecognized fifth input: ' num2str(varargin{6}) '.'])
	end

	%Define Objective Function
	obj_fcn = [];

	%Update the number of inputs that have been evaluated.
	nargin_evaltd = 6;

case 'Min_M1'

	if nargin < 5
		error('Not enough inputs for Min_M1 mode.')
	end

	M2 = varargin{3};
	M3 = varargin{4}
	if iscell(varargin{5})
		L = varargin{5};
	elseif isscalar(varargin{5})
		L = {ones(1,varargin{5})};
	elseif isnumeric(varargin{5})
		L_in = varargin{5};
		L = {};
		for ind = 1:size(L_in,1)
			L{ind} = L_in(ind,:);
		end
	else
		error(['Unrecognized fifth input: ' num2str(varargin{5}) '.']);
	end

	M1 = sdpvar(1,1,'full');

	obj_fcn = -M1;
	obj_constrs = [obj_constrs,M1 >= 0];

	%Update the number of inputs that have been evaluated.
	nargin_evaltd = 5;

	error('This problem is nonlinear and does not currently have a reasonable solution in our framework.')

case 'Min_M2'

	if nargin < 5
		error('Not enough inputs for Min_M2 mode.')
	end

	M1 = varargin{3};
	M3 = varargin{4};
	if iscell(varargin{5})
		L = varargin{5};
	elseif isscalar(varargin{5})
		L = {ones(1,varargin{5})};
	elseif isnumeric(varargin{5})
		L_in = varargin{5};
		L = {};
		for ind = 1:size(L_in,1)
			L{ind} = L_in(ind,:);
		end
	else
		error(['Unrecognized fifth input: ' num2str(varargin{5}) '.']);
	end

	M2 = sdpvar(1,1,'full');
	obj_fcn = M2;

	%Update the number of inputs that have been evaluated.
	nargin_evaltd = 5;

case 'Min_M3'

	%Process Inputs
	if nargin < 5
		error('Not enough inputs for Min_M1 mode.')
	end

	%Save Inputs
	M1 = varargin{3};
	M2 = varargin{4};

	if iscell(varargin{5})
		L = varargin{5};
	elseif isscalar(varargin{5})
		L = {ones(1,varargin{5})};
	elseif isnumeric(varargin{5})
		L_in = varargin{5};
		L = {};
		for ind = 1:size(L_in,1)
			L{ind} = L_in(ind,:);
		end
	else
		error(['Unrecognized fifth input: ' num2str(varargin{5}) '.']);
	end

	%Create Optimization Objective
	M3 = sdpvar(1,1,'full');
	obj_fcn = M3;
	obj_constrs = obj_constrs + [M3 >= 0];

	%Update the number of inputs that have been evaluated.
	nargin_evaltd = 5;

otherwise
	error(['Unrecognized String: ' str_in ] )
end

%If the number of arguments evaluated is greater than 0 (valid)
%and LESS than the maximum. Look for other options to tweak.
while (nargin_evaltd > 0) & (nargin > nargin_evaltd)
	switch varargin{nargin_evaltd+1}
		case 'verbosity'
			verbosity = varargin{nargin_evaltd+2};
			nargin_evaltd = nargin_evaltd + 2;
		otherwise
			error(['Unrecognized option: ' varargin{nargin_evaltd} ])
	end

end

if ~exist('verbosity')
	verbosity = 1;
end
ops = sdpsettings('verbose',verbosity);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Optimization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

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

	[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(ad,T_i,'missing',find(L{pattern_ind} == 0)-1);

	positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T_i
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
	end

	noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ ad.eta_w * ones(2*wd*T_i,1) ; ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M2 * ones(2*n*T_i,1) - [eye(n*T_i);-eye(n*T_i)]*sel_influenced_states*H0*r{pattern_ind} ];
	noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ ad.eta_w * ones(2*wd*T_i,1) ; ad.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M3 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T_i,T_i)*H0*r{pattern_ind} ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T_i
		pre_xi = [ pre_xi ; ad.A^i];
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
		[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(ad,T_i,'missing',find(L{pattern_ind} == 0)-1);

		Q_set{pattern_ind} = value(Q{pattern_ind});
		r_set{pattern_ind} = value(r{pattern_ind});
		F_set{pattern_ind} = value( (pinv(value(eye(size(Q{pattern_ind},1)) + Q{pattern_ind}*Cm0*H0)) ) * Q{pattern_ind});
		u0_set{pattern_ind} = value( pinv(value(eye(size(Q{pattern_ind},1)) + Q{pattern_ind}*Cm0*H0)) * r{pattern_ind} );
		
		%Fix up F and u0 to avoid NaN
		F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
		% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

		%Create Function Outputs
		opt_out.Q_set = Q_set;
		opt_out.r_set = r_set;

		contr = FHAE_pb(L,F_set,u0_set);

	end
end

switch varargin{2}
case 'Max_M1'
	opt_out.M1 = value(M1);
case 'Min_M2'
	opt_out.M2 = value(M2);
case 'Min_M3'
	opt_out.M3 = value(M3);
end

end
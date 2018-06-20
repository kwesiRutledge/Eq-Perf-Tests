% tricky_counterexample.m

clear all;
close all;
clc;

%% Import Functions

try
	FSA({1:10}',... %'
		{1},[0;1],[ [1:10] ; mod([1:10],2) ]', ... %'
		[ [1:10] [1:10] ; ones(1,10) zeros(1,10) ; [1:10] mod([1:10]+1,10)+1 ]');
catch
	add('functions/')
end

try
	Polyhedron
catch
	error('MPT3 is not added to path.')
end

%% Constants
n = 3;
m = n;

M1 = 1;
M2 = 2;
eta.w = 0.1;
eta.v = 0.2;

A_eigval = [0.6; 1.1; 0.4];
sp_evec = [1;0.26;1];
sp_evec = sp_evec / norm(sp_evec,2);

A_eigvec = eye(n);%[1 1 0; 0 0.2 0; 0 1 1];
A_eigvec(:,2) = sp_evec;

sys.A = A_eigvec*diag(A_eigval)*inv(A_eigvec);

sys.B = eye(n);
sys.B_w = eye(n);
sys.C = [0 1 0; 0 0 1];

dyn = Aff_Dyn(sys.A,sys.B,zeros(n,1),sys.C,eta.w,eta.v);

experim_enable = [	false;
					true];

%% Attempt a Filter Design for simple problem.
L=[1,1,1,1];
for t = 1:8
	L = [];
	L = ones(1,t);
	[ opt_out, contr ] = free_rec_design_pb( dyn , 'Min_M3' , M1 , M2 , L )
	exp1.results{t}.opt_out = opt_out;
	exp1.results{t}.contr = contr;
	M3(t) = opt_out.M3;
end


figure;
plot(M3)
xlabel('Time Horizon T')
ylabel('Estimation Error')

disp('The Problem is configured such that, with no memory,')
disp('I believe that an estimation error increase will always occur in the first time step.')
disp('From here, I suspect that without memory we will make an improper assumption and ')

%% Problem That Should Not Be Feasible

disp('Challenge Problem')
disp(' ')
disp('Now, using the problem without memory we will see if recovery is feasible on the hard')
disp('pattern:')
disp(' [ R-1-> X -1-> X -1-> R -1-> R ] ')
disp(' ')

L = [1];
[ opt_out, contr ] = free_rec_design_pb( dyn , 'Min_M3' , M3(3) , M2 , L )

disp('Okay!')
disp('With the tricky problem it could not recover well when given the correct hypercube size.')
disp('But the hypercube is an over approximation, let''s use MPT3 to give a better constraint.')

n = size(dyn.A,1);
m = size(dyn.B,2);
p = size(dyn.C,1);
wd = size(dyn.B_w,2);
vd = size(dyn.C_v,2);

%New Inputs
M1 = M3(3);
%M2 unchanged
L = {[1]};

%Create Optimization Objective
M3 = sdpvar(1,1,'full');
obj_fcn = M3;
obj_constrs = [];
obj_constrs = obj_constrs + [M3 >= 0];

verbosity = 1;
ops = sdpsettings('verbose',verbosity);

%Select matrix
select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Polyhedral Set where xi(t_0) really lies %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_disturb_box = Polyhedron('A',[eye(3*wd);-eye(3*wd)],'b',dyn.eta_w*ones(2*(3*wd),1));
v_disturb_box = Polyhedron('A',[eye(3*vd);-eye(3*vd)],'b',dyn.eta_v*ones(2*3*vd,1));
xi_disturb_box = Polyhedron('A',[eye(n);-eye(n)],'b',M1*ones(2*n,1));

%Get the Q and r that correspond to M3(3)
Q_set0 = exp1.results{3}.opt_out.Q_set;
r_set0 = exp1.results{3}.opt_out.r_set;

[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(dyn,3,'missing',find([1,1,1] == 0)-1);

xi_t_polyh = (S0+S0*Q_set0{1}*Cm0*S0)*w_disturb_box+...
				S0*Q_set0{1}*v_disturb_box+...
				(eye(4*n)+S0*Q_set0{1}*Cm0)*[eye(3);dyn.A;dyn.A^2;dyn.A^3]*xi_disturb_box;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Optimization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimization Variables
% ++++++++++++++++++++++
dyn.x0 	= sdpvar(n,1,'full');

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

	[S0,H0,Cm0,xi0m,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(dyn,T_i,'missing',find(L{pattern_ind} == 0)-1);

	positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

	%Select all influenced states
	sel_influenced_states = [];
	for i = 1 : T_i
		sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
	end

	noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ dyn.eta_w * ones(2*wd*T_i,1) ; dyn.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M2 * ones(2*n*T_i,1) - [eye(n*T_i);-eye(n*T_i)]*sel_influenced_states*H0*r{pattern_ind} ];
	noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ dyn.eta_w * ones(2*wd*T_i,1) ; dyn.eta_v * ones(2*p*T_i,1) ; M1 * ones(2*n,1) ] <= M3 * ones(2*n,1) - [eye(n);-eye(n)]*select_m(T_i,T_i)*H0*r{pattern_ind} ];

	%Dual relationship to design variables
	pre_xi = [];
	for i = 0:T_i
		pre_xi = [ pre_xi ; dyn.A^i];
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
		[S0,H0,Cm0,~,B_w_big,C_v_big] = create_skaf_n_boyd_matrices(dyn,T_i,'missing',find(L{pattern_ind} == 0)-1);

		Q_set{pattern_ind} = value(Q{pattern_ind});
		r_set{pattern_ind} = value(r{pattern_ind});
		F_set{pattern_ind} = value( (inv(value(eye(size(S0,2)) + Q{pattern_ind}*Cm0*H0)) ) * Q{pattern_ind});
		u0_set{pattern_ind} = value( inv(value(eye(size(S0,2)) + Q{pattern_ind}*Cm0*H0)) * r{pattern_ind} );
		
		%Fix up F and u0 to avoid NaN
		F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
		% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

		%Create Function Outputs
		opt_out.Q_set = Q_set;
		opt_out.r_set = r_set;

		contr = FHAE_pb(L,F_set,u0_set);

	end
end


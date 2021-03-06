function [results] = observer_comparison41(varargin)
	%observer_comparison41.m
	%	Finding minimal Zonotopes around the Lane Keeping System
	%	when the disturbance set is considered to be a zonotope

	%%%%%%%%%%%%%%%%%%%%%
	%% Experiment 41.1 %%
	%%%%%%%%%%%%%%%%%%%%%

	%For two zonotopes with dim=num_g identify the smallest enclosing zonotope.

	%% Constants

	%Define the in-body (the zonotope we aim to enclose)
	Z_in.G = [ [1;2;1]*(1/norm(5)) [1;1;1]*(1/norm(2)) [0;2;2]*(1/norm(8)) ];
    %Z_in.G = [ [1;2;1]*(1/norm(5)) [1;1;1]*(1/norm(3)) [0;2;2]*(1/norm(8)) [-1;0;1]*(1/norm(2)) ];
	Z_in.c = zeros(3,1);
	Z_in.dim = 3;
	Z_in.num_g = size(Z_in.G,2);

	Z_circ.G = [ [2;3]*(1/norm([2;3])) [1;-3]*(1/norm([1;-3]))];
	%Z_circ.G = [ [2;3]*(1/norm([2;3])) [1;-3]*(1/norm([1;-3])) [1;1]*norm([1;1]) ];
	Z_circ.c = zeros(2,1);
	Z_circ.dim = 2;
	Z_circ.num_g = size(Z_circ.G,2);

	to_poly = @(G,c) Polyhedron('A', [pinv(G);-pinv(G)] , 'b' , [pinv(G)*c + ones(size(G,2),1);-pinv(G)*c + ones(size(G,2),1)] );

	figure;
	subplot(1,2,1)
	plot(to_poly(Z_in.G,Z_in.c))
	title('''Inbody'' in Experiment 41.1')

	subplot(1,2,2)
	plot(to_poly(Z_circ.G,Z_circ.c))
	title('Circumbody in Experiment 41.1')

	figure;
	subplot(1,2,1)
	temp_poly = to_poly_v2(Z_in.G,Z_in.c);
	plot(temp_poly.slice(3,0))
	title('The slice of in-body that we will contain')

	subplot(1,2,2)
	plot([eye(2) zeros(2,1)]*temp_poly)
	title('Projection of in-body into 2d')

	%Define Projection
	A = [eye(Z_circ.dim) zeros(Z_circ.dim,1)];
	b = ones(Z_circ.dim,1); %zeros(Z_circ.dim,1);

	%Attempt to find the minimum volume zonotope that includes the projected cube

	%% Optimization

	theta_2 = sdpvar(1,1,'full');

	Nu1 = sdpvar(Z_in.dim, Z_circ.num_g, 'full')
	Nu2 = sdpvar(Z_circ.dim, Z_circ.num_g, 'full');
	Lambda1 = sdpvar( 2*Z_in.num_g , Z_circ.num_g , 'full' );

	constr_contain1 = [ Nu1'*Z_in.c - Nu2'*(b-Z_circ.c) + Lambda1'*ones(2*Z_in.num_g,1) <= theta_2*ones(Z_circ.num_g,1) ];
	constr_contain1 = constr_contain1 + [ [zeros(Z_circ.num_g,Z_in.dim+Z_in.num_g) eye(Z_circ.num_g)] == Nu1'*[eye(Z_in.dim), -Z_in.G, zeros(Z_in.dim,Z_circ.num_g)] + ...
																											Nu2'*[A zeros(Z_circ.dim,Z_in.num_g) -Z_circ.G] + ...
																											Lambda1'*[zeros(Z_in.num_g*2,Z_in.dim) [eye(Z_in.num_g); -eye(Z_in.num_g) ] zeros(Z_in.num_g*2,Z_circ.num_g) ] ];

	Nu3 = sdpvar(Z_in.dim, Z_circ.num_g, 'full');
	Nu4 = sdpvar(Z_circ.dim, Z_circ.num_g, 'full');
	Lambda2 = sdpvar( 2*Z_in.num_g , Z_circ.num_g , 'full' );

	constr_contain2 = [Nu3'*Z_in.c - Nu4'*(b-Z_circ.c) + Lambda2'*ones(2*Z_in.num_g,1) <= theta_2*ones(Z_circ.num_g,1)];
	constr_contain2 = constr_contain2 + [ - [zeros(Z_circ.num_g,Z_in.dim+Z_in.num_g) eye(Z_circ.num_g)] == Nu3'*[eye(Z_in.dim) -Z_in.G zeros(Z_in.dim,Z_circ.num_g)] + ...
																											Nu4'*[A zeros(Z_circ.dim,Z_in.num_g) -Z_circ.G] + ...
																											Lambda2'*[zeros(Z_in.num_g*2,Z_in.dim) [eye(Z_in.num_g); -eye(Z_in.num_g) ] zeros(Z_in.num_g*2,Z_circ.num_g) ] ];


	pos_constr = [Lambda1 >= 0] + [Lambda2 >= 0];
	%Call Optimizer

	ops = sdpsettings('verbose',1,'debug',1);
	optim1 = optimize(constr_contain1+constr_contain2+pos_constr, theta_2 , ops);

	%% Analyze Results
	figure;
	hold on;
	plot(to_poly_v2(Z_circ.G*value(theta_2),Z_circ.c),'Color','g')
	plot(A*temp_poly+b)

	res1 = value(theta_2);

	pause;

	clear Nu1 Nu2 Nu3 Nu4 Lambda1 Lambda2 

	%Attempt to duplicate the results using a simplified function.

	%Call optimizer
	CG = constr_gen();
	[dual_vars,constrs] = CG.get_zonot_inclusion_constr(Z_in,Z_circ,A,b,theta_2);
	optim2 = optimize( constrs , theta_2, ops);

	res2 = value(theta_2);

	pause;

	%%%%%%%%%%%%%%%%%%%%%
	%% Experiment 41.2 %%
	%%%%%%%%%%%%%%%%%%%%%

	%% Constants

	dt = 0.01;

	[~,lk_dyn] = get_lk_aff_dyn(dt);

	%Missing Data Parameters
	word_len = 6; window_len = 1;
	L = {};
	for L_ind = 2:1+window_len
		L{end+1} = ones(1,word_len);
		L{end}(L_ind) = 0;
	end
	L{end+1} = ones(1,word_len);

% 	L = {ones(1,word_len)};

	%Simulation constants
	num_runs = 100;

	%target sets
	con = load('data/system_examples/lk_inv.mat');
	lk_inv_set = con.LK_inv;
	lk_cis_zono.c = zeros(size(con.dyn_lk.A,1),1);
	temp = pinv(lk_inv_set.A);
	lk_cis_zono.G = temp(:,[1:2:end]);

	start_set = lk_inv_set*(1/8);

	disp('Converted LK Invariant Set into a Zonotope Representation.')

	%Noise Constants
	eta_w = 0.05; eta_v = 0.02;

	%Create Dynamics
	n_x = size(con.dyn_lk.A,1);
	C = [1,0,0,0;0,0,1,0];
	n_y = size(C,1);

	temp_sys1 = ss(con.dyn_lk.A,con.dyn_lk.B,C,0);
	temp_sys2 = ss(con.dyn_lk.A,[1;0;0;0],C,0);

	temp_dsys1 = c2d(temp_sys1,dt);
	temp_dsys2 = c2d(temp_sys2,dt);

	lk_dyn_disc = Aff_Dyn(	temp_dsys1.A,temp_dsys1.B,zeros(n_x,1),C,...
							eta_w,eta_v,...
							temp_dsys2.B, eye(n_y));

	% simple_LK = Aff_Dyn(	con.dyn_lk.A,con.dyn_lk.B,con.dyn_lk.F,eye(n_x),...
	% 						eta_w,eta_v, ...
	% 						con.dyn_lk.B, eye(n_y) );

	%System Parameters
	A = [ zeros(2,1) [1;-20] ];
	n_x = size(A,1);

	B = [ 0; 1];
	F = zeros(n_x,1);
	
	C = eye(2);%[1,0];
	n_y = size(C,1);

	temp_sys = ss(A,B,C,0);
	dt = 0.05;
	temp_dsys = c2d(temp_sys,dt);

	temp_sys2 = ss(A,[1;0],C,0);
	temp_dsys2 = c2d(temp_sys2,dt);

	eta_w = 0.05; eta_v = 0.1;

	simple_LK = Aff_Dyn(	temp_dsys.A,temp_dsys.B,F,C,...
							eta_w,eta_v, ...
							temp_dsys2.B, eye(n_y) );

	disp('Created Aff_Dyn object.')

	lane_width = 2*0.9;

	%Desired Constants
	M1 = 0.3;

	%%%%%%%%%%%%%%%
	%% Synthesis %%
	%%%%%%%%%%%%%%%

	ad = lk_dyn_disc;
	Z1 = lk_cis_zono;
	Z2 = lk_cis_zono;
	%Z_target = Polyhedron('ub',M1*ones(1,n_x),'lb',-M1*ones(1,n_x));
	target_set = Polyhedron('ub',M1*ones(1,size(ad.A,1)),'lb',-M1*ones(1,size(ad.A,1)));

	alpha1 = 1;

	%+++++++++++++++++++
	%Synthesis Constants
	n = size(ad.A,1);
	m = size(ad.B,2);
	p = size(ad.C,1);
	wd = size(ad.B_w,2);
	vd = size(ad.C_v,2);

	num_g1 = size(Z1.G,2);
	dim_Z1 = size(Z1.G,1);
	num_g2 = size(Z2.G,2);
	dim_Z2 = size(Z2.G,1);
	num_cis_ineqs = size(lk_inv_set.A,1);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%+++++++++++++++++++++++++++++
	%Create Optimization Variables
	mu2 = sdpvar(1,1,'full');
	ad.x0 = sdpvar(n,1,'full');

	max_T_i = -1;
	for pattern_ind = 1 : length(L)
		
	end

	disp('Optimization variables created.')

	%++++++++++++++++++
	%Create Constraints
	obj_constrs = [];
	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = []; 
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = [];
	ineq_constrs = []; eq_constrs = [];

	for pattern_ind = 1 : length(L)

        %% Defining Constants
        % +++++++++++++++++++
        
		%Define Pattern-Based Constants
		T_i = length(L{pattern_ind});

		%Define Zonotope Templates
		nu_polyt = Polyhedron(	'A', [	[eye(T_i*wd);-eye(T_i*wd)], zeros(2*T_i*wd,num_g1+T_i*vd) ;
										zeros(2*T_i*vd,T_i*wd), [eye(T_i*vd);-eye(T_i*vd)], zeros(2*T_i*vd,num_g1);
										zeros(2*dim_Z1,T_i*(wd+vd)) [pinv(Z1.G);-pinv(Z1.G)] ] , ...
								'b', [	ad.eta_w*ones(2*T_i*wd,1);
										ad.eta_v*ones(2*T_i*vd,1);
										alpha1*ones(dim_Z1,1) + pinv(Z1.G)*Z1.c;
										alpha1*ones(dim_Z1,1) - pinv(Z1.G)*Z1.c] );

		nu_zonot.G = [	ad.eta_w*eye(T_i*wd) zeros(T_i*wd,T_i*vd+num_g1) ;
						zeros(T_i*vd,T_i*wd) ad.eta_v*eye(T_i*vd) zeros(T_i*vd,num_g1) ;
						zeros(dim_Z1,T_i*(wd+vd)) Z1.G ];

		nu_zonot.c = [	zeros(T_i*(wd+vd),1); Z1.c ];

		%Target Polyt
		temp_gens = {};
		for i = 1:T_i
			temp_gens{i} = Z2.G;
		end
		tube_zonot.G = blkdiag(temp_gens{:});

		tube_zonot.c = repmat(Z2.c,T_i,1);

        % Creating Variables
        % ++++++++++++++++++
        
        T_i = length(L{pattern_ind});
		% Feedback Variables
		Q{pattern_ind} = sdpvar(m*T_i,p*T_i,'full');
		r{pattern_ind} = sdpvar(m*T_i,1,'full');

		% Dual Variables
		Pi_1{pattern_ind} = sdpvar(2*size(nu_zonot.G,1), size(nu_zonot.G,1)+size(tube_zonot.G,1),'full');
        Pi_2{pattern_ind} = sdpvar(2*size(nu_zonot.G,1), size(nu_zonot.G,1)+size(tube_zonot.G,1),'full');
		Pi_3{pattern_ind} = sdpvar(2*size(nu_zonot.G,1), size(nu_zonot.G,1)+size(Z2.G,1),'full');

		Lambda_1{pattern_ind} = sdpvar(size(nu_zonot.G,1)+size(tube_zonot.G,1), size(Pi_1{pattern_ind},2),'full');
        Lambda_2{pattern_ind} = sdpvar(size(nu_zonot.G,1)+size(tube_zonot.G,1), size(Pi_1{pattern_ind},2),'full');
		Lambda_3{pattern_ind} = sdpvar(dim_Z2,T_i*(wd+vd)+dim_Z1,'full');

		%Zonotope Variables
		theta1{pattern_ind} = sdpvar(num_g1*T_i,1,'full');
		theta2{pattern_ind} = sdpvar(num_g2*T_i,1,'full');

		%Find the maximum T_i
		if T_i > max_T_i
			max_T_i = T_i;
		end
        
		% Creating Constraints
		% ++++++++++++++++++++

        ad_missing = ad;
        ad_missing.C = zeros(size(ad.C));
        ad_missing.C_v = zeros(size(ad.C_v));
        
		[H0,S0,Cm0,J,f_bar,B_w_big,C_v_big] = get_mpc_matrices([ad,ad_missing],L{pattern_ind}+1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		%Select all influenced states
		sel_influenced_states = [];
		for i = 1 : T_i
			sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
		end

		%Create variables to make this look like experiment 39
		G{pattern_ind} = [ 	(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*B_w_big ...
							S0*Q{pattern_ind}*C_v_big ...
							(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*J ];

		clear A b n1 n2

		A = sel_influenced_states*G{pattern_ind};

		b = sel_influenced_states*((eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar + S0*r{pattern_ind});

		n1 = size(nu_zonot.G,2); dim_nu = size(nu_zonot.G,1);
		n2 = size(tube_zonot.G,2); dim_tube = size(tube_zonot.G,1);

		ineq_constrs = ineq_constrs + [ Pi_1{pattern_ind}'*ones(size(Pi_1{pattern_ind},1),1) + Lambda_1{pattern_ind}'* [nu_zonot.c; -b + tube_zonot.c ] <= mu2*ones(size(tube_zonot,2),1) ];
		ineq_constrs = ineq_constrs + [ Pi_2{pattern_ind}'*ones(size(Pi_2{pattern_ind},1),1) + Lambda_2{pattern_ind}'* [nu_zonot.c; -b + tube_zonot.c ] <= mu2*ones(size(tube_zonot,2),1) ];

		eq_constrs = eq_constrs + [ [zeros(n1,n*T_i), eye(n1), zeros(n1,num_g2*T_i); zeros(n1,n*T_i), -eye(n1), zeros(n1,num_g2*T_i)]'*Pi_1{pattern_ind} + [eye(size(A,2)),-nu_zonot.G,zeros(dim_nu,n2); A,zeros(size(A,1),n1),tube_zonot.G]'*Lambda_1{pattern_ind} == [zeros(n2,dim_tube),zeros(n2,n1),eye(n2)]'];

		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * nu_polyt.b <= mu2 * ones(2*n*T_i,1) - [eye(n*T_i);-eye(n*T_i)]*sel_influenced_states*H0*r{pattern_ind} ];
		noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * nu_polyt.b <= target_set.b - target_set.A*select_m(T_i,T_i)*H0*r{pattern_ind} ];

        b = select_m(T_i,T_i)*((eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar + S0*r{pattern_ind});
        
        ineq_constrs = ineq_constrs + [ Pi_3{pattern_ind}'*ones(size(Pi_3{pattern_ind},1),1) + Lambda_3{pattern_ind}'* [nu_zonot.c; -b + Z2.c ] <= ones(size(tube_zonot,2),1) ];
        
		%Dual relationship to design variables
		pre_xi = [];
		for i = 0:T_i
			pre_xi = [ pre_xi ; ad.A^i];
		end

		bounded_disturb_matrix = [ [ eye(wd*T_i) ; -eye(wd*T_i) ] zeros(2*wd*T_i,vd*T_i+n) ;
									zeros(2*vd*T_i,wd*T_i) [ eye(vd*T_i) ; -eye(vd*T_i) ] zeros(2*vd*T_i,n) ;
									zeros(num_cis_ineqs,(vd+wd)*T_i) start_set.A ];

		dual_equal_constrs = dual_equal_constrs + [Pi_1{pattern_ind} * nu_polyt.A == [eye(n*T_i); -eye(n*T_i)]*sel_influenced_states*G{pattern_ind} ];
		dual_equal_constrs = dual_equal_constrs + [Pi_2{pattern_ind} * nu_polyt.A == target_set.A*select_m(T_i,T_i)*G{pattern_ind}];

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

	disp('Constraints Created.')
	
	%Optimization
	obj_fcn = mu2;
	ops = sdpsettings('verbose',1,'debug',1);
	optim0 = optimize(positive_constr+noise_constrs+dual_equal_constrs+l_diag_constr+shared_Q_constrs+shared_r_constrs+obj_constrs, ...
						obj_fcn, ...
						ops);

	disp('Optimization Function Returned.')

    opt_out = optim0;
    if opt_out.problem ~= 0
        oc39_contr1 = [];
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
            F_set{pattern_ind} = value( (inv(value(eye(size(H0,2)) + Q{pattern_ind}*Cm0*H0)) ) * Q{pattern_ind});
            u0_set{pattern_ind} = value( inv(value(eye(size(H0,2)) + Q{pattern_ind}*Cm0*H0)) * r{pattern_ind} );

            %Fix up F and u0 to avoid NaN
            F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
            % u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;

        end

        %Create Function Outputs
        opt_out.Q_set = Q_set;
        opt_out.r_set = r_set;

        oc39_contr1 = FHAE_pb(L,F_set,u0_set);
    end
    
    M2 = value(mu2);
	%[ oc37_opt1 , oc37_contr1 ] = simple_LK.eq_rec_design_pb( 'Min_M2' , M1 , L );

	%%%%%%%%%%%%%
	%% Results %%
	%%%%%%%%%%%%%

	%Simulate
	[ ctrl_sim1, ~ ] = oc39_contr1.simulate_n_runs( ad , M1 , num_runs );

	ctrl_sim1_mod = [];
	for ind = 1:length(ctrl_sim1)
		ctrl_sim1_mod(:,:,ind) = ad.C*ctrl_sim1{ind};
	end

	figure;
	hold on;
	bar([0:word_len]*dt,[M1 ones(1,word_len-1)*M2 M1],1,'w')
	bar([0:word_len]*dt,-1*[M1 ones(1,word_len-1)*M2 M1],1,'w')
	plot([0:word_len]*dt,(lane_width/2)*ones(1,word_len+1),'k','LineWidth',2)
	plot([0:word_len]*dt,-(lane_width/2)*ones(1,word_len+1),'k','LineWidth',2)
	plot([0:word_len]*dt,zeros(1,word_len+1),'k--','LineWidth',2)


	for sim_num = 1:10
		plot([0:word_len]*dt,ctrl_sim1_mod(1,:,sim_num))
	end

	axis([0,word_len*dt,-M2,M2])

	% Saving
    optim0.M2= M2;
    
	results.exp1.sys = simple_LK;
	results.exp1.opt = optim0;
	results.exp1.contr = oc39_contr1;
end
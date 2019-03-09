function [results] = observer_comparison39(varargin)
	%observer_comparison39.m
	%	Finding minimal Zonotopes around the Lane Keeping System

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

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

	%+++++++++++++++++++
	%Synthesis Constants
	n = size(ad.A,1);
	m = size(ad.B,2);
	p = size(ad.C,1);
	wd = size(ad.B_w,2);
	vd = size(ad.C_v,2);

	num_g = size(lk_cis_zono.G,2);
	num_cis_ineqs = size(lk_inv_set.A,1);

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	%+++++++++++++++++++++++++++++
	%Create Optimization Variables
	mu2 = sdpvar(1,1,'full');
	theta1 = sdpvar(num_g,1,'full');
	theta2 = sdpvar(num_g,1,'full');
	ad.x0 = sdpvar(n,1,'full');

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

	disp('Optimization variables created.')

	%++++++++++++++++++
	%Create Constraints
	obj_constrs = [];
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

		noise_constrs = noise_constrs + [ Pi_1{pattern_ind} * [ ad.eta_w * ones(2*wd*T_i,1) ; ad.eta_v * ones(2*vd*T_i,1) ; start_set.b ] <= mu2 * ones(2*n*T_i,1) - [eye(n*T_i);-eye(n*T_i)]*sel_influenced_states*H0*r{pattern_ind} ];
		noise_constrs = noise_constrs + [ Pi_2{pattern_ind} * [ ad.eta_w * ones(2*wd*T_i,1) ; ad.eta_v * ones(2*vd*T_i,1) ; start_set.b ] <= M1*ones(2*n,1) - [eye(n);-eye(n)]*select_m(T_i,T_i)*H0*r{pattern_ind} ];

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
									zeros(num_cis_ineqs,(vd+wd)*T_i) start_set.A ];

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
	[ ctrl_sim1, ~ ] = oc39_contr1.simulate_n_runs( ad , M1 , num_runs )

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
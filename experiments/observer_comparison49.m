function [results] = observer_comparison49(varargin)
	%observer_comparison49.m
	%
	%Description:
	%	The objective of this is to prototype the use of the linear encoding for Zonotope containment from Sadraddini.
	%
	%Assumptions:

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dt = 0.1;
	[~,lk_ad] = get_lk_aff_dyn(dt);
	
	n = size(lk_ad.A,1);
	m = size(lk_ad.B,2);
	p = size(lk_ad.C,1);
	wd = size(lk_ad.B_w,2);
	vd = size(lk_ad.C_v,2);

	results.ad = lk_ad;

	L = {ones(1,6)};
	mu1 = 1;
	mu2 = 10;
	mu3 = 3;

	%Create Zonotope Shapes
	hypercube_bds = [0.5,0.5,0.05,0.05];
	%Z1 = Zonotope(diag(hypercube_bds),zeros(n,1)); P_M1 = Polyhedron('lb',-hypercube_bds,'ub',hypercube_bds);
	Z1 = Zonotope(eye(n),zeros(n,1)); P_M1 = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));

	Z2 = Z1;
	Z2.G = 10*Z1.G;

	G0 = randi(100,n);
	for col_ind = 1:size(G0,2)
		G0(:,col_ind) = (1/norm(G0(:,col_ind)))*G0(:,col_ind);
	end
	%Z0 = Zonotope(G0,zeros(n,1));
	%Z0 = Zonotope([roty(30),zeros(3,1);zeros(1,3),1],zeros(n,1));
	Z0 = Zonotope(diag(hypercube_bds),zeros(n,1));

	results.Z0 = Z0;
	results.Z1 = Z1;
	results.Z2 = Z2;

	% Create Constraint Generator
	cg = constr_gen();

	%Optimization options
	ops = sdpsettings('verbose',1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Testing Zonotope Inclusion Constraint %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	unit_cube = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));

	alpha0 = 0.001;
	[dual_vars,incl_constr] = cg.create_sadraddini_AH_inclusion_constr(	Z1.c,Z1.G,unit_cube.A,unit_cube.b,...
																		Z0.c,alpha0*Z0.G,unit_cube.A,unit_cube.b);

	%Verify if Z1 is included in Z0. (It should not be here)
	optim0 = optimize(incl_constr,[],ops);
	results.incl_tests.opt_out0 = optim0;

	alpha0 = 5;
	[dual_vars,incl_constr] = cg.create_sadraddini_AH_inclusion_constr(	Z1.c,Z1.G,unit_cube.A,unit_cube.b,...
																		Z0.c,alpha0*Z0.G,unit_cube.A,unit_cube.b);

	%Verify if Z1 is included in Z0. (It should be here)
	optim1 = optimize(incl_constr,[],ops);
	results.incl_tests.opt_out1 = optim1;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Equalized Recovery with Zonotopes %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear mu2
	mu2 = sdpvar(1,1,'full');

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	obj_fcn = mu2;

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

	ad_prime = lk_ad;
	ad_prime.C = zeros(size(lk_ad.C));

	ad_arr = [ad_prime,lk_ad];

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = []; obj_constrs = [];

	for pattern_ind = 1 : length(L)
		sigma_i = L{pattern_ind};
		T_i = length(sigma_i);
		% Creating Constraints
		% ++++++++++++++++++++

		[H0,S0,Cm0,J0,f_bar,B_w_big,C_v_big] = get_mpc_matrices(ad_arr,'word',L{pattern_ind}+1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		P_wT = 1; P_vT = 1;
		for t_idx = 1:T_i
			P_wT = P_wT*ad_arr(sigma_i(t_idx)).P_w;
			P_vT = P_vT*ad_arr(sigma_i(t_idx)).P_v;
		end

		bounded_disturb_matrix = [ 	P_wT.A zeros(size(P_wT.A,1),vd*T_i+n) ;
									zeros(size(P_vT.A,1),size(P_wT.A,2)) P_vT.A zeros(size(P_vT.A,1),n) ;
									zeros(size(P_M1.A,1),(vd+wd)*T_i) P_M1.A ];

		G{pattern_ind} = [ 	(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*B_w_big ...
							S0*Q{pattern_ind}*C_v_big ...
							(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*J0 ];

		[extra_dual_vars,z_incl_constr] = cg.create_sadraddini_AH_inclusion_constr(	S0*r{pattern_ind}+(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar, ...
																					G{pattern_ind} , ...
																					bounded_disturb_matrix,[ P_wT.b ; P_vT.b ; P_M1.b ], ...
																					kron(ones(T_i+1,1),Z0.c),kron(eye(T_i+1),Z0.G), ...
																					[eye((T_i+1)*n);-eye((T_i+1)*n)],mu2*ones(2*(T_i+1)*n,1));


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

	optim3 = optimize(positive_constr+z_incl_constr+l_diag_constr+shared_Q_constrs+shared_r_constrs+obj_constrs, ...
			obj_fcn, ...
			ops);

	optim3.M2 = value(mu2);

	results.lk_tests.new_method.opt_out = optim3;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Compare with result of the function when using hyperbox templates. %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[optim4,contr] = lk_ad.rec_synthesis( 'Equalized' , 'prefix' , 'Feasible Set' , P_M1 , optim3.M2*unit_cube , L );

	results.lk_tests.hyperb_method.opt_out = optim4;
	%The hyperbox problem with identical scaling is not feasible.

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Uze the Zonotope Method to Find a Minimal Hyperbox Solution %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear mu2
	mu2 = sdpvar(1,1,'full');

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

	obj_fcn = mu2;

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

	ad_prime = lk_ad;
	ad_prime.C = zeros(size(lk_ad.C));

	ad_arr = [ad_prime,lk_ad];

	shared_Q_constrs = []; shared_r_constrs = [];
	dual_equal_constrs = [];
	positive_constr = [];
	noise_constrs = [];
	l_diag_constr = []; obj_constrs = [];

	for pattern_ind = 1 : length(L)
		sigma_i = L{pattern_ind};
		T_i = length(sigma_i);
		% Creating Constraints
		% ++++++++++++++++++++

		[H0,S0,Cm0,J0,f_bar,B_w_big,C_v_big] = get_mpc_matrices(ad_arr,'word',L{pattern_ind}+1);

		positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

		P_wT = 1; P_vT = 1;
		for t_idx = 1:T_i
			P_wT = P_wT*ad_arr(sigma_i(t_idx)).P_w;
			P_vT = P_vT*ad_arr(sigma_i(t_idx)).P_v;
		end

		bounded_disturb_matrix = [ 	P_wT.A zeros(size(P_wT.A,1),vd*T_i+n) ;
									zeros(size(P_vT.A,1),size(P_wT.A,2)) P_vT.A zeros(size(P_vT.A,1),n) ;
									zeros(size(P_M1.A,1),(vd+wd)*T_i) P_M1.A ];

		G{pattern_ind} = [ 	(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*B_w_big ...
							S0*Q{pattern_ind}*C_v_big ...
							(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*J0 ];

		[extra_dual_vars,z_incl_constr] = cg.create_sadraddini_AH_inclusion_constr(	S0*r{pattern_ind}+(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar, ...
																					G{pattern_ind} , ...
																					bounded_disturb_matrix,[ P_wT.b ; P_vT.b ; P_M1.b ], ...
																					kron(ones(T_i+1,1),Z1.c),kron(eye(T_i+1),Z1.G), ...
																					[eye((T_i+1)*n);-eye((T_i+1)*n)],mu2*ones(2*(T_i+1)*n,1));


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

	optim5 = optimize(positive_constr+z_incl_constr+l_diag_constr+shared_Q_constrs+shared_r_constrs+obj_constrs, ...
			obj_fcn, ...
			ops);

	optim5.M2 = value(mu2);

	results.lk_tests.box_template.opt_out = optim5;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Plot all 3 Polyhedron %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Z3 = Z0;
	Z3.G = optim3.M2 * Z3.G;

	Z5 = Z1;
	Z5.G = optim5.M2 * Z5.G;

	A = [1,0,0,0;0,0,1,0];

	figure;
	hold on;
	plot(A*Z5.to_poly())
	plot(A*Z3.to_poly(),'Color','magenta')

	% legend('Hypercube Template','Zonotope Template')
	saveas(gcf,'results/nahs2019/4d_lk_minM2_comparison.png')

end
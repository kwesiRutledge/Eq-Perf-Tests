function [opt_out, contr] = rec_synthesis( varargin )
%Description:
%	Searches for a set of feedback gains which can satisfy Recovery specifications (Equalized or Free)
%	
%Usage:
%	[opt_out, contr] = rec_synthesis( ad , rec_type , prefix_flag, prob_type )
%	[opt_out, contr] = rec_synthesis( ad , 'Equalized' , 'prefix' , 'Feasible Set' , M1, M2, L )
%	[opt_out, contr] = rec_synthesis( ad , 'Equalized' , 'prefix' , 'Feasible Set' , P1, P2, L )
%	[opt_out, contr] = rec_synthesis( ad , 'Equalized' , 'prefix' , 'Feasible Set' , P1, P2, L , 'Pu' , Pu)
%	[opt_out, contr] = rec_synthesis( ad , 'Equalized' , 'prefix' , 'Minimize M2' , M1, L )
%	[opt_out, contr] = rec_synthesis( ad , 'Equalized' , 'prefix' , 'Minimize M2' , M1, L , 'Pu', Pu)
%	[opt_out, contr] = rec_synthesis( ad , 'Free' , 'prefix' , 'Feasible Set' , P1, P2, P3, L )
%	[opt_out, contr] = rec_synthesis( ad , 'Free' , 'time' , 'Feasible_set' , P1, P2, P3, L)
%	[opt_out, contr] = rec_synthesis( ad , 'Free' , 'prefix' , 'Min_M2' , P1, P3, L)
%
%Inputs:
%	ad 			- Affine Dynamics
%	rec_type 	- A string indicating the type of recovery problem that we want to solve.
%				  Options: 'Equalized' or 'Free'
%	prefix_flag - A string indicating whether or not the controller/estimator designed will be prefix-based or not.
%				  Options: 'prefix' or 'time'
%	prob_type 	- A string indicating the structure of the problem statement.
%				  Options: 'Feasible Set', 'Minimize M2', 'Minimize Z2', 'Minimize M3'
%				  Not all options will be available depending on the choice of Recovery Type.
%				  Example: 'Feasible Set' for an equalized recovery problem, indicates that the inputs to this function
%							will be a set of parameters that define the equalized recovery problem (M1, M2, L) and the 
%							algorithm finds a controller/estimator that satisfies the 3 parameters, or else it returns
%							an infeasible flag.
%				  Example: 'Minimize M3' for a free recovery problem, indicates that the inputs to this function will be
%							all parameters that define the free recovery problem EXCEPT for M3. (i.e. the parameters are
%							M1,M2,M3, and L but the user only needs to provide M1,M2, and L). The algorithm finds a
%							controller/estimator that satisfies the 3 parameters AND minimizes the value of M3 if possible,
%							or else it returns an infeasible flag.
%	M1 			- A scalar value for M1.
%	M2 			- A scalar vakye fir M2.


	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	ad = varargin{1};
	rec_type = varargin{2};
	prefix_flag = varargin{3};
	prob_type = varargin{4};

	if ~(strcmp(rec_type,'Free') || strcmp(rec_type,'Equalized') )
		error('Unrecognized recovery type given.')
	end

	if strcmp(rec_type,'Equalized') && strcmp(prob_type,'Feasible Set')
		M1 = varargin{5};
		M2 = varargin{6};
		L = varargin{7};
		curr_idx = 8;
	elseif strcmp(rec_type,'Equalized') && strcmp(prob_type,'Minimize M2')
		M1 = varargin{5};
		L = varargin{6};
		curr_idx = 7;
	elseif strcmp(rec_type,'Equalized') && strcmp(prob_type,'Minimize Z2')
		M1 = varargin{5};
		Z1 = varargin{6};
		Z2 = varargin{7};
		L = varargin{8};
		curr_idx = 9;
	elseif strcmp(rec_type,'Free') && strcmp(prob_type,'Feasible Set')
		M1 = varargin{5};
		M2 = varargin{6};
		M3 = varargin{7};
		L = varargin{8};
		curr_idx = 9;
	elseif strcmp(rec_type,'Free') && strcmp(prob_type,'Minimize M2')
		M1 = varargin{5};
		M3 = varargin{6};
		L = varargin{7};
		curr_idx = 8;
	elseif strcmp(rec_type,'Free') && strcmp(prob_type,'Minimize M3')
		M1 = varargin{5};
		M2 = varargin{6};
		L = varargin{7};
		curr_idx = 8;
	end

	if curr_idx < nargin
		while curr_idx < nargin
			switch varargin{curr_idx}
				case 'Pu'
					Pu = varargin{curr_idx+1};
					curr_idx = curr_idx + 2;
				case 'System Type'
					sys_type = varargin{curr_idx+1};
					curr_idx = curr_idx + 2;
				case 'debug'
					debug_flag = varargin{curr_idx+1};
					curr_idx = curr_idx + 2;
				otherwise
					error('Unexpected input.')
			end
		end
	end

	if ~exist('sys_type')
		sys_type = 'Missing Data';
	end

	if ~exist('debug_flag')
		debug_flag = 1;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	switch sys_type
		case 'Missing Data'
			n = size(ad.A,1);
			m = size(ad.B,2);
			p = size(ad.C,1);
			wd = size(ad.B_w,2);
			vd = size(ad.C_v,2);

			%Create 'Switching' System for Missing Data and Not Modes
			ad_prime = ad;
			ad_prime.C = zeros(size(ad.C));
			ad_prime.C_v = zeros(size(ad.C_v));
			ad_arr = [ad_prime,ad];

			%Update the language so that it does not use 'zeros'
			for sig_idx = 1:length(L)
				L{sig_idx} = L{sig_idx} + 1;
			end

		case 'Switched'
			%If the string tells you that this is a switched system then the input 'ad' should
			%be an array.
			ad_arr = ad;

			%Check the number of modes in ad_arr and verify that there are no elements in L
			ad_len = length(ad);
			for sig_idx = 1:length(L)
				sigma_i = L{sig_idx};
				if ~all(sigma_i <= ad_len)
					error('There are more modes referenced in the language than there are provided in the Aff_Dyn array.')
				end
			end

			%Define Sizing Parameters
			n = size(ad(1).A,1);
			m = size(ad(1).B,2);
			p = size(ad(1).C,1);
			wd = size(ad(1).B_w,2);
			vd = size(ad(1).C_v,2);			

		otherwise
			error('Unexpected System Type. Expecting ''Missing Data'' or ''Switched''.')
	end

	%Select matrix
	select_m = @(t,T_r) [zeros(n,t*n), eye(n), zeros(n,(T_r-t)*n) ];

	unit_box = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));

	ops = sdpsettings('verbose',debug_flag);

	cg = constr_gen();

	%%%%%%%%%%%%%%%%%%%
	%% Optimizations %%
	%%%%%%%%%%%%%%%%%%%

	positive_constr = []; dual_constrs = [];
	l_diag_constr = []; prefix_constrs = [];
	obj_fcn = [];

	switch rec_type
	case 'Equalized'
		%Create Objective
		mu2 = sdpvar(1,1,'full');

		switch prob_type
		case 'Feasible Set'
			obj_fcn = [];
		case 'Minimize M2'
			obj_fcn = mu2;
		end

		if strcmp(prefix_flag,'prefix')

			%Creating Optimization Variables
			max_T_i = -1;
			for pattern_ind = 1 : length(L)
				T_i = length(L{pattern_ind});
				% Feedback Variables
				Q{pattern_ind} = sdpvar(m*T_i,p*T_i,'full');
				r{pattern_ind} = sdpvar(m*T_i,1,'full');

				% Dual Variables
				Pi_1{pattern_ind} = sdpvar(2*n*T_i,2*(wd+vd)*T_i+2*n,'full');
				Pi_2{pattern_ind} = sdpvar(2*n,2*(wd+vd)*T_i+2*n,'full');
				if exist('Pu')
					switch sys_type
						case 'Missing Data'
							if isa(M1,'Polyhedron')
								Pi_3{pattern_ind} = sdpvar(T_i*size(Pu.A,1),T_i*(size(ad.P_w.A,1)+size(ad.P_v.A,1))+size(M1.A,1),'full');
							else
								Pi_3{pattern_ind} = sdpvar(T_i*size(Pu.A,1),T_i*(size(ad.P_w.A,1)+size(ad.P_v.A,1))+2*n,'full');
							end
						case 'Switched'
							%Collect an array of the P_w
							P_w_dims = []; P_v_dims = [];
							for sys_idx = 1:length(ad)
								P_w_dims = [P_w_dims; size(ad(sys_idx).P_w.A)];
								P_v_dims = [P_v_dims; size(ad(sys_idx).P_v.A)];
							end

							Pi_3_dims = [0,0];
							for symbol_idx = 1:length(L{pattern_ind})
								Pi_3_dims(2) = Pi_3_dims(2) + P_v_dims(L{pattern_ind}(symbol_idx));
								Pi_3_dims(2) = Pi_3_dims(2) + P_w_dims(L{pattern_ind}(symbol_idx));
							end

							if isa(M1,'Polyhedron')
								Pi_3{pattern_ind} = sdpvar(T_i*size(Pu.A,1),Pi_3_dims(2)+size(M1.A,1),'full');
							else
								Pi_3{pattern_ind} = sdpvar(T_i*size(Pu.A,1),Pi_3_dims(2)+2*n,'full');
							end
					end
				end

				%Find the maximum T_i
				if T_i > max_T_i
					max_T_i = T_i;
				end
			end
		
			for pattern_ind = 1 : length(L)

				positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

				switch prob_type
				case 'Feasible Set'
					dual_constrs = dual_constrs + cg.get_er_constr(ad_arr,L{pattern_ind}, ...
																	Pi_1{pattern_ind},Pi_2{pattern_ind}, ...
																	Q{pattern_ind},r{pattern_ind}, ...
																	'Feasible Set', M1 , M2 );
				case 'Minimize M2'
					dual_constrs = dual_constrs + cg.get_er_constr(ad_arr,L{pattern_ind}, ...
																	Pi_1{pattern_ind},Pi_2{pattern_ind}, ...
																	Q{pattern_ind},r{pattern_ind}, ...
																	'Minimize M2', M1 , mu2 );
				case 'Minimize Z2'
					dual_constrs = dual_constrs + cg.get_er_constr(ad_arr,L{pattern_ind}, ...
																	Pi_1{pattern_ind},Pi_2{pattern_ind}, ...
																	Q{pattern_ind},r{pattern_ind}, ...
																	'Minimize Z2', M1 , Z1, Z2, mu2 );
				otherwise
					error('Unrecognized problem type. Cannot create dual constraints.')
				end

				if exist('Pu')
					%If there exist limits on the input that we can give, then introduce that constraint as well.
					dual_constrs = dual_constrs + cg.get_input_constr(ad_arr,L{pattern_ind}, ...
																		Pi_3{pattern_ind}, Q{pattern_ind}, r{pattern_ind}, ...
																		M1,Pu);
					positive_constr = positive_constr + [Pi_3{pattern_ind} >= 0];
				end

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
					prefix_constrs = prefix_constrs +  [Q{pattern_ind}( [1:ind_identical*m] , [1:ind_identical*p] ) == Q{patt_i}( [1:ind_identical*m] , [1:ind_identical*p] )];
					prefix_constrs = prefix_constrs +  [r{pattern_ind}( [1:ind_identical*m] ) == r{patt_i}([1:ind_identical*m]) ];
				end

			end

		elseif strcmp(prefix_flag,'time')

			warning('This flag combination is untested!')

			%Create Objective
			mu2 = sdpvar(1,1,'full');

			switch prob_type
			case 'Feasible Set'
				obj_fcn = [];
			case 'Minimize M2'
				obj_fcn = mu2;
			end

			T = length(L{1});
			L_star = 2*ones(1,T);
		    for sig_i = 1:length(L)
				L_star = bitand(L_star-1,L{sig_i}-1);
				L_star = L_star+1;
			end

			Pi_1 = sdpvar(2*n*T,2*(wd+vd)*T+2*n,'full');
			Pi_2 = sdpvar(2*n,2*(wd+vd)*T+2*n,'full');
			if exist('Pu')
				if isa(M1,'Polyhedron')
					Pi_3 = sdpvar(T*size(Pu.A,1),T*(size(ad.P_w.A,1)+size(ad.P_v.A,1))+size(M1.A,1),'full');
				else
					Pi_3 = sdpvar(T*size(Pu.A,1),T*(size(ad.P_w.A,1)+size(ad.P_v.A,1))+2*n,'full');
				end
			end

			% Feedback Variables
			Q = sdpvar(m*T,p*T,'full');
			r = sdpvar(m*T,1,'full');

			obj_constrs = [];
			obj_fcn = [];

			dual_equal_constrs = [];
			positive_constr = [];
			noise_constrs = [];
			l_diag_constr = [];

			
			% Creating Constraints
			% ++++++++++++++++++++

			positive_constr = positive_constr + [ Pi_1 >= 0, Pi_2 >= 0 ];

			switch prob_type
				case 'Feasible Set'
					dual_constrs = dual_constrs + cg.get_er_constr(ad_arr,L_star, ...
																	Pi_1,Pi_2, Q, r, ...
																	'Feasible Set', M1 , M2 );
				case 'Minimize M2'
					dual_constrs = dual_constrs + cg.get_er_constr(ad_arr,L_star, ...
																	Pi_1,Pi_2, Q,r, ...
																	'Minimize M2', M1 , mu2 );
				otherwise
					error('Unrecognized problem type. Cannot create dual constraints.')
				end

				if exist('Pu')
					%If there exist limits on the input that we can give, then introduce that constraint as well.
					dual_constrs = dual_constrs + cg.get_input_constr(ad_arr,L_star, ...
																		Pi_3, Q, r, ...
																		M1,Pu);
					positive_constr = positive_constr + [Pi_3 >= 0];
				end

			%Lower Diagonal Constraint
			for bl_row_num = 1 : T-1
				l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*m+1:bl_row_num*m], [bl_row_num*p+1:end] ) == 0 ];
			end
		else
			error('Unrecognized value for prefix_flag. Can only choose: ''prefix'' or ''time''.')
		end
	case 'Free'
		%Create Objective
		mu2 = sdpvar(1,1,'full');
		mu3 = sdpvar(1,1,'full');

		switch prob_type
		case 'Feasible Set'
			obj_fcn = [];
		case {'Minimize M2','Minimize Z2'}
			obj_fcn = mu2;
		case 'Minimize M3'
			obj_fcn = mu3;
		end

		if strcmp(prefix_flag,'prefix')

			%Creating Optimization Variables
			max_T_i = -1;
			for pattern_ind = 1 : length(L)
				T_i = length(L{pattern_ind});
				% Feedback Variables
				Q{pattern_ind} = sdpvar(m*T_i,p*T_i,'full');
				r{pattern_ind} = sdpvar(m*T_i,1,'full');

				% Dual Variables
				Pi_1{pattern_ind} = sdpvar(2*n*T_i,2*(wd+vd)*T_i+2*n,'full');
				Pi_2{pattern_ind} = sdpvar(2*n,2*(wd+vd)*T_i+2*n,'full');
				if exist('Pu')
					if isa(M1,'Polyhedron')
						Pi_3{pattern_ind} = sdpvar(T_i*size(Pu.A,1),T_i*(size(ad.P_w.A,1)+size(ad.P_v.A,1))+size(M1.A,1),'full');
					else
						Pi_3{pattern_ind} = sdpvar(T_i*size(Pu.A,1),T_i*(size(ad.P_w.A,1)+size(ad.P_v.A,1))+2*n,'full');
					end
				end

				%Find the maximum T_i
				if T_i > max_T_i
					max_T_i = T_i;
				end
			end
		
			for pattern_ind = 1 : length(L)

				positive_constr = positive_constr + [ Pi_1{pattern_ind} >= 0, Pi_2{pattern_ind} >= 0 ];

				switch prob_type
				case 'Feasible Set'
					dual_constrs = dual_constrs + cg.get_fr_constr(ad_arr,L{pattern_ind}, ...
																	Pi_1{pattern_ind},Pi_2{pattern_ind}, ...
																	Q{pattern_ind},r{pattern_ind}, ...
																	'Feasible Set', M1 , M2 , M3 );
				case 'Minimize M2'
					dual_constrs = dual_constrs + cg.get_fr_constr(ad_arr,L{pattern_ind}, ...
																	Pi_1{pattern_ind},Pi_2{pattern_ind}, ...
																	Q{pattern_ind},r{pattern_ind}, ...
																	'Minimize M2', M1 , mu2 , M3 );
				case 'Minimize M3'
					dual_constrs = dual_constrs + cg.get_fr_constr(ad_arr,L{pattern_ind}, ...
																	Pi_1{pattern_ind},Pi_2{pattern_ind}, ...
																	Q{pattern_ind},r{pattern_ind}, ...
																	'Minimize M3', M1 , M2 , mu3 );
				otherwise
					error('Unrecognized problem type. Cannot create dual constraints.')
				end

				if exist('Pu')
					%If there exist limits on the input that we can give, then introduce that constraint as well.
					dual_constrs = dual_constrs + cg.get_input_constr(ad_arr,L{pattern_ind}, ...
																		Pi_3{pattern_ind}, Q{pattern_ind}, r{pattern_ind}, ...
																		M1,Pu);
					positive_constr = positive_constr + [Pi_3{pattern_ind} >= 0];
				end

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
					prefix_constrs = prefix_constrs +  [Q{pattern_ind}( [1:ind_identical*m] , [1:ind_identical*p] ) == Q{patt_i}( [1:ind_identical*m] , [1:ind_identical*p] )];
					prefix_constrs = prefix_constrs +  [r{pattern_ind}( [1:ind_identical*m] ) == r{patt_i}([1:ind_identical*m]) ];
				end
			end

		elseif strcmp(prefix_flag,'time')

			warning('This flag combination is untested!')

			%Create Objective
			mu2 = sdpvar(1,1,'full');
			mu3 = sdpvar(1,1,'full');

			switch prob_type
			case 'Feasible Set'
				obj_fcn = [];
			case 'Minimize M2'
				obj_fcn = mu2;
			case 'Minimize M3'
				obj_fcn = mu3;
			end

			T = length(L{1});
			L_star = ones(1,T);
		    for sig_i = 1:length(L)
				L_star = bitand(L_star,L{sig_i});
			end

			Pi_1 = sdpvar(2*n*T,2*(wd+vd)*T+2*n,'full');
			Pi_2 = sdpvar(2*n,2*(wd+vd)*T+2*n,'full');
			if exist('Pu')
				if isa(M1,'Polyhedron')
					Pi_3 = sdpvar(T*size(Pu.A,1),T*(size(ad.P_w.A,1)+size(ad.P_v.A,1))+size(M1.A,1),'full');
				else
					Pi_3 = sdpvar(T*size(Pu.A,1),T*(size(ad.P_w.A,1)+size(ad.P_v.A,1))+2*n,'full');
				end
			end

			% Feedback Variables
			Q = sdpvar(m*T,p*T,'full');
			r = sdpvar(m*T,1,'full');
			
			% Creating Constraints
			% ++++++++++++++++++++

			positive_constr = positive_constr + [ Pi_1 >= 0, Pi_2 >= 0 ];

			switch prob_type
				case 'Feasible Set'
					dual_constrs = dual_constrs + cg.get_fr_constr(ad_arr,L_star, ...
																	Pi_1,Pi_2, Q, r, ...
																	'Feasible Set', M1 , M2 , M3 );
				case 'Minimize M2'
					dual_constrs = dual_constrs + cg.get_fr_constr(ad_arr,L_star, ...
																	Pi_1,Pi_2, Q,r, ...
																	'Minimize M2', M1 , mu2 , M3 );

				case 'Minimize M3'
					dual_constrs = dual_constrs + cg.get_fr_constr(ad_arr,L_star, ...
																	Pi_1,Pi_2, Q,r, ...
																	'Minimize M3', M1 , M2 , mu3  );

				otherwise
					error('Unrecognized problem type. Cannot create dual constraints.')
				end

				if exist('Pu')
					%If there exist limits on the input that we can give, then introduce that constraint as well.
					dual_constrs = dual_constrs + cg.get_input_constr(ad_arr,L_star, ...
																		Pi_3, Q, r, ...
																		M1,Pu);
					positive_constr = positive_constr + [Pi_3 >= 0];
				end

			%Lower Diagonal Constraint
			for bl_row_num = 1 : T-1
				l_diag_constr = l_diag_constr + [ Q(	[(bl_row_num-1)*m+1:bl_row_num*m], [bl_row_num*p+1:end] ) == 0 ];
			end
		end

	otherwise
		error('We currently do not support any additional problem types.')
	end

	constrs = positive_constr  + dual_constrs + l_diag_constr + prefix_constrs ;

	%% Call Optimizer %%
	opt_out = optimize(constrs,obj_fcn,ops);

	%%%%%%%%%%%%%%%%%%
	%% Make Outputs %%
	%%%%%%%%%%%%%%%%%%

	if strcmp(prefix_flag,'prefix')
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
				[H0,S0,Cm0,J0,f_bar,B_w_big,~] = get_mpc_matrices(ad_arr,'word',L{pattern_ind});

				Q_set{pattern_ind} = value(Q{pattern_ind});
				r_set{pattern_ind} = value(r{pattern_ind});
				F_set{pattern_ind} = value( (pinv(value(eye(size(Q{pattern_ind},1)) + Q{pattern_ind}*Cm0*S0)) ) * Q{pattern_ind});
				u0_set{pattern_ind} = value( pinv(value(eye(size(Q{pattern_ind},1)) + Q{pattern_ind}*Cm0*S0)) * r{pattern_ind} );
				
				%Fix up F and u0 to avoid NaN
				F_set{pattern_ind}( isnan(F_set{pattern_ind}) ) = 0;
				% u0_set{pattern_ind}( isnan(u0_set{pattern_ind}) ) = 0;
			end

			%Create Function Outputs
			opt_out.Q_set = Q_set;
			opt_out.r_set = r_set;

			for sig_idx = 1:length(L)
				L{sig_idx} = L{sig_idx} - 1;
			end

			contr = FHAE_pb(L,F_set,u0_set);

		end
	elseif strcmp(prefix_flag,'time')
		[H0,S0,Cm0,J0,f_bar,B_w_big,~] = get_mpc_matrices(ad_arr,'word',L_star);

		opt_out.Q = value(Q);
		opt_out.Q(isnan(value(Q))) = 0;
		opt_out.r = value(r);
		contr.F   = value( (pinv(value(eye(size(Q,1)) + opt_out.Q*Cm0*H0)) ) * opt_out.Q);
		contr.u0  = value( pinv(value(eye(size(Q,1)) + opt_out.Q*Cm0*H0)) * r );
	else
		error('Unrecognized value for prefix_flag. Can only choose: ''prefix'' or ''time''.')
	end

	switch prob_type
	case 'Minimize M2'
		opt_out.M2 = value(mu2);
	case 'Minimize M3'
		opt_out.M3 = value(mu3);
	end

end
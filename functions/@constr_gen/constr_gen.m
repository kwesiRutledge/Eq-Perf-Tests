classdef constr_gen
	%Description:
	%	This class abstracts away how we make constraints for the equalized recovery problems.
	properties

	end

	methods
		%Constructor
		function [generator] = constr_gen()
			%Description:
			%	The constructor for this class does nothing.

			disp('constr_gen() does not initialize any member variables.')
			% generator = [];

		end

		function [ dual_vars , constraint ] =  get_zonot_inclusion_constr(obj,Z_in,Z_circ,A,b,theta_2)
			%get_zonot_inclusion_constr.m
			%Description:
			%	Creates a zonotope inclusion constraint based on the desired property:
			%		A*Z_in + b \subseteq Z_circ
			%	This uses the interpretation of inclusion that our group developed in the Jupyter Notebook XXX.ipynb.
			%
			%Inputs:
			%	Z_in

			%Create Constraint for Upper Inclusion (theta_2 respects upper bound)
			Nu1 = sdpvar(Z_in.dim, Z_circ.num_g, 'full');
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

			%Results
			dual_vars = {Nu1,Nu2,Nu3,Nu4,Lambda1,Lambda2};
			constraint = pos_constr + constr_contain1 + constr_contain2;

		end

		function [ constraints ] = get_er_constr(varargin)
			%get_er_constr()
			%Description:
			%	Uses the appropriate dual variables (Pi1 and Pi2) along with the controller parameters
			%	(Q and r) to create the constraints corresponding to Equalized Recovery.
			%
			%Usage:
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,in_str,P_M1,P_M2)
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Feasible Set',M1,M2)
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Feasible Set',P_M1,P_M2)
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Minimize M2',M1,mu2)
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Minimize Z2',P_M1,Z1,Z2,mu2)
			%
			%Inputs:
			%	Q:	The "Q" variable in Q-Parameterization. Together with r, this defines a set of feedbacks
			%		(F,u_0) which can be determined by a nonlinear change of variables.
			%	r:	The "r" variable in Q-Parameterization. See note on Q.

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			cg = varargin{1};
			ad_arr = varargin{2};
			sigma_i = varargin{3};
			Pi1 = varargin{4};
			Pi2 = varargin{5};
			Q = varargin{6};
			r = varargin{7};
			in_str = varargin{8};

			switch in_str
			case 'Feasible Set'
				M1 = varargin{9};
				M2 = varargin{10};
			case 'Minimize M2'
				M1 = varargin{9};
				mu2 = varargin{10};
			case 'Minimize Z2'
				P_M1 = varargin{9};
				Z1 = varargin{10};
				Z2 = varargin{11};
				mu2 = varargin{12};
			otherwise
				error('Unrecognized input string/type of equalized recovery problem.')
			end

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			n = size(ad_arr(1).A,1);
			vd = size(ad_arr(1).C_v,2);
			wd = size(ad_arr(1).B_w,2);

			unit_box = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));

			if strcmp(in_str,'Feasible Set')
				if isa(M1,'Polyhedron') && isa(M2,'Polyhedron')
					P_M1 = M1;
					P_M2 = M2;
	            elseif isscalar(M1) && isscalar(M2)
					P_M1 = M1 * unit_box;
					P_M2 = M2 * unit_box;
				else
					error('Unrecognized input types for M1 and M2.')
				end
			elseif strcmp(in_str,'Minimize M2')
				if isa(M1,'Polyhedron')
					P_M1 = M1;
					P_M2 = unit_box;
				elseif isscalar(M1)
					P_M1 = M1 * unit_box;
					P_M2 = unit_box;
				end
			elseif strcmp(in_str,'Minimize Z2')
				if ~(isa(Z1,'Zonotope') && isa(Z2,'Zonotope'))
					error('Need both inputs to be zonotopes.')
				end
			end
				
			T_i = length(sigma_i);

			%Select matrix
			select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

			sel_influenced_states = [];
			prod_M2 = 1;
			for i = 1 : T_i
				%Select all influenced states
				sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
				%Create product for M2
				prod_M2 = prod_M2*P_M2;
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Creating Constraints %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%

			constraints = [];

			%Get Special Matrices
			[H0,S0,Cm0,J0,f_bar,B_w_big,C_v_big] = get_mpc_matrices(ad_arr,'word',sigma_i);
			P_wT = 1; P_vT = 1;
			for t_idx = 1:T_i
				P_wT = P_wT*ad_arr(sigma_i(t_idx)).P_w;
				P_vT = P_vT*ad_arr(sigma_i(t_idx)).P_v;
			end

			G = [ 	(eye(n*(T_i+1))+S0*Q*Cm0)*H0*B_w_big ...
					S0*Q*C_v_big ...
					(eye(n*(T_i+1))+S0*Q*Cm0)*J0 ];

			bounded_disturb_matrix = [ 	P_wT.A, zeros(size(P_wT.A,1),vd*T_i+n) ;
										zeros(size(P_vT.A,1),size(P_wT.A,2)), P_vT.A, zeros(size(P_vT.A,1),n) ;
										zeros(size(P_M1.A,1),(vd+wd)*T_i) P_M1.A ];
			switch in_str
			case 'Feasible Set'
				constraints = constraints + [ Pi1 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= prod_M2.b - prod_M2.A*sel_influenced_states*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];
				constraints = constraints + [ Pi2 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= P_M1.b - P_M1.A*select_m(T_i,T_i)*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];

				constraints = constraints + [Pi1 * bounded_disturb_matrix == prod_M2.A*sel_influenced_states*G ];
				constraints = constraints + [Pi2 * bounded_disturb_matrix == P_M1.A*select_m(T_i,T_i)*G];
			case 'Minimize M2'
				constraints = constraints + [ Pi1 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= mu2*prod_M2.b - prod_M2.A*sel_influenced_states*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];
				constraints = constraints + [ Pi2 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= P_M1.b - P_M1.A*select_m(T_i,T_i)*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];

				constraints = constraints + [Pi1 * bounded_disturb_matrix == prod_M2.A*sel_influenced_states*G ];
				constraints = constraints + [Pi2 * bounded_disturb_matrix == P_M1.A*select_m(T_i,T_i)*G];
			case 'Minimize Z2'
				[~,temp_constrs] = cg.create_sadraddini_AH_inclusion_constr(	S0*r{pattern_ind}+(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar, ...
																				G{pattern_ind} , ...
																				bounded_disturb_matrix,[ P_wT.b ; P_vT.b ; P_M1.b ], ...
																				kron(ones(T_i+1,1),Z2.c),kron(eye(T_i+1),Z2.G), ...
																				[eye((T_i+1)*n);-eye((T_i+1)*n)],mu2*ones(2*(T_i+1)*n,1));
				constraints = constraints + temp_constrs;

				[~,temp_constrs] = cg.create_sadraddini_AH_inclusion_constr(	select_m(T_i,T_i)*(S0*r{pattern_ind}+(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar), ...
																				select_m(T_i,T_i)*G{pattern_ind} , ...
																				bounded_disturb_matrix,[ P_wT.b ; P_vT.b ; P_M1.b ], ...
																				Z3.c,Z3.G, ...
																				[eye(n);-eye(n)],ones(2*n,1) );
				constraints = constraints + temp_constrs;
			otherwise
				error('Unknown problem type.')
			end

		end

		function [ constraints ] = get_fr_constr(varargin)
			%get_er_constr()
			%Description:
			%	Uses the appropriate dual variables (Pi1 and Pi2) along with the controller parameters
			%	(Q and r) to create the constraints corresponding to Equalized Recovery.
			%
			%Usage:
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,in_str,M1,M2,M3)
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Feasible Set',M1,M2,M3)
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Feasible Set',P_M1,P_M2,P_M3)
			%	constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Minimize M2',M1,mu2,M3)
            %   constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Minimize M3',M1,M2,mu3)
            %   constraints = cg.get_er_constr(ad_arr,sigma_i,Pi1,Pi2,Q,r,'Minimize Z2',Z1,Z2,Z3,mu2)
			%
			%Inputs:
			%	Q:	The "Q" variable in Q-Parameterization. Together with r, this defines a set of feedbacks
			%		(F,u_0) which can be determined by a nonlinear change of variables.
			%	r:	The "r" variable in Q-Parameterization. See note on Q.

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			cg = varargin{1};
			ad_arr = varargin{2};
			sigma_i = varargin{3};
			Pi1 = varargin{4};
			Pi2 = varargin{5};
			Q = varargin{6};
			r = varargin{7};
			in_str = varargin{8};

			switch in_str
			case 'Feasible Set'
				M1 = varargin{9};
				M2 = varargin{10};
				M3 = varargin{11};
			case 'Minimize M2'
				M1 = varargin{9};
				mu2 = varargin{10};
				M3 = varargin{11};
			case 'Minimize M3'
				M1 = varargin{9};
				M2 = varargin{10};
				mu3 = varargin{11};
			case 'Minimize Z2'
				P_M1 = varargin{9};
				Z2 = varargin{10};
				Z3 = varargin{11};
				mu2 = varargin{12};
			otherwise
				error('Unrecognized input string/type of equalized recovery problem.')
			end

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			n = size(ad_arr(1).A,1);
			vd = size(ad_arr(1).C_v,2);
			wd = size(ad_arr(1).B_w,2);

			unit_box = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));

			switch in_str
			case 'Feasible Set'
				if isa(M1,'Polyhedron') && isa(M2,'Polyhedron') && isa(M3,'Polyhedron')
					P_M1 = M1;
					P_M2 = M2;
					P_M3 = M3;
	            elseif isscalar(M1) && isscalar(M2) && isscalar(M3)
					P_M1 = M1 * unit_box;
					P_M2 = M2 * unit_box;
					P_M3 = M3 * unit_box;
				else
					error('Unrecognized input types for M1, M2, and M3.')
				end
			case 'Minimize M2'
				if isscalar(M1) && isscalar(M3)
					P_M1 = M1 * unit_box;
					P_M2 =  unit_box;
					P_M3 = M3 * unit_box;
				else
					error('Unrecognized input types for M1 and M3.')
				end
			case 'Minimize M3'
				if isscalar(M1) && isscalar(M2)
					P_M1 = M1 * unit_box;
					P_M2 = M2 * unit_box;
					P_M3 = unit_box;
				else
					error('Unrecognized input types for M1 and M2.')
				end
			otherwise
				error('Unexpected input string.')
			end
				
			T_i = length(sigma_i);

			%Select matrix
			select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

			sel_influenced_states = [];
			prod_M2 = 1;
			for i = 1 : T_i
				%Select all influenced states
				sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
				%Create product for M2
				prod_M2 = prod_M2*P_M2;
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Creating Constraints %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%

			constraints = [];
            
			%Get Special Matrices
			[H0,S0,Cm0,J0,f_bar,B_w_big,C_v_big] = get_mpc_matrices(ad_arr,'word',sigma_i);
			P_wT = 1; P_vT = 1;
			for t_idx = 1:T_i
				P_wT = P_wT*ad_arr(sigma_i(t_idx)).P_w;
				P_vT = P_vT*ad_arr(sigma_i(t_idx)).P_v;
			end

			G = [ 	(eye(n*(T_i+1))+S0*Q*Cm0)*H0*B_w_big ...
					S0*Q*C_v_big ...
					(eye(n*(T_i+1))+S0*Q*Cm0)*J0 ];

			bounded_disturb_matrix = [ 	P_wT.A zeros(size(P_wT.A,1),vd*T_i+n) ;
										zeros(size(P_vT.A,1),size(P_wT.A,2)) P_vT.A zeros(size(P_vT.A,1),n) ;
										zeros(size(P_M1.A,1),(vd+wd)*T_i) P_M1.A ];
			switch in_str
			case 'Feasible Set'
				constraints = constraints + [ Pi1 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= prod_M2.b - prod_M2.A*sel_influenced_states*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];
				constraints = constraints + [ Pi2 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= P_M3.b - P_M3.A*select_m(T_i,T_i)*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];

				constraints = constraints + [Pi1 * bounded_disturb_matrix == prod_M2.A*sel_influenced_states*G ];
				constraints = constraints + [Pi2 * bounded_disturb_matrix == P_M3.A*select_m(T_i,T_i)*G];
			case 'Minimize M2'
				constraints = constraints + [ Pi1 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= mu2*prod_M2.b - prod_M2.A*sel_influenced_states*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];
				constraints = constraints + [ Pi2 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= P_M3.b - P_M3.A*select_m(T_i,T_i)*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];

				constraints = constraints + [Pi1 * bounded_disturb_matrix == prod_M2.A*sel_influenced_states*G ];
				constraints = constraints + [Pi2 * bounded_disturb_matrix == P_M3.A*select_m(T_i,T_i)*G];
			case 'Minimize M3'
				constraints = constraints + [ Pi1 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= prod_M2.b - prod_M2.A*sel_influenced_states*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];
				constraints = constraints + [ Pi2 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= mu3*P_M3.b - P_M3.A*select_m(T_i,T_i)*(S0*r+(eye(n*(T_i+1))+S0*Q*Cm0)*H0*f_bar) ];

				constraints = constraints + [Pi1 * bounded_disturb_matrix == prod_M2.A*sel_influenced_states*G ];
				constraints = constraints + [Pi2 * bounded_disturb_matrix == P_M3.A*select_m(T_i,T_i)*G];
			case 'Minimize Z2'

				[~,temp_constrs] = cg.create_sadraddini_AH_inclusion_constr(	S0*r{pattern_ind}+(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar, ...
																				G{pattern_ind} , ...
																				bounded_disturb_matrix,[ P_wT.b ; P_vT.b ; P_M1.b ], ...
																				kron(ones(T_i+1,1),Z2.c),kron(eye(T_i+1),Z2.G), ...
																				[eye((T_i+1)*n);-eye((T_i+1)*n)],mu2*ones(2*(T_i+1)*n,1));

				constraints = constraints + temp_constrs;

				[~,temp_constrs] = cg.create_sadraddini_AH_inclusion_constr(	select_m(T_i,T_i)*(S0*r{pattern_ind}+(eye(n*(T_i+1))+S0*Q{pattern_ind}*Cm0)*H0*f_bar), ...
																				select_m(T_i,T_i)*G{pattern_ind} , ...
																				bounded_disturb_matrix,[ P_wT.b ; P_vT.b ; P_M1.b ], ...
																				Z3.c,Z3.G, ...
																				[eye(n);-eye(n)],ones(2*n,1) );
				constraints = constraints + temp_constrs;
			otherwise
				error('Unknown problem type.')
			end

		end

		function [constraints] = get_input_constr(varargin)
			%get_input_constr()
			%Description:
			%
			%Usage:
			%	[constraints] = cg.get_input_constr(ad_arr,sigma_i,Pi3,Q,r,M1,P_u)
			%	[constraints] = cg.get_input_constr(ad_arr,sigma_i,Pi3,Q,r,M1,P_u,'u_d',u_d)
			%
			%Assumptions:
			%	This function assumes that all dynamical systems in ad_arr have the same number of inputs.

			%%%%%%%%%%%%%%%%%%%%%%
			%% Input Processing %%
			%%%%%%%%%%%%%%%%%%%%%%

			cg = varargin{1};
			ad_arr = varargin{2};
			sigma_i = varargin{3};
			Pi3 = varargin{4};
			Q = varargin{5};
			r = varargin{6};
			M1 = varargin{7};
			P_u = varargin{8};

			if nargin > 8
				curr_idx = 9;
				while curr_idx <= nargin
					switch varargin{curr_idx}
						case 'u_d'
							u_d = varargin{curr_idx + 1};
							
							%Check dimension of u_d
							m = size(ad_arr(1).B,2);
							T = length(sigma_i);
							if size(u_d,1) ~= m*T
								error('u_d is not of the correct dimension.')
							end

							curr_idx = curr_idx+2;
						otherwise
							error('Unrecognized String.')
					end
				end
			else
				m = size(ad_arr(1).B,2);
				T = length(sigma_i);
				u_d = zeros(m*T,1);
			end

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			n = size(ad_arr(1).A,1);
			m = size(ad_arr(1).B,2);
			vd = size(ad_arr(1).C_v,2);
			wd = size(ad_arr(1).B_w,2);

			unit_box = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));

			if isa(M1,'Polyhedron')
				P_M1 = M1;
            elseif isscalar(M1)
				P_M1 = M1 * unit_box;
			else
				error('Unrecognized input types for M1.')
			end
				
			T_i = length(sigma_i);

			% %Select matrix
			% select_m = @(t,T_r) [zeros(n,t*n) eye(n) zeros(n,(T_r-t)*n) ];

			% sel_influenced_states = [];
			%prod_M2 = 1;
			P_uT = 1; P_wT = 1; P_vT = 1;
			for t_idx = 1 : T_i
				%Select all influenced states
				% sel_influenced_states = [ sel_influenced_states ; select_m(i,T_i) ];
				%Create product for M2
				% prod_M2 = prod_M2*P_M2;
				P_uT = P_uT*P_u;
				P_wT = P_wT*ad_arr(sigma_i(t_idx)).P_w;
				P_vT = P_vT*ad_arr(sigma_i(t_idx)).P_v;
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%%
			%% Creating Constraints %%
			%%%%%%%%%%%%%%%%%%%%%%%%%%

			constraints = [];

			%Get Special Matrices
			[H0,S0,Cm0,J0,f_bar,B_w_big,C_v_big] = get_mpc_matrices(ad_arr,'word',sigma_i);

			constraints = constraints + [ Pi3 * [ P_wT.b ; P_vT.b ; P_M1.b ] <= P_uT.b - P_uT.A*(r+Q*Cm0*H0*f_bar+u_d) ];

			G_tilde = [ Q*Cm0*H0*B_w_big, Q*C_v_big, Q*Cm0*J0 ];
			bounded_disturb_matrix = [ 	P_wT.A zeros(size(P_wT.A,1),vd*T_i+n) ;
										zeros(size(P_vT.A,1),size(P_wT.A,2)) P_vT.A zeros(size(P_vT.A,1),n) ;
										zeros(size(P_M1.A,1),(vd+wd)*T_i) P_M1.A ];

			constraints = constraints + [ Pi3 * bounded_disturb_matrix == P_uT.A*G_tilde];

		end

	end
end
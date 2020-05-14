function [results] = observer_comparison72( varargin )
	%observer_comparison71.m
	%Description:
	%	Testing various systems with the algorithms that have been developed so far.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	% if nargin >= 1
	% 	load_data_flag = varargin{1};
	% end

	% if nargin >= 3
	% 	c_sq.dim_x = varargin{2};
	% 	c_sq.dim_y = varargin{3};
	% end

	% if nargin >= 4
	% 	verbosity = varargin{4};
	% end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison72';
	save_file_name = [ 'results/' test_name '_results.mat'];

	dim = 3;

	PA = Polyhedron('lb',-ones(1,dim),'ub',ones(1,dim));

	A = PA.A(:,[1:2]);
	Ap = PA.A(:,3);
	a = -PA.b;

	for row_idx = [0:(2^dim)-1]
		
		%Create a binary vector to determine which value goes in each rows.
		temp_bin_vector = NaN(1,dim);
		for place_idx = 1:dim
			temp_bin_vector(dim+1-place_idx) = floor(row_idx/(2^place_idx));
		end

		B_tot(row_idx+1,:) = temp_bin_vector*1 + (~temp_bin_vector)*(-1);
	end

	B = B_tot(:,[1:2]);
	Bp = B_tot(:,3);

	b = -ones(size(B,1),1);

	results.A = A;
	results.Ap = Ap;
	results.a = a;

	results.B = B;
	results.Bp = Bp;
	results.b = b;

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	disp(['Beginning ' test_name '.' ])

	disp('1. Analyzing the kernel of B-prime.')
	kerB = null(Bp');
	kerB

	results.kerB = kerB;

	disp('  - Is there a nonempty intersection between this and R_+^l?')

	nonempty_intersection_flag = any( all( kerB' >= 0 ) );
	if nonempty_intersection_flag
		disp('    + Yes!')
	else
		disp('    + No!')
	end

	%% Create Optimization?
	disp('2. Attempting the full, bilinear optimization.')

	dim_struct.k = size(a,1);
	dim_struct.d = size(A,2);
	dim_struct.m = size(Ap,2);
	dim_struct.l = size(b,1);
	dim_struct.n = size(Bp,2);

	eps0 = 10^(-3);

	x = sdpvar(dim_struct.d,1);
	y = sdpvar(dim_struct.m,1);
	z = sdpvar(dim_struct.l,1);

	objective = [];

	projA_constr = [a + A*x + Ap*y >= 0];
	z_domain_constr = [Bp'*z == 0] + [ones(dim_struct.l,1)'*z == 1] + [z >= 0];
	containment_violation_constr = [z'*(b + B*x) <= -eps0];

	ops0 = sdpsettings('verbose',1);

	diagnostics = optimize(	projA_constr+z_domain_constr+containment_violation_constr, ...
							objective, ...
							ops0 );
	
	results.exp2.diagnostics = diagnostics;
	if diagnostics.problem ~= 0
		disp('  + The optimization was not feasible!')
	else
		disp('  + The optimization was feasible!')
	end

	%% Sufficient Optimization Condition %%

	clear objective containment_violation_constr diagnostics x y z

	disp('3. Attempting the full, bilinear optimization.')

	dim_struct.k = size(a,1);
	dim_struct.d = size(A,2);
	dim_struct.m = size(Ap,2);
	dim_struct.l = size(b,1);
	dim_struct.n = size(Bp,2);

	eps0 = 10^(-3);

	x = sdpvar(dim_struct.d,1);
	y = sdpvar(dim_struct.m,1);
	z = sdpvar(dim_struct.l,1);

	objective = [];

	projA_constr = [a + A*x + Ap*y >= 0];
	z_domain_constr = [Bp'*z == 0] + [ones(dim_struct.l,1)'*z == 1] + [z >= 0];
	containment_violation_constr = [b + B*x <= -eps0];

	ops0 = sdpsettings('verbose',1);

	diagnostics = optimize(	projA_constr+z_domain_constr+containment_violation_constr, ...
							objective, ...
							ops0 );
	
	results.exp3.diagnostics = diagnostics;
	if diagnostics.problem ~= 0
		disp('  + The optimization was not feasible!')
	else
		disp('  + The optimization was feasible!')
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	save( save_file_name ,'results')

end
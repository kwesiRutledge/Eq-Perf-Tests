function [results] = observer_comparison75( varargin )
	%observer_comparison75.m
	%Description:
	%	Comparing the method for detecting projection inclusion that I have against sadra's condition.
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

	test_name = 'observer_comparison75';
	save_file_name = [ 'results/' test_name '_results.mat'];

	dim = 2;
	projected_dims = 1;

	R_pd = zeros(length(projected_dims),dim);
	for dim_idx = [1:length(projected_dims)]
		R_pd(dim_idx,projected_dims(dim_idx)) = 1;
	end

	num_comparisons = 100;

	verbosity = 0;
	update_freq = 50;

	results.parameters.dim = dim;
	results.parameters.num_comparisons = num_comparisons;
	results.parameters.R_pd = R_pd;

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	disp(['Beginning ' test_name '.' ])
	disp(' ')

	%% Construct num_comparisons Polyhedra pairs %%

	disp(['Creating ' num2str(num_comparisons) ' Polytope pairs.' ])
	disp(' ')

	tic;
	PolytopePairs = {};
	while length(PolytopePairs) < num_comparisons

		%% Randomly Create the two Polyhedra.
		H_x = unifrnd(-1,1,2*dim,dim);
		h_x = ones(2*dim,1);

		H_y = unifrnd(-1,1,2*dim,dim);
		h_y = ones(2*dim,1);

		P_x = Polyhedron('A',H_x,'b',h_x);
		P_y = Polyhedron('A',H_y,'b',h_y);

		%% Check to see that P_x is full dimensional and nonempty.
		P_x_is_acceptable = P_x.isBounded() && (~P_x.isEmptySet());
		P_y_is_acceptable = P_y.isBounded() && (~P_y.isEmptySet());

		if P_x_is_acceptable && P_y_is_acceptable
			PolytopePairs{length(PolytopePairs)+1} = [ P_x , P_y ];
		else
			continue;
		end
	end
	results.timing.creating_pairs = toc;
	disp(['- It took ' num2str(results.timing.creating_pairs) ' seconds to create all Polytope Pairs.'])

	%% Use Sadraddini's Condition to Identify All Pairs that Contain One Another %%

	disp('Using Sadraddini''s Condition to verify if the projection polytope')
	disp('contains the projection of the other.')

	tic;

	cg = constr_gen(0);
	ops0 = sdpsettings('verbose',verbosity,'cachesolvers',1);

	feas_list = [];
	for pair_idx = [1:num_comparisons]
		P_x = PolytopePairs{pair_idx}(1);
		P_y = PolytopePairs{pair_idx}(2);

		[dual_variables,sadra_constrs] = cg.create_sadraddini_AH_inclusion_constr( ...
											zeros(length(projected_dims),1), R_pd , P_x.A, P_x.b, ...
											zeros(length(projected_dims),1), R_pd , P_y.A, P_y.b);

		%Optimize
		diagnostics = optimize(sadra_constrs,[],ops0);
		feas_list = [feas_list; (diagnostics.problem == 0) ];

		%Send Updates
		if mod(pair_idx,update_freq) == 0
			disp(['- Completed ' num2str(pair_idx) ' tests.'])
		end
	end
	results.timing.sadra_checks = toc;
	disp(['- It took ' num2str(results.timing.sadra_checks) ' seconds to check Sadra''s Containment Condition'])
	results.sadra_feas_list = feas_list;
	disp(' ')

	%% Use Our Condition to Identify All Pairs that Contain One Another %%

	disp('Using my condition to verify if the projected polytope contains the')
	disp('contains the projection of the other.')

	tic;

	feas_list = [];
	range_checklist = [];
	for pair_idx = [1:num_comparisons]
		% Collect Polytopes
		P_x = PolytopePairs{pair_idx}(1);
		P_y = PolytopePairs{pair_idx}(2);

		H_x = P_x.A;
		h_x = P_x.b;
		H_y = P_y.A;
		h_y = P_y.b;

		% Partition H_y
		H_y_prime = H_y(:,[1:projected_dims]);
		H_y_dblprime = H_y(:,[projected_dims+1:end]);

		% Verify the range of H_y
		H_yp_subset = H_y_prime';
		H_yp_rest = H_y_dblprime';
		col_count = 0;
		for col_idx = 1:size(H_yp_subset,2)
			col_out = H_yp_subset(:,col_idx);
			temp_H = H_yp_subset(:,[1:size(H_yp_subset,2)] ~= col_idx);	

			[xi,fvali,exitflagi,outputi] = constrained_nonnegative_ls( temp_H , -col_out , ...
																	[] , [] , ...
																	H_yp_rest , zeros(size(H_yp_rest,1),1) , ...
																	'verbosity', 0 );
			fval_list(col_idx) = fvali;
			if fvali < 1e-4
				col_count = col_count + 1;
			end
		end
		range_checklist = [ range_checklist ; (col_count == size(H_yp_subset,2) ) ];

		% Create Optimization Variables and Constraints
		q_x = size(H_x,1);

		[E,S,V] = svd(H_y_dblprime,0);
		P = E*((E'*E)^(-1))*E';

		Lambda = sdpvar( size(P,1) , q_x , 'full');

		nonneg_constr = [ Lambda >= 0 ];

		inclusion_constrs = [ Lambda*H_x == P*H_y_prime * R_pd ] + ...
							[ Lambda*h_x <= P*h_y ];

		% Optimize
		diagnostics = optimize(nonneg_constr+inclusion_constrs,[],ops0);
		feas_list = [feas_list; (diagnostics.problem == 0) ];

		%Send Updates
		if mod(pair_idx,update_freq) == 0
			disp(['- Completed ' num2str(pair_idx) ' tests.'])
		end

	end
	results.timing.my_checks = toc;
	disp(['- It took ' num2str(results.timing.my_checks) ' seconds to check My Containment Condition'])
	results.my_feas_list = feas_list;
	results.range_checklist = range_checklist;
	
	disp(' ')

	%% Use A Second Condition to Identify All Pairs that Contain One Another %%

	disp('Using an exhaustive condition to verify if the projected polytope contains the')
	disp('contains the projection of the other.')

	tic;

	feas_list = [];
	for pair_idx = [1:num_comparisons]
		% Collect Polytopes
		P_x = PolytopePairs{pair_idx}(1);
		P_y = PolytopePairs{pair_idx}(2);

		H_x = P_x.A;
		h_x = P_x.b;
		H_y = P_y.A;
		h_y = P_y.b;

		% Partition H_y
		H_y_prime = H_y(:,[1:projected_dims]);
		H_y_dblprime = H_y(:,[projected_dims+1:end]);

		% Create Optimization Variables and Constraints
		q_x = size(H_x,1);

		[q_y,n_y] = size(H_y);
		[q_x,n_x] = size(H_x);

		Gamma0 	= sdpvar(n_y,n_x,'full');
		beta0 	= sdpvar(n_y,1,'full');
		Lambda0 = sdpvar(q_y,q_x,'full');

		dual_vars = {Gamma0,beta0,Lambda0};

		% Create Constraint

		constrs = [];

		constrs = constrs + [Lambda0 >= 0];

		% constrs = constrs 	+ [ R_pd == R_pd*Gamma0 ] ...
		% 					+ [ Lambda0*H_x == (R_pd*H_y')'*R_pd*pinv(H_y) * H_y * Gamma0 ] ...
		% 					+ [ Lambda0*h_x <= (R_pd*H_y')'*R_pd*pinv(H_y) * h_y ];

		for hyperpl_idx = 1:size(H_y,1)
			c = R_pd*H_y(hyperpl_idx,:)';
			
			% disp(all(c'*R_pd*pinv(H_y) >= 0))

			% constrs = constrs + [ R_pd == R_pd*Gamma0 ] ...
			% 				+ [ Lambda0(hyperpl_idx,:)*H_x == c'*R_pd*pinv(H_y) * H_y * Gamma0 ] ...
			% 				+ [ Lambda0(hyperpl_idx,:)*h_x <= c'*R_pd*pinv(H_y) * h_y ];

			[u,~,exitflag,output] = nonnegative_ls(H_y',R_pd'*c );

			constrs = constrs + [ R_pd == R_pd*Gamma0 ] ...
							+ [ Lambda0(hyperpl_idx,:)*H_x == u'*H_y * Gamma0 ] ...
							+ [ Lambda0(hyperpl_idx,:)*h_x <= u'* h_y ];

		end

		% Optimize
		diagnostics = optimize(constrs,[],ops0);
		feas_list = [feas_list; (diagnostics.problem == 0) ];

		% Send Updates
		if mod(pair_idx,update_freq) == 0
			disp(['- Completed ' num2str(pair_idx) ' tests.'])
		end

	end
	results.timing.my_checks2 = toc;
	disp(['- It took ' num2str(results.timing.my_checks2) ' seconds to check My Second Containment Condition'])
	results.my_feas_list2 = feas_list;
	
	disp(' ')

	%%%%%%%%%%%%%%%%%%%%%
	%% Analyze Results %%
	%%%%%%%%%%%%%%%%%%%%%

	matching_data = (results.my_feas_list == results.sadra_feas_list);
	disp(['The checks matched ' num2str(sum(matching_data)) ' out of ' num2str(length(matching_data)) ' times.'])

	% Plot interesting mismatches
	interesting_idx = find(matching_data == 0,1);

	disp([' - For example consider case #' num2str(interesting_idx) ])
	disp([' - Result of Sadra''s check is ' num2str(results.sadra_feas_list(interesting_idx)) ])
	disp([' - Result of my check is ' num2str(results.my_feas_list(interesting_idx)) ])

	P_x = PolytopePairs{interesting_idx}(1);
	P_y = PolytopePairs{interesting_idx}(2);

	figure;
	hold on;
	plot(P_x,'Color','magenta')
	plot(P_y,'Color','cyan')


	%Checking how correct the algorithm was on polytopes with the desired
	range_satisfying_data = ( results.my_feas_list( range_checklist == 1 ) == results.sadra_feas_list( range_checklist == 1 ) );
	disp(['The checks matched (on all instances with proper range) ' ...
			num2str( sum(range_satisfying_data) ) ' out of ' num2str(length(range_satisfying_data)) ...
			' times.'  ])


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%save([ save_file_name '.mat'])

end
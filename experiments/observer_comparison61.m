function [results] = observer_comparison60( varargin )
	%observer_comparison60.m
	%Description:
	%	Testing various systems with the algorithms that have been developed so far.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin >= 1
		load_data_flag = varargin{1};
	end

	if nargin >= 3
		c_sq.dim_x = varargin{2};
		c_sq.dim_y = varargin{3};
	end


	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	if ~exist('c_sq')
		c_sq.dim_x = 2;
		c_sq.dim_y = 2;
	end
	
	dt = 0.1;

	in_sys = get_consensus_dyn(c_sq.dim_x,c_sq.dim_y,dt);

	n_x = size(in_sys.Dyn(1).A,1);
	n_u = size(in_sys.Dyn(1).B,2);

	eta_u = 0.5;
	Pu = Polyhedron('lb',-eta_u*ones(1,n_u),'ub',eta_u*ones(1,n_u));

	eta_x = 1.0;
	Px0 = Polyhedron('lb',-eta_x*ones(1,n_x),'ub',eta_x*ones(1,n_x));	

	results.params.sys = in_sys;
	results.params.Pu = Pu;
	results.params.Px0 = Px0;

	verbosity = 1; %Verbosity of Functions. Gives debugging info

	%Data file parameters.
	save_file_name = 'data/oc61_interm_results.mat';
	if ~exist('load_data_flag')
		load_data_flag = false;
	end
	run_opt_flag = false;

	if load_data_flag
		load(save_file_name)
		load_data_flag = true;
	end

	cg = constr_gen();

	ops = sdpsettings('verbose',verbosity);

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	disp('Test 1: Do a comparison between the time it takes to run the consistency set operation for two different words q1 and q2 and compare that with the language {q1,q2}.')
	disp(' ')
	q1 = in_sys.L.words{1};
	q2 = in_sys.L.words{2};
	q3 = in_sys.L.words{3};

	L2 = Language({q1,q2});
	disp('Created words and sublanguage.')

	t = 1;

	timing_at_t = [];
	for t = 1:5
		tic;
		[C_q1, ~] = in_sys.consistent_set(t,Language({q1}),Pu,Px0);
		time_for.C_q1 = toc;

		tic;
		[C_q2, ~] = in_sys.consistent_set(t,Language({q2}),Pu,Px0);
		time_for.C_q2 = toc;

		tic;
		[C_q3, ~] = in_sys.consistent_set(t,Language({q3}),Pu,Px0);
		time_for.C_q3 = toc;

		tic;
		[C_L2, ~] = in_sys.consistent_set(t,in_sys.L,Pu,Px0);
		time_for.C_L = toc;

		tic;
		C_L2_prime = C_q1.intersect(C_q2);
		C_L_prime = C_L2_prime.intersect(C_q3);
		time_for.C_L_prime = toc;

		%Text Update
		disp(['Completed timing at t = ' num2str(t) '.'])

		%Save timing data to an array of structs
		timing_at_t = [timing_at_t,time_for];
	end

	results.timing_at_t = timing_at_t;

	disp('Observe that the use of intersetions instead of the fancy call to consistent_set is faster.')
	disp('Probably because the large A matrix in the Ax <= b description of the Polyhedron is processed by MPT3.')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Simulate the System With this Form of Controller %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% figure;
	% plot3(x(1,:),x(2,:),x(3,:))

	% figure;
	% plot(x(:,1),x(:,2))

	% T = size(x,2);
	% figure;
	% hold on;
	% plot([1:T],x(1,:))
	% plot([1:T],x(2,:))
	% plot([1:T],x(3,:))

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
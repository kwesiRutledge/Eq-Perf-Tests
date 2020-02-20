function [results] = observer_comparison64( varargin )
	%observer_comparison64.m
	%Description:
	%	Testing the generation of a cover using the modified post algorithm.
	%	We will compute the belief graph 3 different times.
	%	- 

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin >= 4
		test_flags.simple = varargin{1};
		test_flags.simple_accelerated = varargin{2};
		test_flags.no_proj = varargin{3};
		test_flags.no_proj_accel = varargin{4};
	end

	if nargin >= 5
		test_flags.skip_unobs = varargin{5};
	end

	if ~exist('test_flags')
		test_flags.simple = false;
		test_flags.simple_accelerated = true;
		test_flags.no_proj = false;
		test_flags.no_proj_accel = false;
		test_flags.skip_unobs = true;
	end

	if ~isfield(test_flags,{'skip_unobs'})
		test_flags.skip_unobs = false;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;

	L1 = Language([1,1,2,3],[1,2,2,4],[1,2,4,3]);

	%Create a simple Language Constrainted Switching System
	A1 = eye(dim);
	B1 = eye(dim);
	C1 = eye(dim);
	f1 = [0;1];
	
	eta_v = 0.2; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));


	ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

	f2 = [1;0];
	f3 = -f1;
	f4 = -f2;

	lcss2 = [	ad1,...
				Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

	% Select Matrix
	RT = @(n,t,T) [zeros(n,n*t) eye(n) zeros(n,n*(T-t))];

	Hc = {}; Sc = {}; Jc = {}; fc = {};

	for word_idx = 1:length(L1.words)
		[Hc{word_idx},Sc{word_idx},Cc{word_idx},Jc{word_idx},fc{word_idx}] = get_mpc_matrices(lcss2,'word',L1.words{word_idx});

		%Create Cartesian Product of Polyhedra
		P_wT{word_idx} = 1; P_vT{word_idx} = 1; P_uT{word_idx} = 1;
		for sym_idx = 1:length(L1.words{word_idx})
			P_wT{word_idx} = P_wT{word_idx} * lcss2(L1.words{word_idx}(sym_idx)).P_w;
			P_vT{word_idx} = P_vT{word_idx} * lcss2(L1.words{word_idx}(sym_idx)).P_v;
			P_uT{word_idx} = P_uT{word_idx} * P_u;
		end
	end

	T = L1.find_longest_length();

	%%%%%%%%%%%%%%%%%%%%
	%% Plot Behaviors %%
	%%%%%%%%%%%%%%%%%%%%

	figure;
	for seq_ind = 1:L1.cardinality()
		
		X = (Hc{seq_ind}*(P_wT{seq_ind}+fc{seq_ind})) + Jc{seq_ind}*P_x0;
		X.computeVRep;
		Y = Cc{seq_ind}*X + P_vT{seq_ind};

		subplot(1,length(Hc),seq_ind)
		hold on;
		plot(RT(dim,0,T-1)*Y,'Color','white')
		plot(RT(dim,1,T-1)*Y,'Color','cyan')
		plot(RT(dim,2,T-1)*Y,'Color','magenta')
		plot(RT(dim,3,T-1)*Y,'Color','orange')
		%title(['Trajectory of Word ' num2str(seq_ind)])
		xlabel('$x_1(t)$','Interpreter','latex')
		ylabel('$x_2(t)$','Interpreter','latex')
		axis([-2 3 -2 3])
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph Using the Simple Algorithm %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if test_flags.simple
		disp('- BeliefGraph Approach #1')
		disp('  + No flags used.')

		timer_start = tic;
		bg = BeliefGraph(lcss2,L1, P_u, P_x0,'fb_method','output','verbosity',1);
		timing_info.create_bg_simple = toc(timer_start);

		results.bg_simple = bg;

		figure;
		bg.plot();

		disp(['  + Completed in ' num2str(timing_info.create_bg_simple) ' seconds.' ])
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph Using the Accelerated Algorithm %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if test_flags.simple_accelerated

		disp('- BeliefGraph Approach #2')
		disp('  + Accelerate Flag is True.')

		timer_start = tic;
		bg = BeliefGraph(lcss2,L1, P_u, P_x0,'fb_method','output','accel_flag',true,'verbosity',2);
		timing_info.create_bg_accelerated = toc(timer_start);

		results.bg_accelerated = bg;

		figure;
		bg.plot();

		disp(['  + Completed in ' num2str(timing_info.create_bg_accelerated) ' seconds.' ])
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph Using the No-Projection Algorithm %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if test_flags.no_proj
		disp('- BeliefGraph Approach #3')
		disp('  + No projections used.')
		disp('  + No acceleration used.')

		timer_start = tic;
		bg = BeliefGraph(lcss2,L1, P_u, P_x0,'fb_method','output','accel_flag',false,'use_proj_flag',false,'verbosity',0);
		timing_info.create_bg_no_proj = toc(timer_start);

		results.bg_no_proj = bg;

		figure;
		bg.plot();

		disp(['  + Completed in ' num2str(timing_info.create_bg_no_proj) ' seconds.' ])
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph Using the No-Projection Algorithm with Acceleration %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if test_flags.no_proj_accel
		disp('- BeliefGraph Approach #4')
		disp('  + No projections used.')
		disp('  + Acceleration is used.')

		timer_start = tic;
		bg = BeliefGraph(lcss2,L1, P_u, P_x0,'fb_method','output','accel_flag',true,'use_proj_flag',false,'verbosity',0);
		timing_info.create_bg_no_proj_accel = toc(timer_start);

		results.bg_no_proj_accel = bg;

		figure;
		bg.plot();

		disp(['  + Completed in ' num2str(timing_info.create_bg_no_proj_accel) ' seconds.' ])
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create Belief Graph Using the Standard Algorithm without Unobservability Checks %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if test_flags.skip_unobs
		disp('- BeliefGraph Approach #5')
		disp('  + Projection is used.')
		disp('  + No unobservability checks are used.')

		timer_start = tic;
		bg = BeliefGraph(lcss2,L1, P_u, P_x0,'fb_method','output','verbosity',2,'use_unobservability_checks',false);
		timing_info.create_bg_no_proj_accel = toc(timer_start);

		results.bg_skip_unobs_check = bg;

		figure;
		bg.plot();

		disp(['  + Completed in ' num2str(timing_info.create_bg_no_proj_accel) ' seconds.' ])
	end

	%%%%%%%%%%%%%
	%% Results %%
	%%%%%%%%%%%%%
	
	results.System = lcss2;
	results.TimingInfo = timing_info;
end
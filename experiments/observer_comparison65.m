function [results] = observer_comparison65( varargin )
	%observer_comparison65.m
	%Description:
	%	Testing the generation of a cover using the modified post algorithm.

	disp('observer_comparison65() started.')

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if ~exist('verbosity')
		verbosity = 1;
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;
	ss_side_length = 10;
	state_space = Polyhedron('lb',-(ss_side_length/2)*ones(1,dim),'ub',(ss_side_length/2)*ones(1,dim));

	cutout = Polyhedron('lb',-ones(1,dim),'ub',ones(1,dim));

	state_space_partition = state_space \ cutout;
	p1 = Polyhedron('A',[1,0],'b',-1);
	p2 = Polyhedron('A',[-1,0],'b',-1);
	p3 = Polyhedron('A',[0,1],'b',-1);
	p4 = Polyhedron('A',[0,-1],'b',-1);

	state_space_overlaps = PolyUnion([	p1.intersect(state_space), ...
										p2.intersect(state_space), ...
										p3.intersect(state_space), ...
										p4.intersect(state_space)])

	target_set = cutout + [-1;-3];

	disp('- Constants Loaded.')

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	for poly_idx = 1:length(state_space_overlaps.Set)
		q_y(poly_idx) = size(state_space_overlaps.Set(poly_idx).A,1);
	end

	q_x = size(target_set.A,1);
	
	%% Start Optimization

	select_bin = binvar(length(state_space_overlaps.Set),1);
	for lambda_num = 1:length(q_y)
		Lambda{lambda_num} = sdpvar(q_y(lambda_num),q_x,'full');
	end

	%Create Constraints
	constraints = [];
	for constraint_idx = 1:length(q_y)
		P_Y = state_space_overlaps.Set(poly_idx);
		P_X = target_set;

		constraints = constraints + [ select_bin(constraint_idx)* Lambda{constraint_idx}*P_X.A == select_bin(constraint_idx)* P_Y.A] + ...
									[ select_bin(constraint_idx)* Lambda{constraint_idx}*P_X.b <= select_bin(constraint_idx)* P_Y.b];
	end

	sum_constraint = [sum(select_bin) == 1];

	constraints = constraints + sum_constraint;

	%Optimization
	ops = sdpsettings('verbose',verbosity);
	optim0 = optimize(constraints,[],ops);

	%%%%%%%%%%%%%
	%% Results %%
	%%%%%%%%%%%%%
	
	results.StateSpace = state_space;
	results.Partition = state_space_partition;
	results.OverlappingCover = state_space_overlaps;
	results.Cutout = cutout;
	results.TargetSet = target_set;

	results.Test1.OptimizationData = optim0;
	results.Test1.BinaryVariables = value(select_bin);
end
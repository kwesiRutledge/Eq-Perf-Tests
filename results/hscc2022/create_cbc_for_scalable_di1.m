% create_cbc_for_scalable_di1.m
% Description:
%	Attempts to find a controller for the scalable double integrator system from Liren's paper.

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 5;

us_dim = 0;

target_lb = [ -62 , -10 , -65 , -10 , -65 , -10 , -ones(1,us_dim)*10^(4) ];
target_ub = [ -50 , 10 , -50 , 10 , -45 , 10 , ones(1,us_dim)*10^(4) ];
new_target = Polyhedron('lb',target_lb,'ub',target_ub);

eta_u = 10;
[ lcsas0 , P_u , Pw1 , x0 , P_target ] = get_hidden_uncontrollable_subspace_lcsas1(	'UncontrollableSubspaceDimension',us_dim, ...
																					'TimeHorizon',TimeHorizon,...
																					'eta_u',eta_u, ...
																					'P_target',new_target, ...
																					'scl',0.05);

%% Synthesis %%

[ scalable_controller , info ] = lcsas0.FindAdaptiveControllerWithMStar( P_target , ...
																	'SearchStrategy' , 'DescendingCardinality', ...
																	'DoOptimizationPruningWhere' , 'DuringSearch' , ...
																	'RemoveBilinearityInInputConstraints', true , ...
																	'RemoveBilinearityInReachabilityConstraints', true , ...
																	'LinearizeBilinearContainment', true )

%% Visualizing %%

successfullySolved = (info.About(end).problem == 0);
if successfullySolved

	results.SimulationData = [];

	axis_limits = [inf,-inf,inf,-inf,inf,-inf];
	eps1 = 5;

	figure;

	hold on;
	plot3(lcsas0.X0.V(1),lcsas0.X0.V(3),lcsas0.X0.V(5))
	plot(P_target.projection([1,3,5]),'alpha',0.2)

	for simulation_index = 1:10
		[ x_0_t, u_0_tm1 , y_0_t , sig ] = scalable_controller.simulate_1run();

		% for t = 0:TimeHorizon
		% 	x_t = x_0_t(:,t+1);
		% 	scatter(x_t(1),x_t(2))
		% end
		plot3(x_0_t(1,:),x_0_t(3,:),x_0_t(5,:))

		axis_limits(1) = min([x_0_t(1,:),axis_limits(1)]);
		axis_limits(2) = max([x_0_t(1,:),axis_limits(2)]);
		axis_limits(3) = min([x_0_t(3,:),axis_limits(3)]);
		axis_limits(4) = max([x_0_t(3,:),axis_limits(4)]);
		axis_limits(5) = min([x_0_t(5,:),axis_limits(5)]);
		axis_limits(6) = max([x_0_t(5,:),axis_limits(6)]);
		
		%Save data
		results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

	end
	
	axis(axis_limits + eps1*[-1,1,-1,1,-1,1])

	saveas(gcf,'images/scalable_example1_runs','epsc')
	saveas(gcf,'images/scalable_example1_runs','png')

end

%% Save Data

save(['data/scalable_example1_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
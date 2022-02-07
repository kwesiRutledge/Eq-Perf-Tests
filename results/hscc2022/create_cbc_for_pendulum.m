% create_cbc_for_pendulum.m
% Description:
%	

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 5;
new_target = Polyhedron('lb',[-0.2,-2],'ub',[0.2,2]);
eta_u = 10;
[ lcsas0 , ~ , ~ , eta_v , eta_x0 , P_target ] = get_uncertain_thick_pendulum_lcsas('TimeHorizon',TimeHorizon,'P_target',new_target,'eta_u',eta_u);

%% Synthesis %%

[ pendulum_controller , info ] = lcsas0.FindConsistentBeliefController( P_target , ...
																	'SearchStrategy' , 'DescendingCardinality', ...
																	'DoOptimizationPruningWhere' , 'DuringSearch' , ...
																	'RemoveBilinearityInInputConstraints', true , ...
																	'RemoveBilinearityInReachabilityConstraints', true , ...
																	'LinearizeBilinearContainment', true )

%% Visualizing %%

successfullySolved = (info.About(end).problem == 0);
if successfullySolved

	results.SimulationData = [];

	axis_limits = zeros(1,4);
	eps1 = 0.2;

	figure;

	hold on;
	plot(lcsas0.X0)
	plot(P_target)

	for simulation_index = 1:10
		[ x_0_t, u_0_tm1 , y_0_t , sig ] = pendulum_controller.simulate_1run();

		% for t = 0:TimeHorizon
		% 	x_t = x_0_t(:,t+1);
		% 	scatter(x_t(1),x_t(2))
		% end
		plot(x_0_t(1,:),x_0_t(2,:))

		axis_limits(1) = min([x_0_t(1,:),axis_limits(1)]);
		axis_limits(2) = max([x_0_t(1,:),axis_limits(2)]);
		axis_limits(3) = min([x_0_t(2,:),axis_limits(3)]);
		axis_limits(4) = max([x_0_t(2,:),axis_limits(4)]);
		

		%Save data
		results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

	end
	
	axis(axis_limits + eps1*[-1,1,-1,1])

	saveas(gcf,'images/pendulum_runs','epsc')
	saveas(gcf,'images/pendulum_runs','png')

end

%% Save Data

save(['data/pendulum_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
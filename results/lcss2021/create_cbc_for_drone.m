% create_cbc_for_drone.m
% Description:
%	

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 8;
[ lcsas0 , x0 , TimeHorizon , P_target ] = get_differently_loaded_drone_lcsas('TimeHorizon',TimeHorizon);

%% Synthesis %%

[ drone_controller , info ] = lcsas0.FindConsistentBeliefController( P_target , ...
																	'SearchStrategy' , 'AscendingCardinality', ...
                                                                    'GurobiNodeLimit' , 10^6, ...
																	'DoOptimizationPruningWhere' , 'DuringSearch' , ...
																	'RemoveBilinearityInInputConstraints', true , ...
																	'RemoveBilinearityInReachabilityConstraints', true , ...
																	'LinearizeBilinearContainment', false )

%% Visualizing %%

results.SimulationData = [];

figure;

hold on;
plot(lcsas0.X0)
plot(P_target)

for simulation_index = 1:10
	[ x_0_t, u_0_tm1 , y_0_t , sig ] = drone_controller.simulate_1run();

	% for t = 0:TimeHorizon
	% 	x_t = x_0_t(:,t+1);
	% 	scatter(x_t(1),x_t(2))
	% end
	plot(x_0_t(1,:),x_0_t(2,:))
	

	%Save data
	results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

end
axis([-0.5,12.5,-6,6.5])

saveas(gcf,'images/drone_runs','epsc')
saveas(gcf,'images/drone_runs','png')

save(['data/drone_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
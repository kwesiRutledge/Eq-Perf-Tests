%create_cbc_for_toy2.m
%Description:
%	Synthesizes a consistent belief controller for toy problem 2. Similar rotation

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 4;
[ lcsas0 , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_similar_rotation_lcsas('TimeHorizon',TimeHorizon,'eta_u',10);

%% Synthesis %%

[ toy2_controller , info ] = lcsas0.FindConsistentBeliefController( P_target , 'SearchStrategy' , 'AscendingCardinality' , ...
																				'DoOptimizationPruningWhere' , 'DuringSearch', ...
																				'GurobiNodeLimit' , 3*10^5 , ...
																				'RemoveBilinearityInInputConstraints', true , ...
																				'RemoveBilinearityInReachabilityConstraints', true );

%% Visualizing %%

results.SimulationData = [];

figure;

hold on;
plot(lcsas0.X0)
plot(P_target)

for simulation_index = 1:15
	[ x_0_t, u_0_tm1 , y_0_t , sig ] = toy2_controller.simulate_1run();

	% for t = 0:TimeHorizon
	% 	x_t = x_0_t(:,t+1);
	% 	scatter(x_t(1),x_t(2))
	% end
	plot(x_0_t(1,:),x_0_t(2,:))
	

	%Save data
	results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

end
axis([-0.5,12.5,-0.5,8.5])

saveas(gcf,'images/toy2_runs','epsc')
saveas(gcf,'images/toy2_runs','png')

save(['data/toy2_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
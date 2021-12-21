%create_cbc_for_toy2.m
%Description:
%	Synthesizes a consistent belief controller for toy problem 2. Similar rotation

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

% TimeHorizon = 4;
% eta_u = 4;
TimeHorizon = 6;
eta_u = 2;

[ lcsas0 , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_similar_rotation_lcsas('TimeHorizon',TimeHorizon,'eta_u',eta_u);

%% Synthesis %%

[ toy2_controller , info ] = lcsas0.FindConsistentBeliefController( P_target , 'SearchStrategy' , 'DescendingCardinality' , ...
																				'DoOptimizationPruningWhere' , 'DuringSearch', ...
																				'GurobiNodeLimit' , 2*10^5 , ...
																				'RemoveBilinearityInInputConstraints', true , ...
																				'RemoveBilinearityInReachabilityConstraints', true , ...
																				'LinearizeBilinearContainment', true );

%% Visualizing %%

results.SimulationData = [];
extreme_vector = zeros(1,4);

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
	
	%Update extreme_vector
	extreme_vector(1) = min([x_0_t(1,:),extreme_vector(1)]); % minimum of all first state
	extreme_vector(2) = max([x_0_t(1,:),extreme_vector(2)]); % maximum of all first state
	extreme_vector(3) = min([x_0_t(2,:),extreme_vector(3)]); % minimum of all second state
	extreme_vector(4) = max([x_0_t(2,:),extreme_vector(4)]); % maximum of all second state


	%Save data
	results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

end
axis( extreme_vector + [-0.5,0.5,-0.5,0.5])

saveas(gcf,'images/toy2_runs','epsc')
saveas(gcf,'images/toy2_runs','png')

save(['data/toy2_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
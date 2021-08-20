%create_cbc_for_toy1.m
%Description:
%	

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 4;
[ lcsas0 , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_opposing_rotation_lcsas('TimeHorizon',TimeHorizon);

%% Synthesis %%

[ toy1_controller , info ] = lcsas0.FindConsistentBeliefController( P_target );

%% Visualizing %%

results.SimulationData = [];

figure;

hold on;
plot(lcsas0.X0)
plot(P_target)

for simulation_index = 1:10
	[ x_0_t, u_0_tm1 , y_0_t , sig ] = toy1_controller.simulate_1run();
	% for t = 0:TimeHorizon
	% 	x_t = x_0_t(:,t+1);
	% 	scatter(x_t(1),x_t(2))
	% end
	plot(x_0_t(1,:),x_0_t(2,:))
	toy1_controller.clear_histories()

	%Save data
	results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

end
axis([-3,4.5,-2.5,2.5])

saveas(gcf,'images/toy1_runs','epsc')
saveas(gcf,'images/toy1_runs','png')

save(['data/toy1_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
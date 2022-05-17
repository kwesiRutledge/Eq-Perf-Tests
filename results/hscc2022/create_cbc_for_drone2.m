% create_cbc_for_drone2.m
% Description:
%	Attempts to find a controller for the simplified drone system example.

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 10; x0 = [0;0];
X_Target = Polyhedron('lb',[2,2],'ub',[3,3])
[ lcsas0 , x0 , TimeHorizon , P_target ] = get_tweaked_turn_drone_lcsas('TimeHorizon',TimeHorizon, 'x0',x0);

%% Synthesis %%

[ drone_controller , info ] = lcsas0.FindAdaptiveControllerWithMStar( P_target );

if strcmp(info.Message,'Solved!')

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
	axis([-0.5,5.5,-0.5,5.5])
	xlabel('$z$','Interpreter','latex')
	ylabel('$\dot{z}$','Interpreter','latex')

	saveas(gcf,'images/drone_turning_runs','epsc')
	saveas(gcf,'images/drone_turning_runs','png')
else
	disp('The drone controller was not successfully created.')
end

save(['data/drone_turning_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
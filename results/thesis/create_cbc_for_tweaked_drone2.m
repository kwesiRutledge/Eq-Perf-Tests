% create_cbc_for_tweaked_drone2.m
% Description:
%	Attempts to find a controller for the simplified drone system example
%	where there is uncertainty in whether or not the velocity vector input is rotated.
%	This example finds a controller that satisfies the \bigcirc^T X_T task.
% 	Successfully used on hardware.

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 8;
x0 = [0;0];
dt = 0.5;
X_Target = Polyhedron('lb',[1.2,1.0],'ub',[1.8,1.4])
eta_u = 0.5;

[ lcsas0 , x0 , TimeHorizon , P_target ] = get_tweaked_turn_drone_lcsas('TimeHorizon',TimeHorizon, 'x0',x0, ...
																		'eta_u', eta_u, 'X_target', X_Target, ...
																		'dt' , dt , 'eta_w' , 0.1 , ...
																		'theta2' , pi/16 );

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
	axis([-0.5,2.0,-0.5,2.0])
	xlabel('$z$','Interpreter','latex')
	ylabel('$\dot{z}$','Interpreter','latex')

	saveas(gcf,'images/drone_turning_runs','epsc')
	saveas(gcf,'images/drone_turning_runs','png')

	%Save controller data somewhere.
	drone_controller.deconstruct_and_save_to(['data/turning_controller2_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat'])

else
	disp('The drone controller was not successfully created.')
end

save(['data/drone_turning_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
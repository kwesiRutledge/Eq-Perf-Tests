% hw_kltl_create_cbc_for_tweaked_drone.m
% Description:
%	Attempts to find a controller for the simplified drone system example
%	where there is uncertainty in whether or not the velocity vector input is rotated.
%	This example finds a controller that satisfies the task:
%
%		/^\_{i = 1,2} <> K { \mode_i \} \implies <> X_i
%
% 	Targeted to run on hardware.
%
%Usage:
%
%
%Assumptions:
%

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 4;
x0 = [0;0];
dt = 1.0;
% X_Target = Polyhedron('lb',[1.2,1.0],'ub',[1.8,1.4])
X_Target = Polyhedron('lb',[1.2,0.0],'ub',[1.8,0.4])
eta_u = 0.5;

[ lcsas0 , x0 , TimeHorizon , P_target ] = get_tweaked_turn_drone_lcsas('TimeHorizon',TimeHorizon, 'x0',x0, ...
																		'eta_u', eta_u, 'X_target', X_Target, ...
																		'dt' , dt , 'eta_w' , 0.25 , ...
																		'theta2' , pi/16 );

% P_Targets = [ P_target , P_target - [0.7;0] ];

%% Synthesis %%

[ drone_controller , info ] = lcsas0.FindAdaptiveController_AlwaysHidenInfo( P_target );

date_string = datestr(now,'ddmmmyyyy-HHMM');
save(['data/kltl-drone_turning_data-hidden_info-' date_string '.mat' ])

if strcmp(info.Message,'Solved!')

	%% Visualizing %%

	results.SimulationData = [];

	figure;

	hold on;
	plot(P_target,'color','red','alpha',0.5)
% 	plot(P_targets(2),'color','magenta','alpha',0.5)

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

	axis([-0.5,2.5,-0.5,2.5])
	xlabel('$z$','Interpreter','latex')
	ylabel('$\dot{z}$','Interpreter','latex')

	saveas(gcf,'images/kltl-drone_turning_runs-hidden_info','epsc')
	saveas(gcf,'images/kltl-drone_turning_runs-hidden_info','png')

	%Save controller data somewhere.
	drone_controller.deconstruct_and_save_to(['data/kltl_turning_controller_data_' date+stromg '.mat'])

else
	disp('The drone controller was not successfully created.')
end

save(['data/kltl-drone_turning_data-hidden_info-' date_string '.mat' ])
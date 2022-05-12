% create_cbc_for_drone.m
% Description:
%	

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

%TimeHorizon = 8;
TimeHorizon = 8; x0 = [1;0]; X_target = Polyhedron('lb',[2],'ub',[3]) * Polyhedron('lb',-6,'ub',6);
%TimeHorizon = 16; x0 = [1;2]; X_target = Polyhedron('lb',[3],'ub',[9]) * Polyhedron('lb',-10,'ub',10);
% [ lcsas0 , x0 , TimeHorizon , P_target ] = get_differently_loaded_drone_lcsas('TimeHorizon',TimeHorizon,'m1',0.030,'m2',0.040);
[ lcsas0 , x0 , TimeHorizon , P_target ] = get_differently_loaded_drone_lcsas(	'TimeHorizon',TimeHorizon, ...
																				'm1',1.0,'m2',1.5, ...
																				'x0',x0, ...
																				'X_target' , X_target );

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
	axis([-0.5,12.5,-5.5,5.5])
	xlabel('$z$','Interpreter','latex')
	ylabel('$\dot{z}$','Interpreter','latex')

	saveas(gcf,'images/drone_runs','epsc')
	saveas(gcf,'images/drone_runs','png')
else
	disp('The drone controller was not successfully created.')
end

save(['data/drone_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
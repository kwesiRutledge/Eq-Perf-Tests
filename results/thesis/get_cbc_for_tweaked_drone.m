%get_cbc_for_tweaked_drone.m
%Description:
%	Defines a 

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

TimeHorizon = 10; 
x0 = [0;0];
eta_u = 3;
X_Target = Polyhedron('lb',[2.5,2.5],'ub',[4,4]);

[ lcsas0 , x0 , TimeHorizon , P_target ] = get_tweaked_turn_drone_lcsas('TimeHorizon',TimeHorizon, 'x0',x0, ...
																		 'eta_u', eta_u, 'X_target', X_Target );
P_Targets = [ P_target , P_target - [2.5;0] ];

%% Synthesis %%

[ drone_controller , info ] = lcsas0.FindAdaptiveControllerWithMStar_2Goals( P_Targets );

if strcmp(info.Message,'Solved!')

	%% Visualizing %%

	results.SimulationData = [];

	figure;

	hold on;
	plot(lcsas0.X0)
	plot(P_Targets(1),'color','red','alpha',0.5)
	plot(P_Targets(2),'color','magenta','alpha',0.5)

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
	xlabel('$p_x$','Interpreter','latex','FontSize',20)
	ylabel('$p_y$','Interpreter','latex','FontSize',20)

	legend('','$X_{T}^{(1)}$','$X_{T}^{(2)}$','Interpreter','latex','FontSize',20)

	saveas(gcf,'images/drone_turning_runs','epsc')
	saveas(gcf,'images/drone_turning_runs','png')
else
	disp('The drone controller was not successfully created.')
end

save(['data/drone_turning_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
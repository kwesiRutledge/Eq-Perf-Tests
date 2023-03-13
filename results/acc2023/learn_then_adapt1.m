%get_cbc_for_tweaked_drone.m
%Description:
%	Defines a 

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

rng("default")

TimeHorizon = 9; 
x0 = [0;0];
eta_u = 0.5;
eta_w = 0.1;
X_Target = Polyhedron('lb',[1.2,1.0],'ub',[1.8,1.4]);

[ lcsas0 , x0 , TimeHorizon , P_target ] = get_tweaked_turn_drone_lcsas('TimeHorizon',TimeHorizon, 'x0',x0, ...
																		 'eta_u', eta_u, 'X_target', X_Target, ...
                                                                         'eta_w', eta_w);
P_Targets = [ P_target , P_target - [0.7;0] ];

%% Synthesis %%

[ drone_controller , info ] = lcsas0.FindAdaptiveControllerWithMStar_2Goals( P_Targets );

if strcmp(info.Message,'Solved!')

	%% Visualizing %%

    trajectory_lw = 1.5; % Linewidth of trajectories
    axis_fs = 24; %default 11

	results.SimulationData = [];

	figure;

	hold on;
	plot(lcsas0.X0)
	plot(P_Targets(1),'color','red','alpha',0.5)
	plot(P_Targets(2),'color','magenta','alpha',0.5)

    mode1_count = 0;
    mode2_count = 0;
	for simulation_index = 1:10
		[ x_0_t, u_0_tm1 , y_0_t , sig ] = drone_controller.simulate_1run();

		if sig(1) == 1
            %trajectory_color = 'yellow';
            trajectory_color = "#6200EE";
 
            mode1_count = mode1_count + 1;
        elseif sig(1) == 2
            %trajectory_color = 'blue';
            trajectory_color = "#03DAC6";

            mode2_count = mode2_count + 1;
        else
            error([ 'Unexpected sig value! ' num2str(sig) ])
        end

		plot(x_0_t(1,:),x_0_t(2,:), ...
            'LineWidth',trajectory_lw, ...
            'Color', trajectory_color)
		

		%Save data
		results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

	end
	axis([-0.5,1.8,-0.5,1.4])
	xlabel('$p_x$','Interpreter','latex','FontSize',axis_fs)
	ylabel('$p_y$','Interpreter','latex','FontSize',axis_fs)

	legend('','$X_{T}^{(1)}$','$X_{T}^{(2)}$','Interpreter','latex','FontSize',20)

	saveas(gcf,'images/drone_learn_then_adapt_runs','epsc')
	saveas(gcf,'images/drone_learn_then_adapt_runs','png')
else
	disp('The drone controller was not successfully created.')
end

save(['data/drone_turning_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
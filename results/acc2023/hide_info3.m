%hide_info3.m
%Description:
%	This test will attempt to synthesize a controller that satisfies the hidden mode specification.
%	This example uses a simple 2d drone model with a wacky B matrix.
%   While one input needs to remain "0" to keep th emode of the system
%   "hidden" from outside observers, the other input needs to aggressively
%   act to satisfy a task.
%
%	The task:
%
%		/^\_{i = 1,2} <> K { \mode_i \} \implies <> X_i
%

close all; clear all; clc; 

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

rng("default")

[ lcsas0 , x0 , TimeHorizon , P_target ] = get_2d_lcsas_for_hiding( ...
    'eta_u',2.0, ...
    'W',Polyhedron('lb',[0.8,-0.1],'ub',[1.2,0.1]));

% P_Targets = [ P_target , P_target - [0.7;0] ];

%% Synthesis %%

[ scalar_system_controller , info ] = lcsas0.FindAdaptiveController_AlwaysHidenInfo( P_target , 'RemoveBilinearityInInputConstraints' , true );

date_string = datestr(now,'ddmmmyyyy-HHMM');
save(['data/scalar_sys_data-hidden_info-' date_string '.mat' ])

if strcmp(info.Message,'Solved!')

	%% Visualizing %%

    trajectory_lw = 1.5; % Linewidth of trajectories
    x0_markersize = 108; %72;
    fs = 21;
    legend_fs = 24;
    axis_fs = 24; %default 11

	results.SimulationData = [];

	figure;

	hold on;
	target_plot = plot(P_target,'color','red','alpha',0.5);
% 	plot(P_targets(2),'color','magenta','alpha',0.5)
    ic_plot = scatter(lcsas0.X0.V(1),lcsas0.X0.V(2),x0_markersize,'black','filled');

    mode1_count = 0;
    mode2_count = 0;
	for simulation_index = 1:10
		[ x_0_t, u_0_tm1 , y_0_t , sig ] = scalar_system_controller.simulate_1run();

		if sig(1) == 1
            %trajectory_color = 'yellow';
            % trajectory_color = "#6200EE";
            trajectory_color = [98,0,238]/255;
            trajectory_color(3) = trajectory_color(3) - (40/255) * mode1_count;
            mode1_count = mode1_count + 1;
        elseif sig(1) == 2
            %trajectory_color = 'blue';
            %trajectory_color = "#03DAC6";
            trajectory_color = [3, 218, 198]/255;
            trajectory_color(2) = trajectory_color(2) - (30/255)*mode2_count;

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

	axis([-0.5,3.0,-1.0,3.0]) 
    xticks([0 1 2])

    xlabel('$p^{(x)}$','Interpreter','latex','FontSize',axis_fs)
	ylabel('$p^{(y)}$','Interpreter','latex','FontSize',axis_fs)

    legend( ...
        [target_plot,ic_plot],'$$X_{T}$$','$$x_0$$', ...
        'Interpreter','latex', ...
        'FontSize',legend_fs, ...
        'Location','southeast' ...
        )

	saveas(gcf,'images/kltl-scalar_runs-hidden_info2','epsc')
	saveas(gcf,'images/kltl-scalar_runs-hidden_info2','png')
    

	%Save controller data somewhere.
	scalar_system_controller.deconstruct_and_save_to(['data/scalar_controller2_data_' date_string '.mat'])

else
	disp('The drone controller was not successfully created.')
end

save(['data/scalar_sys_data-hidden_info-' date_string '.mat' ])	
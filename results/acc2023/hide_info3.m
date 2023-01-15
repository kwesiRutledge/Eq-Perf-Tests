%hide_info2.m
%Description:
%	This test will attempt to synthesize a controller that satisfies the hidden mode specification.
%	It doesn't use the drone model, but instead a simple scalar system.
%	The task:
%
%		/^\_{i = 1,2} <> K { \mode_i \} \implies <> X_i
%

close all; clear all; clc; 

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

[ lcsas0 , x0 , TimeHorizon , P_target ] = get_scalar_lcsas_for_hiding2d('eta_u',0.1);

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

	for simulation_index = 1:10
		[ x_0_t, u_0_tm1 , y_0_t , sig ] = scalar_system_controller.simulate_1run();

		% for t = 0:TimeHorizon
		% 	x_t = x_0_t(:,t+1);
		% 	scatter(x_t(1),x_t(2))
		% end
		plot(x_0_t(1,:),x_0_t(2,:),'LineWidth',trajectory_lw)
		

		%Save data
		results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

	end

	axis([-0.5,2.5,-0.5,3.0])
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
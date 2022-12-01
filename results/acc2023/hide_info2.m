%hide_info2.m
%Description:
%	This test will attempt to synthesize a controller that satisfies the hidden mode specification.
%	It doesn't use the drone model, but instead a simple scalar system.
%	The task:
%
%		/^\_{i = 1,2} <> K { \mode_i \} \implies <> X_i
%

%% Add Libraries %%

addpath(genpath('../../functions/'))

%% Constants %%

[ lcsas0 , x0 , TimeHorizon , P_target ] = get_scalar_lcsas_for_hiding();

% P_Targets = [ P_target , P_target - [0.7;0] ];

%% Synthesis %%

[ scalar_system_controller , info ] = lcsas0.FindAdaptiveController_AlwaysHidenInfo( P_target );

date_string = datestr(now,'ddmmmyyyy-HHMM');
save(['data/scalar_sys_data-hidden_info-' date_string '.mat' ])

if strcmp(info.Message,'Solved!')

	%% Visualizing %%

	results.SimulationData = [];

	figure;

	hold on;
	plot(P_target*Polyhedron('lb',TimeHorizon,'ub',TimeHorizon),'color','red','alpha',0.5)
% 	plot(P_targets(2),'color','magenta','alpha',0.5)

	for simulation_index = 1:10
		[ x_0_t, u_0_tm1 , y_0_t , sig ] = scalar_system_controller.simulate_1run();

		% for t = 0:TimeHorizon
		% 	x_t = x_0_t(:,t+1);
		% 	scatter(x_t(1),x_t(2))
		% end
		plot([0:length(x_0_t)-1],x_0_t)
		

		%Save data
		results.SimulationData = [ results.SimulationData ; struct('x_0_t',x_0_t,'u_0_tm1',u_0_tm1) ];

	end

	axis([-0.5,2.5,-0.5,2.5])
	xlabel('$z$','Interpreter','latex')
	ylabel('$\dot{z}$','Interpreter','latex')

	saveas(gcf,'images/kltl-scalar_runs-hidden_info','epsc')
	saveas(gcf,'images/kltl-scalar_runs-hidden_info','png')

	%Save controller data somewhere.
	drone_controller.deconstruct_and_save_to(['data/scalar_controller_data_' date+stromg '.mat'])

else
	disp('The drone controller was not successfully created.')
end

save(['data/scalar_sys_data-hidden_info-' date_string '.mat' ])	
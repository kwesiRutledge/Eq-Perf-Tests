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

% Plot the reachable sets on different figures.

[pH,R] = toy1_controller.plotReachableSets('PlotTarget',true,'TargetSet',P_target);

% Plot the reachable sets on the same figure.

num_sequences = size(toy1_controller.KnowledgeSequences,2);
RColors = {'cyan','magenta'};

[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

figure;

for sequence_index = 1:num_sequences
	sp = subplot(num_sequences,1,sequence_index)
	hold on;

	pH_i = [];

	%Plot initial condition
	pH_i(end+1) = scatter(x0(1),x0(2));

	%Plot Target
	pH_i(end+1) = plot(P_target, 'color', 'yellow' , 'alpha' , 0.2);

	for t = 1:TimeHorizon
		% Plot The Polyhedron For Word 1
		if t==1
			pH_i(end+1) = plot( R{sequence_index}.projection(n_x*t+[1:n_x]), ...
				'color',RColors{sequence_index});
		else
			plot( R{sequence_index}.projection(n_x*t+[1:n_x]), ...
				'color',RColors{sequence_index} );
		end
	end

	axis([-3,4.5,-3,3])

	%Create legend
	legend(sp,'$$x_0$$','$$\mathcal{X}_T$$',['$$R(\mathbf{m}^{(' num2str(sequence_index)' ')})$$'],'Interpreter','latex');

	%Create title
	%title(['Reachable Sets for Mode Sequence #' num2str(sequence_index) ])
end

saveas(gcf,'images/toy1_reachable_sets','epsc')
saveas(gcf,'images/toy1_reachable_sets','png')

save(['data/toy1_data_' datestr(now,'ddmmmyyyy-HHMM') '.mat' ])
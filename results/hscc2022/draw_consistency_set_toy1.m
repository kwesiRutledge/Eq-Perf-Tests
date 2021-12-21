%draw_consistency_set_example3.m
%Description:
%	Drawing the consistency sets / reachable sets for the opposing rotation
%	example.
%

%%%%%%%%%%%%%%%%%%%%%%
%% Input Processing %%
%%%%%%%%%%%%%%%%%%%%%%

%% Constants

TimeHorizon = 4;
[ lcsas0 , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_opposing_rotation_lcsas('TimeHorizon',TimeHorizon);
[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

% Defaults
if ~exist('save_gifs')
	save_gifs = true;
end

settings = struct( ...
    'verbosity', 1 , ...
    'DoOptimizationPruningWhere', 'BeforeSearch' ...
);

%% Create Feasible Beliefs %%

[ unchecked_knowledge_sequences , sequence_construction_history ] = lcsas0.L.create_belief_sequences_of_length(TimeHorizon);
[ possible_subsets_of_paths , ~ , choices_as_binary_flags ] = lcsas0.get_feasible_combinations_of_beliefs( unchecked_knowledge_sequences , ...
																												'verbosity' , settings.verbosity , ... 
																												'SkipOptimizationBasedPruning', strcmp(settings.DoOptimizationPruningWhere,'DuringSearch') );

num_beliefs = size(unchecked_knowledge_sequences,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Reachable Sets using EBS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KnowledgeSequences = unchecked_knowledge_sequences;

IBSCollection={};
IBSCollectionAsPolyhedra = {};
for sequence_index = 1:num_beliefs
	temp_knowl_sequence = KnowledgeSequences(:,sequence_index);

	IBSCollectionOL{sequence_index} = InternalBehaviorSet( lcsas0, temp_knowl_sequence, ...
														'fb_type','state');
	IBSCollection{sequence_index} = InternalBehaviorSet( lcsas0,temp_knowl_sequence, ...
														'fb_type','state', ...
														'OpenLoopOrClosedLoop','Closed' , zeros(n_u*TimeHorizon,n_w*TimeHorizon) , zeros(n_u*TimeHorizon,1));
	
	IBSCollectionAsPolyhedra{sequence_index} = IBSCollection{sequence_index}.ToPolyhedron();
	IBSCollectionOLAsPolyhedra{sequence_index} = IBSCollectionOL{sequence_index}.ToPolyhedron();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Comparison of Open Loop Consistency Sets with Closed Loop Consistency Sets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathIndiciesToVisualize = [1,2];

colorArray = {'cyan','magenta'};

primaryAxis = [-2 5 -2 2];
secondaryAxis = [-2 5 -4 4];

hs = [];
figure;
hold on;

scatter(x0(1),x0(2)); %Plot x0

%plot target
hs(end+1) = plot(P_target,'color','yellow','alpha',0.5);

for t = [TimeHorizon-1:-1:1]

	for path_index = pathIndiciesToVisualize

		% Plot The Polyhedron For Word 1
		hs(end+1) = plot( 	IBSCollectionOLAsPolyhedra{path_index}.projection(n_x*t+[1:n_x]), ...
							'color',colorArray{path_index} );

	end

end


axis(secondaryAxis)
legend(hs,'Target','System 1','System 2')

saveas(gcf,['images/toy1_ol_at_time_' num2str(t)],'epsc')
saveas(gcf,['images/toy1_ol_at_time_' num2str(t)],'png')

hs = [];

figure;
hold on;

scatter(x0(1),x0(2)); %Plot x0
%plot target
hs(end+1) = plot(P_target,'color','yellow','alpha',0.5);

for t = [TimeHorizon-1:-1:1]

	for path_index = pathIndiciesToVisualize

		% Plot The Polyhedron For Word 1
		hs(end+1) = plot( IBSCollectionAsPolyhedra{path_index}.projection(n_x*t+[1:n_x]), ...
			'color',colorArray{path_index});

	end

end

axis(primaryAxis)
legend(hs,'Target','System 1','System 2')

saveas(gcf,['images/toy1_zero_input_at_time_' num2str(t)],'epsc')
saveas(gcf,['images/toy1_zero_input_at_time_' num2str(t)],'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Figure Showing Closed Loop Separation of Consistency Sets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathIndiciesToVisualize = [1,2];

colorArray = {'cyan','magenta','cyan','green'};
CLColorArray = {'blue','lightgreen','blue','lightgreen'};

L = lcsas0.L;

PwT = 1;
for t = 1:TimeHorizon
	PwT = PwT * lcsas0.Dyn(1).P_w;
end

% Use That Gain
[S_w_i,S_u_i,~,J_i,f_bar_i] = get_mpc_matrices(lcsas0,'word',L.words{1});

%K = pob1.F_set{matching_seq_id}; k = pob1.u0_set{matching_seq_id};
K = zeros(n_u*TimeHorizon,n_w*TimeHorizon); k = zeros(n_u*TimeHorizon,1);
K_trimmed = K([1:n_u*TimeHorizon],[1:n_w*TimeHorizon]);
k_trimmed = k([1:n_u*TimeHorizon],1);


S_w_trimmed = S_w_i([1:n_x*(TimeHorizon+1)],[1:n_w*TimeHorizon]);
S_u_trimmed = S_u_i([1:n_x*(TimeHorizon+1)],[1:n_u*TimeHorizon]);
J_trimmed = J_i([1:n_x*(TimeHorizon+1)],:);
f_bar_prime = S_w_trimmed*f_bar_i([1:n_w*TimeHorizon]);

x1 = ( S_w_trimmed + S_u_trimmed * K_trimmed)*PwT + S_u_trimmed*k_trimmed + J_trimmed * Px0 + f_bar_prime ;

[S_w_i,S_u_i,~,J_i,f_bar_i] = get_mpc_matrices(lcsas0,'word',L.words{2});

%K = pob1.F_set{matching_seq_id}; k = pob1.u0_set{matching_seq_id};
K = zeros(n_u*TimeHorizon,n_w*TimeHorizon); k = zeros(n_u*TimeHorizon,1);
K_trimmed = K([1:n_u*t],[1:n_w*TimeHorizon]);
k_trimmed = k([1:n_u*t],1);


S_w_trimmed = S_w_i([1:n_x*(t+1)],[1:n_w*TimeHorizon]);
S_u_trimmed = S_u_i([1:n_x*(t+1)],[1:n_u*t]);
J_trimmed = J_i([1:n_x*(t+1)],:);
f_bar_prime = S_w_trimmed*f_bar_i([1:n_w*TimeHorizon]);

x2 = ( S_w_trimmed + S_u_trimmed * K_trimmed)*PwT + S_u_trimmed*k_trimmed + J_trimmed * Px0 + f_bar_prime ;

% ReachableSet = x1.projection(n_x*t+[1:n_x]);

x = {x1,x2};

hs = [];

figure;
hold on;

scatter(x0(1),x0(2)); %Plot x0

for t = [TimeHorizon:-1:1]

	for path_index = pathIndiciesToVisualize

		% Plot The Polyhedron For Word 1
		if length(hs) < 2
			hs(end+1) = plot( x{path_index}.projection(n_x*t+[1:n_x]), ...
				'color',colorArray{path_index});
		else
			plot( x{path_index}.projection(n_x*t+[1:n_x]), ...
				'color',colorArray{path_index} );
		end

	end

end

%plot target
hs(end+1) = plot(P_target,'color','yellow','alpha',0.5);

axis(primaryAxis)
legend(hs,'System 1','System 2','Target')

saveas(gcf,['images/toy1_zero_input_complete_at_time_' num2str(t)],'epsc')
saveas(gcf,['images/toy1_zero_input_complete_at_time_' num2str(t)],'png')

for path_index = pathIndiciesToVisualize

	hs = [];

	figure;
	hold on;

	scatter(x0(1),x0(2)); %Plot x0

	for t = [TimeHorizon:-1:1]

		% Plot The Polyhedron For Word 1
		if length(hs) < 1
			hs(end+1) = plot( x{path_index}.projection(n_x*t+[1:n_x]), ...
				'color',colorArray{path_index});
		else
			plot( x{path_index}.projection(n_x*t+[1:n_x]), ...
				'color',colorArray{path_index} );
		end

	end

	%plot target
	hs(end+1) = plot(P_target,'color','yellow','alpha',0.5);

	axis(primaryAxis)
	legend(hs,['System ' num2str(path_index) ],'Target')

	saveas(gcf,['images/toy1_zero_input_mode' num2str(path_index) '_complete_at_time_' num2str(t)],'epsc')
	saveas(gcf,['images/toy1_zero_input_mode' num2str(path_index) '_complete_at_time_' num2str(t)],'png')

end

%% Create this as one plot

figure
for path_index = pathIndiciesToVisualize

	hs = [];

	subplot(2,1,path_index)
	hold on;

	scatter(x0(1),x0(2)); %Plot x0

	for t = [TimeHorizon:-1:1]

		% Plot The Polyhedron For Word 1
		if length(hs) < 1
			hs(end+1) = plot( x{path_index}.projection(n_x*t+[1:n_x]), ...
				'color',colorArray{path_index});
		else
			plot( x{path_index}.projection(n_x*t+[1:n_x]), ...
				'color',colorArray{path_index} );
		end

	end

	%plot target
	hs(end+1) = plot(P_target,'color','yellow','alpha',0.5);

	axis(primaryAxis)
	legend(hs,['System ' num2str(path_index) ],'Target')

end

saveas(gcf,['images/toy1_zero_input_v2_complete_at_time_' num2str(t)],'epsc')
saveas(gcf,['images/toy1_zero_input_v2_complete_at_time_' num2str(t)],'png')

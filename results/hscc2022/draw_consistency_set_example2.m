%draw_consistency_set_example2.m
%Description:
%	Based on observer_comparison94.m

%%%%%%%%%%%%%%%%%%%%%%
%% Input Processing %%
%%%%%%%%%%%%%%%%%%%%%%

%% Constants

TimeHorizon = 4;
[ lcsas0 , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_similar_rotation_lcsas('TimeHorizon',TimeHorizon);
[ n_x , n_u , n_y , n_w , n_v ] = lcsas0.Dimensions();

% Defaults
if ~exist('save_gifs')
	save_gifs = true;
end

results.Parameters.LCSAS = lcsas0;

load('data/toy2_data_21Aug2021-1457.mat','toy2_controller')
KnowledgeSequences = toy2_controller.KnowledgeSequences;
num_beliefs = size(toy2_controller.KnowledgeSequences,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Reachable Sets using EBS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IBSCollection={};
IBSCollectionAsPolyhedra = {};
for sequence_index = 1:num_beliefs
	temp_knowl_sequence = KnowledgeSequences(:,sequence_index);

	IBSCollectionOL{sequence_index} = InternalBehaviorSet( lcsas0,temp_knowl_sequence, ...
														'fb_type','state');
	IBSCollection{sequence_index} = InternalBehaviorSet( lcsas0,temp_knowl_sequence, ...
														'fb_type','state', ...
														'OpenLoopOrClosedLoop','Closed' , toy2_controller.K_set{sequence_index} , toy2_controller.k_set{sequence_index});
	
	IBSCollectionAsPolyhedra{sequence_index} = IBSCollection{sequence_index}.ToPolyhedron();
	IBSCollectionOLAsPolyhedra{sequence_index} = IBSCollectionOL{sequence_index}.ToPolyhedron();
end

figure;
hold on;

scatter(x0(1),x0(2)); %Plot x0

for t = 1:TimeHorizon-1

	plot( IBSCollectionAsPolyhedra{1}.projection(n_x*t+[1:n_x]), ...
		'color','salmon' )

	% Plot The Polyhedron For Word 2
	plot( IBSCollectionAsPolyhedra{2}.projection(n_x*t+[1:n_x]), ...
		'color','cyan')

end

% axis([-0.5 11 -0.5 6.5 ])

saveas(gcf,'images/similarRotationSystemCLImage','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Comparison of Open Loop Consistency Sets with Closed Loop Consistency Sets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathIndiciesToVisualize = [3,4];

for path_index = pathIndiciesToVisualize

	figure;
	hold on;

	scatter(x0(1),x0(2)); %Plot x0

	for t = 1:TimeHorizon-1

		% Plot The Polyhedron For Word 2
		plot( IBSCollectionOLAsPolyhedra{path_index}.projection(n_x*t+[1:n_x]), ...
			'color','cyan')
	end

	for t = 1:TimeHorizon-1

		plot( IBSCollectionAsPolyhedra{path_index}.projection(n_x*t+[1:n_x]), ...
			'color','blue' )

	end

	saveas(gcf,['images/example2_olVScl_over_all_t'],'epsc')
	saveas(gcf,['images/example2_olVScl_over_all_t'],'png')

	for t = 1:TimeHorizon-1
		figure;
		hold on;

		% Plot The Polyhedron For Word 2
		plot( IBSCollectionOLAsPolyhedra{path_index}.projection(n_x*t+[1:n_x]), ...
			'color','cyan')

		plot( IBSCollectionAsPolyhedra{path_index}.projection(n_x*t+[1:n_x]), ...
			'color','blue' )

		legend('$\mathcal{C}_{S}(p)$','$\mathcal{C}_{S}^{(cl)}(p)$','Interpreter','latex')

		saveas(gcf,['images/example2_olVScl_at_time_' num2str(t)],'epsc')
		saveas(gcf,['images/example2_olVScl_at_time_' num2str(t)],'png')

	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Figure Showing Closed Loop Separation of Consistency Sets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathIndiciesToVisualize = [3,4];

colorArray = {'cyan','green','cyan','green'};
CLColorArray = {'blue','lightgreen','blue','lightgreen'};

t1 = 2;

figure;
hold on;

for path_index = pathIndiciesToVisualize

	% Plot The Polyhedron For Word 2
	plot( IBSCollectionOLAsPolyhedra{path_index}.projection(n_x*t+[1:n_x]), ...
		'color',colorArray{path_index})

end

for path_index = pathIndiciesToVisualize

	plot( IBSCollectionAsPolyhedra{path_index}.projection(n_x*t+[1:n_x]), ...
			'color',CLColorArray{path_index} )

end

saveas(gcf,['images/example2_two_olVScl_at_time_' num2str(t1)],'epsc')
saveas(gcf,['images/example2_two_olVScl_at_time_' num2str(t1)],'png')
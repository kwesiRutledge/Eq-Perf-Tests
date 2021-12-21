%draw_summary_figure.m
%Description:
%	Drawing the consistency sets / reachable sets for the opposing rotation
%	example with saved data.
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

load('data/toy1_data_30Sep2021-2131.mat','toy1_controller')

KnowledgeSequences = unchecked_knowledge_sequences;

IBSCollection={};
IBSCollectionAsPolyhedra = {};
for sequence_index = 1:2
	temp_knowl_sequence = KnowledgeSequences(:,sequence_index);

	IBSCollectionOL{sequence_index} = InternalBehaviorSet( lcsas0, temp_knowl_sequence, ...
														'fb_type','state');
	IBSCollection{sequence_index} = InternalBehaviorSet( lcsas0,temp_knowl_sequence, ...
														'fb_type','state', ...
														'OpenLoopOrClosedLoop','Closed' , toy1_controller.K_set{sequence_index} , toy1_controller.k_set{sequence_index});
	
	IBSCollectionAsPolyhedra{sequence_index} = IBSCollection{sequence_index}.ToPolyhedron();
	IBSCollectionOLAsPolyhedra{sequence_index} = IBSCollectionOL{sequence_index}.ToPolyhedron();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Figure Showing Comparison of Closed Loop And Open Loop of Consistency Sets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathIndiciesToVisualize = [1,2];

colorArray = {'cyan','magenta','cyan','green'};
CLColorArray = {'blue','lightgreen','blue','lightgreen'};

L = lcsas0.L;

PwT = 1; PuT = 1;
for t = 1:TimeHorizon
	PwT = PwT * lcsas0.Dyn(1).P_w;
	PuT = PuT * lcsas0.U;
end

% Use That Gain
[S_w_i,S_u_i,~,J_i,f_bar_i] = get_mpc_matrices(lcsas0,'word',L.words{1});

%K = pob1.F_set{matching_seq_id}; k = pob1.u0_set{matching_seq_id};
K = toy1_controller.K_set{1}; k = toy1_controller.k_set{1};
K_trimmed = K([1:n_u*TimeHorizon],[1:n_w*TimeHorizon]);
k_trimmed = k([1:n_u*TimeHorizon],1);


S_w_trimmed = S_w_i([1:n_x*(TimeHorizon+1)],[1:n_w*TimeHorizon]);
S_u_trimmed = S_u_i([1:n_x*(TimeHorizon+1)],[1:n_u*TimeHorizon]);
J_trimmed = J_i([1:n_x*(TimeHorizon+1)],:);
f_bar_prime = S_w_trimmed*f_bar_i([1:n_w*TimeHorizon]);

x1 = ( S_w_trimmed + S_u_trimmed * K_trimmed)*PwT + S_u_trimmed*k_trimmed + J_trimmed * Px0 + f_bar_prime ;

[S_w_i,S_u_i,~,J_i,f_bar_i] = get_mpc_matrices(lcsas0,'word',L.words{2});

%K = pob1.F_set{matching_seq_id}; k = pob1.u0_set{matching_seq_id};
K = toy1_controller.K_set{2}; k = toy1_controller.k_set{2};
K_trimmed = K([1:n_u*t],[1:n_w*TimeHorizon]);
k_trimmed = k([1:n_u*t],1);


S_w_trimmed = S_w_i([1:n_x*(t+1)],[1:n_w*TimeHorizon]);
S_u_trimmed = S_u_i([1:n_x*(t+1)],[1:n_u*t]);
J_trimmed = J_i([1:n_x*(t+1)],:);
f_bar_prime = S_w_trimmed*f_bar_i([1:n_w*TimeHorizon]);

x2 = ( S_w_trimmed + S_u_trimmed * K_trimmed)*PwT + S_u_trimmed*k_trimmed + J_trimmed * Px0 + f_bar_prime ;

% ReachableSet = x1.projection(n_x*t+[1:n_x]);

xR = S_w_trimmed * PwT + S_u_trimmed * PuT + J_trimmed * Px0 + f_bar_prime ;

x = {x1,x2};

hs = [];
fs = 15;

figure('DefaultAxesFontSize',fs);
hold on;

scatter(x0(1),x0(2)); %Plot x0

%Plot target set
hs(end+1) =  plot(P_target, 'color','yellow','Alpha',0.5);

% Plot total reachable set
hs(end+1) = plot(xR.projection(n_x*TimeHorizon+[1:n_x]), ...
			'color','cyan','Alpha',0.15);

for t = [TimeHorizon:-1:1]

	% Plot The Polyhedron For Word 1
	if length(hs) < 3
		hs(end+1) = plot( x{1}.projection(n_x*t+[1:n_x]), ...
			'color','cyan','LineStyle',':');
		hs(end+1) = plot( x{2}.projection(n_x*t+[1:n_x]), ...
			'color','magenta','LineStyle',':');
		
	else

		plot( x{1}.projection(n_x*t+[1:n_x]), ...
			'color','cyan','LineStyle',':' );
		plot( x{2}.projection(n_x*t+[1:n_x]), ...
			'color','magenta','LineStyle',':');
	end

end

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')

%axis(secondaryAxis)
legend(hs,'Target','Reachable Set','$$\mathcal{C}(\mathbf{m}^{(1)})$$','$$\mathcal{C}(\mathbf{m}^{(2)})$$','Interpreter','latex')

saveas(gcf,['images/summary_figure'],'epsc')
saveas(gcf,['images/summary_figure'],'png')
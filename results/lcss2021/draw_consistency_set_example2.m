%draw_consistency_set_example2.m
%Description:
%	Based on observer_comparison94.m

%%%%%%%%%%%%%%%%%%%%%%
%% Input Processing %%
%%%%%%%%%%%%%%%%%%%%%%

%% Constants

TimeHorizon = 4;
[ lcsas0 , TimeHorizon , Pu , Pw , x0 , Px0 , P_target ] = get_similar_rotation_lcsas('TimeHorizon',TimeHorizon);

% Defaults
if ~exist('save_gifs')
	save_gifs = true;
end

results.Parameters.LCSAS = lcsas0;

load('data/toy2_data_20Aug2021-1810.mat','toy2_controller')
KnowledgeSequences = toy2_controller.KnowledgeSequences;
num_beliefs = size(toy2_controller.KnowledgeSequences,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Reachable Sets using EBS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IBSCollection={};
IBSCollectionAsPolyhedra = {};
for sequence_index = 1:num_beliefs
	temp_knowl_sequence = KnowledgeSequences(:,sequence_index);

	IBSCollection{sequence_index} = InternalBehaviorSet( lcsas0,temp_knowl_sequence, ...
														'fb_type','state', ...
														'OpenLoopOrClosedLoop','Closed' , toy2_controller.K_set{sequence_index} , toy2_controller.k_set{sequence_index});
	sequence_index
	toy2_controller.K_set{sequence_index} 
	toy2_controller.k_set{sequence_index}
	IBSCollectionAsPolyhedra{sequence_index} = IBSCollection{sequence_index}.ToPolyhedron();
end

figure;
hold on;

scatter(x0(1),x0(2)); %Plot x0

for t = 1:TimeHorizon-1

	plot( IBSCollectionAsPolyhedra{1}.projection(dim_x*t+[1:dim_x]), ...
		'color','salmon' )

	% Plot The Polyhedron For Word 2
	plot( IBSCollectionAsPolyhedra{2}.projection(dim_x*t+[1:dim_x]), ...
		'color','cyan')

end

% axis([-0.5 11 -0.5 6.5 ])

saveas(gcf,'images/similarRotationSystemCLImage','epsc')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Compare Closed Loop Consistency Sets %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ConsistencySets2 = {};
% 
% ad1 = Aff_Dyn( A1-0.2*eye(2) , zeros(size(B1)) , K1 , eye(dim_x) , Pw , Pv );
% ad2 = Aff_Dyn( A2-0.2*eye(2) , zeros(size(B2)) , K2 , eye(dim_x) , Pw , Pv );
% 
% lcsas1 = LCSAS( [ad1,ad2], Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) , ...
% 				'X0' , Polyhedron('lb',x0','ub',x0') , ...
% 				'U'  , Pu  );
% 
% Px0 = Polyhedron('lb',x0','ub',x0');
% 
% for word_index = 1:lcsas1.L.cardinality()
% 	temp_word = lcsas1.L.words{word_index};
% 	temp_single_word_lang = Language(temp_word);
% 	temp_knowl_seq = [lcsas1.L; repmat(temp_single_word_lang,TimeHorizon-1,1)];
% 
% 	ConsistencySets2{word_index} = ExternalBehaviorSet(lcsas1,temp_knowl_seq,'fb_type','state');
% 	CSAsPolyhedron2{word_index} = ConsistencySets2{word_index}.ToPolyhedron();
% end
% 
% results.ConsistencySets2 = ConsistencySets2;
% 
% figure;
% hold on;
% 
% scatter(x0(1),x0(2)); %Plot x0
% 
% for t = 1:TimeHorizon-1
% 
% 	plot( CSAsPolyhedron2{1}.projection(dim_x*t+[1:dim_x]), ...
% 		'color','salmon' )
% 
% 	% Plot The Polyhedron For Word 2
% 	plot( CSAsPolyhedron2{2}.projection(dim_x*t+[1:dim_x]), ...
% 		'color','cyan')
% 
% end
% 
% axis([-0.5 11 -0.5 6.5 ])
% 
% saveas(gcf,'images/similarRotationSystemComparison2','epsc')
% 
% return;
% 
% %Test Where intersections Lie
% intersectionIsNonempty = [];
% for t = 1:TimeHorizon-1
% 	temp_intersection = CSets2{1,t}.intersect(CSets2{2,t});
% 	intersectionIsNonempty = [ intersectionIsNonempty , ~temp_intersection.isEmptySet ];
% end
% 
% results.IntersectionNonemptyVector = intersectionIsNonempty;
% 
% %%%%%%%%%%%%%
% %% Results %%
% %%%%%%%%%%%%%
% 
% function saveToSimpleGIF( TimeHorizonIn , PolyX_History, x0 , X_Target , gifFilename , axis_limits , colorIn , titleIn )
% 	%Description:
% 	%	Saves simple version of the input Polytope to a gif file.
% 
% 	%% Constants
% 
% 	secondsPerImage = 1;
% 
% 	if ~exist('titleIn')
% 		titleIn = 'Temporary GIF Title'
% 	end
% 
% 	dim_x = X_Target.Dim;
% 
% 	%% Algorithm
% 
% 	h = figure;
% 
% 	%filename = 'testAnimated.gif';
% 	for t = 0:TimeHorizonIn
% 
% 	    % Draw plot for current PolyX_History
% 	    hold on;
% 	    plot(X_Target,'Color','White') %Plot Target Set
% 	    if t == 0
% 	    	scatter(x0(1),x0(2))
% 	    else
% 	    	plot(PolyX_History{t}.projection([PolyX_History{t}.Dim-dim_x+1:PolyX_History{t}.Dim]),'color',colorIn)
% 	    end
% 	    hold off;
% 
% 	    title(titleIn)
% 		grid on
% 		axis(axis_limits)
% 
% 	    drawnow 
% 	    
% 	    % Capture the plot as an image 
% 	    frame = getframe(h); 
% 	    im = frame2im(frame); 
% 	    [imind,cm] = rgb2ind(im,256); 
% 	    % Write to the GIF File 
% 	    if t == 0 
% 	        imwrite(imind,cm,gifFilename,'gif', 'Loopcount',inf,'DelayTime',secondsPerImage); 
% 	    else 
% 	        imwrite(imind,cm,gifFilename,'gif','WriteMode','append'); 
% 	    end 
% 
% 	end
% end
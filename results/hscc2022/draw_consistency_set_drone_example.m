%draw_consistency_set_drone_example.m
%Description:
%	Based on observer_comparison94.m

%%%%%%%%%%%%%%%%%%%%%%
%% Input Processing %%
%%%%%%%%%%%%%%%%%%%%%%


%% Constants

TimeHorizon = 5; %x0 = [1;0]; X_target = Polyhedron('lb',[3],'ub',[9]) * Polyhedron('lb',-10,'ub',10);
[ lcsas0 , x0 , TimeHorizon , P_target ] = get_differently_loaded_drone_lcsas(	'TimeHorizon',TimeHorizon, ...
																				'm1',1.0,'m2',1.5 );

axes_lims = [-0.5 11 -3 11 ];

% Defaults
if ~exist('save_gifs')
	save_gifs = true;
end

results.Parameters.LCSAS = lcsas0;

[ dim_x , dim_u , dim_y , dim_w , dim_v ] = lcsas0.Dimensions();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Reachable Sets using EBS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConsistencySets={};
for word_index = 1:lcsas0.L.cardinality()
	temp_word = lcsas0.L.words{word_index};
	temp_single_word_lang = Language(temp_word);
	temp_knowl_seq = [lcsas0.L; repmat(temp_single_word_lang,TimeHorizon-1,1)];

	ConsistencySets{word_index} = ExternalBehaviorSet(lcsas0,temp_knowl_seq,'fb_type','state');
	CSAsPolyhedron{word_index} = ConsistencySets{word_index}.ToPolyhedron();
end

time_reversed = fliplr([1:TimeHorizon-1]);

figure;
hold on;

scatter(x0(1),x0(2)); %Plot x0

for t = 1:TimeHorizon-1

	plot( CSAsPolyhedron{1}.projection(dim_x*t+[1:dim_x]), ...
		'color','salmon' )

	% Plot The Polyhedron For Word 2
	plot( CSAsPolyhedron{2}.projection(dim_x*t+[1:dim_x]), ...
		'color','cyan')

end

axis(axes_lims)

xlabel('$(x_k)_1$','Interpreter','latex','FontSize',20)
ylabel('$(x_k)_2$','Interpreter','latex','FontSize',20)

saveas(gcf,'images/similarRotationSystemComparison1','epsc')

%% Draw GIF %%
%%%%%%%%%%%%%%
eta_w = 0.25;
Pw = lcsas0.Dyn(1).P_w;

% Create PwT;
PwT = {};
for t = 1:TimeHorizon
	PwT{t} = 1;
	for tau = 1:t
		PwT{t} = PwT{t} * Pw;
	end
end

% Create PuT
eta_u = 2*eta_w;
Pu = Polyhedron(...
	'lb',-eta_u*ones(1,dim_u), ...
	'ub', eta_u*ones(1,dim_u));

PuT = {};
for t = 1:TimeHorizon
	PuT{t} = 1;
	for tau = 1:t
		PuT{t} = PuT{t} * Pu;
	end
end

for t = 1:TimeHorizon
	%Plot The Polyhedron.
    word1_prefix = lcsas0.L.words{1}(1:t);
	[Sw,Su,~,J,~] = lcsas0.get_mpc_matrices('word',word1_prefix);
	Pxt{1,t} = Sw * PwT{t} + Su * PuT{t} + J * x0;
	plot( ...
		Pxt{1,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
		'color','salmon' ...
		)

	% Plot The Polyhedron For Word 2
	word2_prefix = lcsas0.L.words{2}(1:t);
	[Sw,Su,~,J,~] = lcsas0.get_mpc_matrices('word',word2_prefix);
	Pxt{2,t} = Sw * PwT{t} + Su * PuT{t} + J * x0;
	plot( ...
		Pxt{2,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
		'color','cyan')

end

% Save GIF for each mode
mode_colors = {'magenta','cyan'};
for mode_index = 1:lcsas0.n_modes
	saveToSimpleGIF( ...
		TimeHorizon , {Pxt{mode_index,:}}, x0 , X_Target , ...
		['images/droneSystem_Mode' num2str(mode_index) '_ReachableSet.gif'] , ...
		temp_axis , ...
		mode_colors{mode_index} , ...
		['Mode ' num2str(mode_index) ' Reachable Set'] )
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare Closed Loop Consistency Sets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConsistencySets2 = {};

ad1 = Aff_Dyn( A1-0.2*eye(2) , zeros(size(B1)) , K1 , eye(dim_x) , Pw , Pv );
ad2 = Aff_Dyn( A2-0.2*eye(2) , zeros(size(B2)) , K2 , eye(dim_x) , Pw , Pv );

lcsas1 = LCSAS( [ad1,ad2], Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) , ...
				'X0' , Polyhedron('lb',x0','ub',x0') , ...
				'U'  , Pu  );

Px0 = Polyhedron('lb',x0','ub',x0');

for word_index = 1:lcsas1.L.cardinality()
	temp_word = lcsas1.L.words{word_index};
	temp_single_word_lang = Language(temp_word);
	temp_knowl_seq = [lcsas1.L; repmat(temp_single_word_lang,TimeHorizon-1,1)];

	ConsistencySets2{word_index} = ExternalBehaviorSet(lcsas1,temp_knowl_seq,'fb_type','state');
	CSAsPolyhedron2{word_index} = ConsistencySets2{word_index}.ToPolyhedron();
end

results.ConsistencySets2 = ConsistencySets2;

figure;
hold on;

scatter(x0(1),x0(2)); %Plot x0

for t = 1:TimeHorizon-1

	plot( CSAsPolyhedron2{1}.projection(dim_x*t+[1:dim_x]), ...
		'color','salmon' )

	% Plot The Polyhedron For Word 2
	plot( CSAsPolyhedron2{2}.projection(dim_x*t+[1:dim_x]), ...
		'color','cyan')

end

axis(axes_lims)

xlabel('$(x_k)_1$','Interpreter','latex')
ylabel('$(x_k)_2$','Interpreter','latex')

saveas(gcf,'images/similarRotationSystemComparison2','epsc')

return;

%Test Where intersections Lie
intersectionIsNonempty = [];
for t = 1:TimeHorizon-1
	temp_intersection = CSets2{1,t}.intersect(CSets2{2,t});
	intersectionIsNonempty = [ intersectionIsNonempty , ~temp_intersection.isEmptySet ];
end

results.IntersectionNonemptyVector = intersectionIsNonempty;

%%%%%%%%%%%%%
%% Results %%
%%%%%%%%%%%%%

function saveToSimpleGIF( TimeHorizonIn , PolyX_History, x0 , X_Target , gifFilename , axis_limits , colorIn , titleIn )
	%Description:
	%	Saves simple version of the input Polytope to a gif file.

	%% Constants

	secondsPerImage = 1;

	if ~exist('titleIn')
		titleIn = 'Temporary GIF Title'
	end

	dim_x = X_Target.Dim;

	%% Algorithm

	h = figure;

	%filename = 'testAnimated.gif';
	for t = 0:TimeHorizonIn

	    % Draw plot for current PolyX_History
	    plot(X_Target,'Color','White') %Plot Target Set
	    hold on;
	    if t == 0
	    	scatter(x0(1),x0(2))
	    else
	    	plot(PolyX_History{t}.projection([PolyX_History{t}.Dim-dim_x+1:PolyX_History{t}.Dim]),'color',colorIn)
	    end
	    hold off;

	    legend('Target','$\mathcal{R}(\sigma)$', ...
    			'Interpreter','latex')

	    title(titleIn)
	    xlabel('$(x_k)_1$','Interpreter','latex','FontSize',20)
		ylabel('$(x_k)_2$','Interpreter','latex','FontSize',20)

		grid on
		axis(axis_limits)

	    drawnow 
	    
	    % Capture the plot as an image 
	    frame = getframe(h); 
	    im = frame2im(frame); 
	    [imind,cm] = rgb2ind(im,256); 
	    % Write to the GIF File 
	    if t == 0 
	        imwrite(imind,cm,gifFilename,'gif', 'Loopcount',inf,'DelayTime',secondsPerImage); 
	    else 
	        imwrite(imind,cm,gifFilename,'gif','WriteMode','append'); 
	    end 

	end
end
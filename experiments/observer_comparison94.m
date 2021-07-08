function [results] = observer_comparison94(varargin)
	%observer_comparison94.m
	%Description:
	%	Second Illustrative Example.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargin >= 1
		save_gifs = varargin{1};
	end


	%% Constants

	twoDRotation = @(theta) [ cos(theta), -sin(theta) ; sin(theta), cos(theta) ];
	dim_x = 2;

	r1 = 10;
	r2 = 11;

	A1 = twoDRotation(pi/12);
	A2 = twoDRotation(pi/12);

	B1 = eye(dim_x);
	B2 = eye(dim_x);

	K1 = -A1*[ 0 ; r1 ]+[0;r1];
	K2 = -A2*[ 0 ; r2 ]+[0;r2];

	eta_w = 0.5;
	Pw = Polyhedron('lb',-eta_w*ones(1,dim_x),'ub',eta_w*ones(1,dim_x));
	Pv = Pw; %We won't use it.

	TimeHorizon = 4;

	% Create PwT;
    PwT = {};
	for t = 1:TimeHorizon
		PwT{t} = 1;
		for tau = 1:t
			PwT{t} = PwT{t} * Pw;
		end
	end

	% Create PuT
	eta_u = 0.5*eta_w;
	Pu = Polyhedron(...
		'lb',-eta_u*ones(1,dim_x), ...
		'ub', eta_u*ones(1,dim_x))
	PuT = {};
    for t = 1:TimeHorizon
		PuT{t} = 1;
		for tau = 1:t
			PuT{t} = PuT{t} * Pu;
		end
	end

	% Create XT
	X_Target = Polyhedron('lb',-eta_w*TimeHorizon*ones(1,dim_x), 'ub', eta_w*TimeHorizon*ones(1,dim_x) ) + [9.5;5];

	% Defaults
	if ~exist('save_gifs')
		save_gifs = true;
	end

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	ad1 = Aff_Dyn( A1 , B1 , K1 , eye(dim_x) , Pw , Pv );
	ad2 = Aff_Dyn( A2 , B2 , K2 , eye(dim_x) , Pw , Pv );

	lcsas0 = LCSAS( [ad1,ad2], Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) );

	results.Parameters.LCSAS = lcsas0;

	%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Show Reachable Sets %%
	%%%%%%%%%%%%%%%%%%%%%%%%%

	x0 = [0;0];
	
	figure;
	hold on;

	scatter(x0(1),x0(2))

	for t = 1:TimeHorizon

		%Plot The Polyhedron.
        word1_prefix = lcsas0.L.words{1}(1:t);
		[Sw,Su,~,J,K_bar] = lcsas0.get_mpc_matrices('word',word1_prefix);
		Pxt{1,t} = Sw * PwT{t} + Su * PuT{t} + J * x0 + Sw*K_bar;
		plot( ...
			Pxt{1,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
			'color','salmon' ...
			)

		% Plot The Polyhedron For Word 2
		word2_prefix = lcsas0.L.words{2}(1:t);
		[Sw,Su,~,J,K_bar] = lcsas0.get_mpc_matrices('word',word2_prefix);
		Pxt{2,t} = Sw * PwT{t} + Su * PuT{t} + J * x0 + Sw*K_bar;
		plot( ...
			Pxt{2,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
			'color','cyan')

    end

    plot(X_Target,'FaceColor','none')

    title('Reachable Sets for Each Mode')
    legend('$x_0$','Mode 1','Mode 2', '', '' , 'Target', ...
    	'Interpreter','latex')

    % Save GIF for each mode
    if save_gifs
	    mode_colors = {'magenta','cyan'};
	    for mode_index = 1:lcsas0.n_modes
	    	saveToSimpleGIF( ...
	    		TimeHorizon , {Pxt{mode_index,:}}, x0 , X_Target , ...
	    		['results/lcss2021/images/similarRotationSystem_Mode' num2str(mode_index) '_ReachableSet.gif'] , ...
	    		[-2,14,-1,10] , ...
	    		mode_colors{mode_index} , ...
	    		['Mode ' num2str(mode_index) ' Reachable Set'] )
	    end
	end

    results.ReachableSetAtTime = Pxt;

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Show Zero Input Sets %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;
	hold on;

	scatter(x0(1),x0(2))

	for t = 1:TimeHorizon

		%Plot The Polyhedron.
        word1_prefix = lcsas0.L.words{1}(1:t);
		[Sw,~,~,J,K_bar] = lcsas0.get_mpc_matrices('word',word1_prefix);
		Pxt2{1,t} = Sw * PwT{t} + J * x0 + Sw*K_bar;
		plot( ...
			Pxt2{1,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
			'color','salmon' ...
			)

		% Plot The Polyhedron For Word 2
		word2_prefix = lcsas0.L.words{2}(1:t);
		[Sw,~,~,J,K_bar] = lcsas0.get_mpc_matrices('word',word2_prefix);
		Pxt2{2,t} = Sw * PwT{t} + J * x0 + Sw*K_bar;
		plot( ...
			Pxt2{2,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
			'color','cyan')

    end

    plot(X_Target,'FaceColor','none')

    title('Zero Input Reachable Sets for Each Mode')
    legend('$x_0$','Mode 1','Mode 2', '', '' , 'Target', ...
    	'Interpreter','latex')

    % Save to Gifs
    if save_gifs
	    for mode_index = 1:lcsas0.n_modes
	    	saveToSimpleGIF( ...
	    		TimeHorizon , {Pxt2{mode_index,:}}, x0 , X_Target , ...
	    		['results/lcss2021/images/similarRotationSystem_Mode' num2str(mode_index) '_ZeroInputSet.gif'] , ...
	    		[-2,14,-1,10] , ...
	    		mode_colors{mode_index} , ...
	    		['Mode ' num2str(mode_index) ' Zero-Input Reachable Set'] )
	    end
	end

    results.ZeroInputTrajectoryTubes = Pxt2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compare Consistency Sets %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    CSets = {};

    Px0 = Polyhedron('lb',x0','ub',x0');

    for mode_index = 1:lcsas0.n_modes
    	for t = 1:TimeHorizon-1
    		[ CSets{mode_index,t} , ~ ] = lcsas0.consistent_set(t,lcsas0.L,Pu,Px0,'fb_method','state');
    	end
    end

    results.CSets = CSets;

    %Test Where intersections Lie
    intersectionIsNonempty = [];
    for t = 1:TimeHorizon-1
    	temp_intersection = CSets{1,t}.intersect(CSets{2,t});
    	intersectionIsNonempty = [ intersectionIsNonempty , ~temp_intersection.isEmptySet ];
    end

    results.IntersectionNonemptyVector = intersectionIsNonempty;

    %%%%%%%%%%%%%
    %% Results %%
    %%%%%%%%%%%%%
   
    
end

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
	    hold on;
	    plot(X_Target,'Color','White') %Plot Target Set
	    if t == 0
	    	scatter(x0(1),x0(2))
	    else
	    	plot(PolyX_History{t}.projection([PolyX_History{t}.Dim-dim_x+1:PolyX_History{t}.Dim]),'color',colorIn)
	    end
	    hold off;

	    title(titleIn)
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
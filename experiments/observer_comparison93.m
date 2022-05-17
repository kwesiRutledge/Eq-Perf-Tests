function [results] = observer_comparison93(varargin)


	%% Constants

	twoDRotation = @(theta) [ cos(theta), -sin(theta) ; sin(theta), cos(theta) ];
	dim_x = 2;

	TimeHorizon = 4;

	A1 = twoDRotation(pi/TimeHorizon);
	A2 = twoDRotation(-pi/TimeHorizon);

	B1 = eye(dim_x);
	B2 = eye(dim_x);

	eta_w = 0.25;
	Pw = Polyhedron('lb',-eta_w*ones(1,dim_x),'ub',eta_w*ones(1,dim_x));
	Pv = Pw; %We won't use it.

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
	X_Target = Polyhedron('lb',-(eta_w*TimeHorizon+sqrt(eta_w))*ones(1,dim_x), 'ub', (eta_w*TimeHorizon+sqrt(eta_w))*ones(1,dim_x) ) + [2.75;0];

	temp_axis = [-3,5,-4,4];

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	ad1 = Aff_Dyn( A1 , B1 , zeros(dim_x,1) , eye(dim_x) , Pw , Pv );
	ad2 = Aff_Dyn( A2 , B2 , zeros(dim_x,1) , eye(dim_x) , Pw , Pv );

	lcsas0 = LCSAS( [ad1,ad2], Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) );

	%% Show Reachable Sets

	x0 = [-1;0];
	
	figure;
	hold on;

	scatter(x0(1),x0(2))

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

    plot(X_Target,'FaceColor','none')

    title('Reachable Sets for Each Mode')
    legend('$x_0$','Mode 1','Mode 2', '', '' , 'Target', ...
    	'Interpreter','latex')

    xlabel('$(x_k)_1$','Interpreter','latex')
	ylabel('$(x_k)_2$','Interpreter','latex')

    % Save GIF for each mode
    mode_colors = {'magenta','cyan'};
    for mode_index = 1:lcsas0.n_modes
    	saveToSimpleGIF( ...
    		TimeHorizon , {Pxt{mode_index,:}}, x0 , X_Target , ...
    		['results/hscc2022/images/pureRotationSystem_Mode' num2str(mode_index) '_ReachableSet.gif'] , ...
    		temp_axis , ...
    		mode_colors{mode_index} , ...
    		['Mode ' num2str(mode_index) ' Reachable Set'] )
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
		[Sw,~,~,J,~] = lcsas0.get_mpc_matrices('word',word1_prefix);
		Pxt2{1,t} = Sw * PwT{t} + J * x0;
		plot( ...
			Pxt2{1,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
			'color','salmon' ...
			)

		% Plot The Polyhedron For Word 2
		word2_prefix = lcsas0.L.words{2}(1:t);
		[Sw,~,~,J,~] = lcsas0.get_mpc_matrices('word',word2_prefix);
		Pxt2{2,t} = Sw * PwT{t} + J * x0;
		plot( ...
			Pxt2{2,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
			'color','cyan')

    end

    plot(X_Target,'FaceColor','none')

    title('Zero Input Reachable Sets for Each Mode')
    legend('$x_0$','Mode 1','Mode 2', '', '' , 'Target', ...
    	'Interpreter','latex')

    xlabel('$(x_k)_1$','Interpreter','latex')
	ylabel('$(x_k)_2$','Interpreter','latex')

    % Save to Gifs
    for mode_index = 1:lcsas0.n_modes
    	saveToSimpleGIF( ...
    		TimeHorizon , {Pxt2{mode_index,:}}, x0 , X_Target , ...
    		['results/hscc2022/images/pureRotationSystem_Mode' num2str(mode_index) '_ZeroInputSet.gif'] , ...
    		temp_axis , ...
    		mode_colors{mode_index} , ...
    		['Mode ' num2str(mode_index) ' Zero-Input Reachable Set'] )
    end

    results.ZeroInputTrajectoryTubes = Pxt2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Nice Solution Inputs Reachable Sets %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u_nice = {};
    u_nice{1,1} = [ 0 ; 0 ];
    u_nice{1,2} = [[ 0 ; 0 ];[0;-eta_u]];
    u_nice{1,3} = [[ 0 ; 0 ];[0;-eta_u];[eta_u;-eta_u]];
    u_nice{1,4} = [[ 0 ; 0 ];[0;-eta_u];[eta_u;-eta_u];[eta_u;0]];

    u_nice{2,1} = [ 0 ; 0 ]; 
    u_nice{2,2} = [[ 0 ; 0 ];[0;eta_u]];
    u_nice{2,3} = [[ 0 ; 0 ];[0;eta_u];[eta_u;eta_u]];
    u_nice{2,4} = [[ 0 ; 0 ];[0;eta_u];[eta_u;eta_u];[eta_u;0]];

    figure;
	hold on;

	scatter(x0(1),x0(2))

	for t = 1:TimeHorizon

		%Plot The Polyhedron.
        word1_prefix = lcsas0.L.words{1}(1:t);
		[Sw,Su,~,J,~] = lcsas0.get_mpc_matrices('word',word1_prefix);
		Pxt3{1,t} = Sw * PwT{t} + Su * u_nice{1,t} + J * x0;
		plot( ...
			Pxt3{1,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
			'color','salmon' ...
			)

		% Plot The Polyhedron For Word 2
		word2_prefix = lcsas0.L.words{2}(1:t);
		[Sw,Su,~,J,~] = lcsas0.get_mpc_matrices('word',word2_prefix);
		Pxt3{2,t} = Sw * PwT{t} + Su * u_nice{2,t} + J * x0;
		plot( ...
			Pxt3{2,t}.projection([size(Sw,1)-dim_x+1:size(Sw,1)]), ...
			'color','cyan')

    end

    plot(X_Target,'FaceColor','none')

    title('Nice Input Reachable Sets for Each Mode')
    legend('$x_0$','Mode 1','Mode 2', '', '' , 'Target', ...
    	'Interpreter','latex')
    
    xlabel('$(x_k)_1$','Interpreter','latex')
	ylabel('$(x_k)_2$','Interpreter','latex')

    % Save to Gifs
    for mode_index = 1:lcsas0.n_modes
    	saveToSimpleGIF( ...
    		TimeHorizon , {Pxt3{mode_index,:}}, x0 , X_Target , ...
    		['results/hscc2022/images/pureRotationSystem_Mode' num2str(mode_index) '_NiceInputSet.gif'] , ...
    		temp_axis , ...
    		mode_colors{mode_index} , ...
    		['Mode ' num2str(mode_index) ' Nice-Input Reachable Set'] )
    end

    results.NiceInputTrajectoryTubes = Pxt3;

    %% Results
    
    results.Parameters.LCSAS = lcsas0;
    
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
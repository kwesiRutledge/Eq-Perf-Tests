%draw_consistency_set_example1.m
%Description:
%	Based on observer_comparison94.m

%%%%%%%%%%%%%%%%%%%%%%
%% Input Processing %%
%%%%%%%%%%%%%%%%%%%%%%


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
	'ub', eta_u*ones(1,dim_x));
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

x0 = [0;0];

ad1 = Aff_Dyn( A1 , B1 , K1 , eye(dim_x) , Pw , Pv );
ad2 = Aff_Dyn( A2 , B2 , K2 , eye(dim_x) , Pw , Pv );

lcsas0 = LCSAS( [ad1,ad2], Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) , ...
				'X0' , Polyhedron('lb',x0','ub',x0') , ...
				'U'  , Pu  );

results.Parameters.LCSAS = lcsas0;

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

axis([-0.5 11 -0.5 6.5 ])

saveas(gcf,'images/similarRotationSystemComparison1','epsc')

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

axis([-0.5 11 -0.5 6.5 ])

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
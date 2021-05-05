function [lcsas_out,TimeHorizon,Pu,Pw,x0,Px0,X_Target] = get_opposing_rotation_lcsas(varargin)
	%Description:
	%	This system was created in observer_comparison93 and is an example of Subproblem 3.
	%	It should be possible to find a simple constant control law which satisfies the reachability problem.
	%	This is used to test a simple idea for CDC 2021.
	%
	%Usage:
	%	lcsas0 = get_opposing_rotation_lcsas()
	%	[lcsas_out,TimeHorizon,Pu,Pw,x0,Px0,X_Target] = get_opposing_rotation_lcsas()

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[TimeHorizon,x0,eta_w,eta_u,X_Target] = handle_similar_rotation_inputs(varargin{:});

	% By default,
	%	TimeHorizon = 4
	%	x0 = [-1;0];
	%	eta_w = 0.25;
	%	eta_v = 2*eta_w = 0.5
	%	X_Target = Polyhedron('lb',-(eta_w*TimeHorizon+eta*sqrt(2))*ones(1,dim_x), 'ub', (eta_w*TimeHorizon+eta_w*sqrt(2))*ones(1,dim_x) ) + [2.75;0];

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	twoDRotation = @(theta) [ cos(theta), -sin(theta) ; sin(theta), cos(theta) ];
	dim_x = 2;

	A1 = twoDRotation(pi/TimeHorizon);
	A2 = twoDRotation(-pi/TimeHorizon);

	B1 = eye(dim_x);
	B2 = eye(dim_x);

	% eta_w = 0.25;
	Pw = Polyhedron('lb',-eta_w*ones(1,dim_x),'ub',eta_w*ones(1,dim_x));
	Pv = Pw; %We won't use it.

	% eta_u = 2*eta_w;
	Pu = Polyhedron(...
		'lb',-eta_u*ones(1,dim_x), ...
		'ub', eta_u*ones(1,dim_x));

	Px0 = Polyhedron('lb',x0','ub',x0');

	% Create XT

	%X_Target = Polyhedron('lb',-(eta_w*TimeHorizon+eta*sqrt(2))*ones(1,dim_x), 'ub', (eta_w*TimeHorizon+eta_w*sqrt(2))*ones(1,dim_x) ) + [2.75;0];

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	ad1 = Aff_Dyn( A1 , B1 , zeros(dim_x,1) , eye(dim_x) , Pw , Pv );
	ad2 = Aff_Dyn( A2 , B2 , zeros(dim_x,1) , eye(dim_x) , Pw , Pv );

	lcsas_out = LCSAS( [ad1,ad2], Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) , 'X0' , Px0 , 'U' , Pu );

end

function [ TimeHorizon_out , x0 , eta_w , eta_u , X_Target ] = handle_similar_rotation_inputs( varargin )

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim_x = 2;

	%%%%%%%%%%%%%%%%%%%%
	%% Process Inputs %%
	%%%%%%%%%%%%%%%%%%%%

	if nargin > 0
		argidx = 1;
		while argidx <= nargin
			switch varargin{argidx}
				case 'TimeHorizon'
					TimeHorizon_out = varargin{argidx+1};
					argidx = argidx + 2;
				case 'x0'
					x0 = varargin{argidx+1};
					argidx = argidx + 2;
				case 'eta_w'
					eta_w = varargin{argidx+1};
					argidx = argidx+2;
				case 'eta_v'
					eta_v = varargin{argidx+1};
					argidx = argidx + 2;
				case 'X_Target'
					X_Target = varargin{argidx+1};
					argidx = argidx + 2;
				otherwise
					error(['Unexpected input given to the function get_cdc2021_similar_rotation_example(): ' varargin{argidx}])
			end
		end
	end

	%%%%%%%%%%%%%%%%%%
	%% Set Defaults %%
	%%%%%%%%%%%%%%%%%%

	if ~exist('TimeHorizon_out')
		TimeHorizon_out = 4;
	end

	if ~exist('x0')
		x0 = [-1;0];
	end

	if ~exist('eta_w')
		eta_w = 0.25;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Defining Dependent Variables %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	eta_u = 2*eta_w;

	target_width = (eta_w*TimeHorizon_out+eta_w*sqrt(2));

	% Compute Center after TimeHorizon_out rotations.
	target_center = [1/sqrt(2);1/sqrt(2)];
	target_center = [ 0 ; norm(target_center)] + [0;eta_u];
	target_center = norm(target_center)*[1/sqrt(2);1/sqrt(2)] + eta_u*ones(dim_x,1);
	target_center = [norm(target_center);0] + [eta_u;0];

	X_Target = Polyhedron('lb',-(target_width)*ones(1,dim_x), 'ub', (target_width)*ones(1,dim_x) ) + target_center;


end
function [lcsas_out,TimeHorizon,Pu,Pw,x0,Px0,X_Target] = get_cdc2021_similar_rotation_example(varargin)
	%Description:
	%	This system was created in observer_comparison94 and is an example of Subproblem 3.
	%	It should be possible to generate nonempty consistency sets for the behavior of this system
	%	under a variety of controllers.
	%
	%Usage:
	%	

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[TimeHorizon,x0] = handle_similar_rotation_inputs(varargin{:});

	% By default,
	%	TimeHorizon = 4
	%	x0 = [0;0];

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

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

	% TimeHorizon = 4;

	% Create PuT
	eta_u = 0.5*eta_w;
	Pu = Polyhedron(...
		'lb',-eta_u*ones(1,dim_x), ...
		'ub', eta_u*ones(1,dim_x))

	Px0 = Polyhedron('lb',x0','ub',x0');

	% Create XT
	X_Target = Polyhedron('lb',-eta_w*TimeHorizon*ones(1,dim_x), 'ub', eta_w*TimeHorizon*ones(1,dim_x) ) + [9.5;5];

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	ad1 = Aff_Dyn( A1 , B1 , K1 , eye(dim_x) , Pw , Pv );
	ad2 = Aff_Dyn( A2 , B2 , K2 , eye(dim_x) , Pw , Pv );

	lcsas_out = LCSAS( [ad1,ad2], Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) , 'X0' , Px0 );

end

function [ TimeHorizon_out , x0 ] = handle_similar_rotation_inputs( varargin )

	%%%%%%%%%%%%%%%%%%
	%% Set Defaults %%
	%%%%%%%%%%%%%%%%%%

	TimeHorizon_out = 4;
	x0 = [0;0];

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
				otherwise
					error(['Unexpected input given to the function get_cdc2021_similar_rotation_example(): ' varargin{argidx}])
			end
		end
	end

end
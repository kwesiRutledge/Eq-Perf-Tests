function [ lcsas , x0 , TimeHorizon , X_target ] = get_differently_loaded_drone_lcsas(varargin)
	%get_differently_loaded_drone_lcsas.m
	%Description:
	%	Creates the LCSAS representation of a one-dimensional drone with an unknown mass.
	%
	%Usage:
	%	[ lcsas , x0 , X_target ] = get_differetly_loaded_drone_lcsas()
	%	[ lcsas , x0 , X_target ] = get_differetly_loaded_drone_lcsas('dt',1)
	%	[ lcsas , x0 , X_target ] = get_differetly_loaded_drone_lcsas('TimeHorizon',10)

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ dt , m1 , m2 , TimeHorizon , x0 , X_target , U ] = ip_get_differently_loaded_drone_lcsas(varargin{:});

	g = 10; %m/s^2 % 9.8

	n_x = 2;
	n_w = 1;

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	A1 = [ 1 , dt ; 0 , 1 ];
	B1 = [ 0 ; dt/(m1) ];
	f1 = [ 0 ; -dt * g ];

	eta_w = 0.5;
	W1 = Polyhedron('lb',-eta_w,'ub',eta_w);

	ad1 = Aff_Dyn(	A1 , B1 , f1 , zeros(1,n_x), ...
					W1 , W1 , ...
					B1 , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	A2 = [ 1 , dt ; 0 , 1 ];
	B2 = [ 0 ; dt/(m2) ];
	f2 = [ 0 ; -dt * g ];

	% eta_w = 0.5;
	W2 = Polyhedron('lb',-eta_w,'ub',eta_w);

	ad2 = Aff_Dyn(	A2 , B2 , f2 , zeros(1,n_x), ...
					W2 , W2 , ...
					B2 , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	lcsas = LCSAS( ...
		[ad1,ad2] , ...
		Language([ones(1,TimeHorizon)],[2*ones(1,TimeHorizon)]), ...
		'X0', Polyhedron('lb',x0,'ub',x0) , ...
		'U' , U );

end

function [ dt , m1 , m2 , TimeHorizon , x0 , X_target , U ] = ip_get_differently_loaded_drone_lcsas(varargin)
	%Description:
	%	Process the inputs

	%%%%%%%%%%%%%%
	%% Defaults %%
	%%%%%%%%%%%%%%

	drone_settings = struct( ...
		'dt' , 0.1 , ...
		'm1' , 10 , ... %kilograms
		'm2' , 15 , ... %kilograms
		'TimeHorizon' , 8 , ...
		'x0' , [ 1 ; 0 ] , ...
		'X_target' , Polyhedron('lb',2,'ub',3) * Polyhedron('lb',-1,'ub',1) , ...
		'U' , Polyhedron('lb',-2,'ub',2) );

	%%%%%%%%%%%%%%%%
	%% Processing %%
	%%%%%%%%%%%%%%%%

	argin_index = 1;
	while argin_index < nargin
		switch varargin{argin_index}
			case 'dt'
				drone_settings.dt = varargin{argin_index+1};
				argin_index = argin_index + 2;
			case 'TimeHorizon'
				drone_settings.TimeHorizon = varargin{argin_index+1};
				argin_index = argin_index + 2;
			case 'x0'
				drone_settings.x0 = varargin{argin_index+1};
				argin_index = argin_index + 2;
			case 'X_target'
				drone_settings.X_target = varargin{argin_index+1};
				if ~isa(drone_settings.X_target,'Polyhedron')
					error('Input X_target must be a Polyhedron object.')
				end
				argin_index = argin_index + 2;
			otherwise
				error(['Unexpected input to get_differently_loaded_drone_lcsas(): ' varargin{argin_index} ])
		end

	end

	%%%%%%%%%%%%%%%%%%%%
	%% Create Outputs %%
	%%%%%%%%%%%%%%%%%%%%

	dt = drone_settings.dt;
	m1 = drone_settings.m1;
	m2 = drone_settings.m2;
	TimeHorizon = drone_settings.TimeHorizon;
	x0 = drone_settings.x0;
	X_target = drone_settings.X_target;
	U = drone_settings.U;

end
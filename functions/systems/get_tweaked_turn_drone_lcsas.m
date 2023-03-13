function [ lcsas , x0 , TimeHorizon , X_target ] = get_tweaked_turn_drone_lcsas(varargin)
	%get_tweaked_turn_drone_lcsas.m
	%Description:
	%	Creates the LCSAS representation of a drone in two-dimensional world where the turning may or may not be
	%	tweaked.
	%
	%Usage:
	%	[ lcsas , x0 , X_target ] = get_tweaked_turn_drone_lcsas()
	%	[ lcsas , x0 , X_target ] = get_tweaked_turn_drone_lcsas('dt',1)
	%	[ lcsas , x0 , X_target ] = get_tweaked_turn_drone_lcsas('TimeHorizon',10)

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ dt , theta1 , theta2 , TimeHorizon , x0 , X_target , U , eta_w ] = ip_get_tweaked_turn_drone_lcsas(varargin{:});

	n_x = 2;
	n_u = 2;
	n_w = 2;

	RotationMatrix = @(theta_in) [ cos(theta_in), - sin(theta_in) ; sin(theta_in) , cos(theta_in) ];

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	A = eye(n_x);

	A1 = A;
	B1 = RotationMatrix(theta1) * dt * eye(n_x);
	f1 = zeros(n_x,1);

	W1 = Polyhedron('lb',-eta_w*ones(1,n_w),'ub',eta_w*ones(1,n_w));

	ad1 = Aff_Dyn(	A1 , B1 , f1 , zeros(1,n_x), ...
					W1 , W1 , ...
					B1 , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	A2 = A;
	B2 = RotationMatrix(theta2) * dt * eye(n_x);
	f2 = zeros(n_x,1);

	% eta_w = 0.5;
	W2 = Polyhedron('lb',-eta_w*ones(1,n_w),'ub',eta_w*ones(1,n_w));

	ad2 = Aff_Dyn(	A2 , B2 , f2 , zeros(1,n_x), ...
					W2 , W2 , ...
					B2 , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	lcsas = LCSAS( ...
		[ad1,ad2] , ...
		Language(ones(1,TimeHorizon),2*ones(1,TimeHorizon)), ...
		'X0', Polyhedron('lb',x0','ub',x0') , ...
		'U' , U );

end

function [ dt , theta1 , theta2 , TimeHorizon , x0 , X_target , U , eta_w ] = ip_get_tweaked_turn_drone_lcsas(varargin)
	%Description:
	%	Process the inputs

	%%%%%%%%%%%%%%
	%% Defaults %%
	%%%%%%%%%%%%%%

	drone_settings = struct( ...
		'dt' , 0.1 , ...
		'theta1' , 0 , ... %radians
		'theta2' , pi/16 , ... %radians
		'TimeHorizon' , 10 , ...
		'x0' , [ 0 ; 0 ] , ...
		'X_target' , Polyhedron('lb',1.0,'ub',3.0) * Polyhedron('lb',1.0,'ub',3.0) , ...
		'eta_u' , 3  , ... %m/s
		'eta_w' , 0.1 ); %

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
            case 'theta1'
                drone_settings.theta1 = varargin{argin_index+1};
                argin_index = argin_index + 2;
			case 'theta2'
				drone_settings.theta2 = varargin{argin_index+1};
				argin_index = argin_index + 2;
			case 'eta_u'
				drone_settings.eta_u = varargin{argin_index+1};
				argin_index = argin_index + 2;
			case 'eta_w'
				drone_settings.eta_w = varargin{argin_index+1};
				argin_index = argin_index + 2;
			otherwise
				error(['Unexpected input to get_differently_loaded_drone_lcsas(): ' varargin{argin_index} ])
		end

	end

	%%%%%%%%%%%%%%%%%%%%
	%% Create Outputs %%
	%%%%%%%%%%%%%%%%%%%%

	dt = drone_settings.dt;
	theta1 = drone_settings.theta1;
	theta2 = drone_settings.theta2;
	TimeHorizon = drone_settings.TimeHorizon;
	x0 = drone_settings.x0;
	X_target = drone_settings.X_target;
	U = Polyhedron('lb',-drone_settings.eta_u*ones(1,2),'ub',drone_settings.eta_u*ones(1,2));
	eta_w = drone_settings.eta_w;

end
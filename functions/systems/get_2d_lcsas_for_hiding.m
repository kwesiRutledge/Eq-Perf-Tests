function [ lcsas , x0 , TimeHorizon , X_target ] = get_2d_lcsas_for_hiding(varargin)
	%get_scalar_lcsas_for_hiding.m
	%Description:
	%	Creates the LCSAS representation of a scalar system and its reachability task.
	%	It should be possible to satisfy the task without revealing which mode the system
	%	is in.
	%
	%Usage:
	%	[ lcsas , x0 , X_target ] = get_tweaked_turn_drone_lcsas()
	%	[ lcsas , x0 , X_target ] = get_tweaked_turn_drone_lcsas('TimeHorizon',10)

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	[ B1 , B2 , TimeHorizon , x0 , X_target , U , W ] = ip_get_tweaked_turn_drone_lcsas(varargin{:});

	n_x = 2;
	n_u = 2;
	n_w = 2;

	RotationMatrix = @(theta_in) [ cos(theta_in), - sin(theta_in) ; sin(theta_in) , cos(theta_in) ];

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	A = eye(n_x);

	A1 = A; f1 = zeros(n_x,1); W1 = W; %B2 = 1

	ad1 = Aff_Dyn(	A1 , B1 , f1 , zeros(1,n_x), ...
					W1 , W1 , ...
					eye(n_w) , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	A2 = A; f2 = zeros(n_x,1);	W2 = W; % B2 = 2;

	ad2 = Aff_Dyn(	A2 , B2 , f2 , zeros(1,n_x), ...
					W2 , W2 , ...
					eye(n_w) , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	lcsas = LCSAS( ...
		[ad1,ad2] , ...
		Language(ones(1,TimeHorizon),2*ones(1,TimeHorizon)), ...
		'X0', Polyhedron('lb',x0','ub',x0') , ...
		'U' , U );

end

function [ B1 , B2 , TimeHorizon , x0 , X_target , U , W ] = ip_get_tweaked_turn_drone_lcsas(varargin)
	%Description:
	%	Process the inputs

	%%%%%%%%%%%%%%
	%% Defaults %%
	%%%%%%%%%%%%%%

	system_settings = struct( ...
		'b1' , eye(2) , ...
		'b2' , [-1, 0; 0, 1] , ... 
		'TimeHorizon' , 2 , ...
		'x0' , [0;0] , ...
		'X_target' , Polyhedron('lb',[1.5,1.3],'ub',[2.5,2.3]) , ...
		'eta_u' , 0.2  , ... %m/s
		'W' , Polyhedron('lb',[0.8,0.8],'ub',[1.2,1.2] ) );

	n_x = 2;
	n_u = 2;
	n_w = 2;

	%%%%%%%%%%%%%%%%
	%% Processing %%
	%%%%%%%%%%%%%%%%

	argin_index = 1;
	while argin_index < nargin
		switch varargin{argin_index}
			case 'TimeHorizon'
				system_settings.TimeHorizon = varargin{argin_index+1};
				argin_index = argin_index + 2;
			case 'x0'
				system_settings.x0 = varargin{argin_index+1};
				argin_index = argin_index + 2;
			case 'X_target'
				system_settings.X_target = varargin{argin_index+1};
				if ~isa(system_settings.X_target,'Polyhedron')
					error('Input X_target must be a Polyhedron object.')
				end
				argin_index = argin_index + 2;
            case 'b1'
                system_settings.b1 = varargin{argin_index+1};
                argin_index = argin_index + 2;
			case 'b2'
				system_settings.b2 = varargin{argin_index+1};
				argin_index = argin_index + 2;
			case 'eta_u'
				system_settings.eta_u = varargin{argin_index+1};
				argin_index = argin_index + 2;
            case 'W'
                system_settings.W = varargin{argin_index+1};
                argin_index = argin_index + 2;
			otherwise
				error(['Unexpected input to get_differently_loaded_drone_lcsas(): ' varargin{argin_index} ])
		end

	end

	%%%%%%%%%%%%%%%%%%%%
	%% Create Outputs %%
	%%%%%%%%%%%%%%%%%%%%

	B1 = system_settings.b1;
	B2 = system_settings.b2;
	TimeHorizon = system_settings.TimeHorizon;
	x0 = system_settings.x0;
	X_target = system_settings.X_target;
	U = Polyhedron('lb',-system_settings.eta_u*ones(1,n_u),'ub',system_settings.eta_u*ones(1,n_u));
	W = system_settings.W;

end
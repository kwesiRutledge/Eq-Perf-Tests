function [ lcsas_out , P_u , Pw1 , eta_v , x0 , P_target ] = get_uncertain_thick_pendulum_lcsas( varargin )
	%Description:
	%	Create a one-dimensional system to quickly synthesize belief trees.
	%
	%Usage:
	%	[lcsas_out,P_u] = get_uncertain_thick_pendulum_lcsas()
	%	[lcsas_out,P_u,Pw1,Pw2,eta_v,eta_x0,P_target] = get_uncertain_thick_pendulum_lcsas()

	%% Input Processing %%

	[ lcsas_settings , TimeHorizon , dt ] = get_uncertain_thick_pendulum_lcsas_ip( varargin{:} );	

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	itp0 = InvertedThickPendulum();

	n_x = 2; 
	n_y = 1; n_v = 1;


	n_w1 = 1;
	Pw1 = Polyhedron('lb',-lcsas_settings.eta_w*ones(1,n_w1),'ub',lcsas_settings.eta_w*ones(1,n_w1));
	
	% n_w2 = 2;
	% Pw1 = Polyhedron('lb',-lcsas_settings.eta_w*[0.1,1],'ub',lcsas_settings.eta_w*[0.1,1]);

	eta_v = 0.0;
	Pv = Polyhedron('lb',-eta_v*ones(1,n_v),'ub',eta_v*ones(1,n_v));

	x0 = itp0.x();
	P_x0 = Polyhedron('lb',x0','ub',x0');
	P_u = Polyhedron('lb',-lcsas_settings.eta_u,'ub',lcsas_settings.eta_u);


	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	% For normal itp, create linearized A,B matrices
	[A,B,K] = itp0.LinearizedDiscreteDynamicsAbout(itp0.x,0,dt);
	Bw = B;
	ad1 = Aff_Dyn(A,B,K,ones(n_y,n_x),Pw1,Pv,Bw,zeros(n_y,n_v));

	% For changed CoM, consider the controller.
	itp0.CoMx_rel = -0.25;
	[A,B,K] = itp0.LinearizedDiscreteDynamicsAbout(itp0.x,0,dt);

	ad2 = Aff_Dyn(A,B,K,ones(n_y,n_x),Pw1,Pv,Bw,zeros(n_y,n_v));

	lcsas_out = LCSAS( [ad1,ad2] , Language(1*ones(1,TimeHorizon),2*ones(1,TimeHorizon)) , 'X0' , P_x0, 'U' , P_u );

	P_target = lcsas_settings.P_target;

end

function [ lcsas_settings , TimeHorizon , dt ] = get_uncertain_thick_pendulum_lcsas_ip( varargin )
	%Description:
	%	Parses the potential nonstandard inputs to the lcsas.

	%% Default

	lcsas_settings = struct( ...
		'eta_w', 0.1 , ...
		'eta_u', 5.0 , ...
		'TimeHorizon', 3 , ...
		'P_target', Polyhedron('lb',[-0.5,-10],'ub',[0.5,10]) , ...
		'dt', 0.1 ...
		);

	%% Input Processing %%

	argument_index = 1;
	while nargin >= argument_index
		switch varargin{argument_index}
		case 'TimeHorizon'
			lcsas_settings.TimeHorizon = varargin{argument_index+1};
			argument_index = argument_index + 2;
		case 'dt'
			lcsas_settings.dt = varargin{argument_index+1};
			argument_index = argument_index + 2;
		case 'P_target'
			lcsas_settings.P_target = varargin{argument_index+1};
			argument_index = argument_index + 2;
		case 'eta_u'
			lcsas_settings.eta_u = varargin{argument_index+1};
			argument_index = argument_index + 2;
		otherwise
		 	error(['Unexpeced input to get_uncertain_thick_pendulum_lcsas: ' varargin{argument_index} ])
		end
	end


	%% Creating Some Outputs %%

	TimeHorizon = lcsas_settings.TimeHorizon;
	dt 			= lcsas_settings.dt;

end
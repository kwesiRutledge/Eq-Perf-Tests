function [ lcsas , x0 , TimeHorizon , X_target ] = get_1d_lcsas_hscc1(varargin)
	%get_1d_lcsas_hscc1.m
	%Description:
	%	Creates a simple, one-dimensional LCSAS to illustrate the concept of "differentiating no matter what".
	%
	%Usage:
	%	[ lcsas , x0 , TimeHorizon , X_target ] = get_1d_lcsas_hscc1()
	%	[ lcsas , x0 , TimeHorizon , X_target ] = get_1d_lcsas_hscc1()

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	n_x = 1;
	n_w = 1;

	eta_w = 0.1;
	f = 2*eta_w;
	eta_u = 1.25*(f+eta_w);

	TimeHorizon = 5;

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	A1 = 1;
	B1 = 1;
	f1 = f;

	W1 = Polyhedron('lb',-eta_w,'ub',eta_w);

	ad1 = Aff_Dyn(	A1 , B1 , f1 , zeros(1,n_x), ...
					W1 , W1 , ...
					B1 , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	A2 = 1;
	B2 = 1;
	f2 = -f;

	% eta_w = 0.5;
	W2 = Polyhedron('lb',-eta_w,'ub',eta_w);

	ad2 = Aff_Dyn(	A2 , B2 , f2 , zeros(1,n_x), ...
					W2 , W2 , ...
					B2 , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%
	x0 = 0;
	U = Polyhedron('lb',-eta_u,'ub',eta_u);

	lcsas = LCSAS( ...
		[ad1,ad2] , ...
		Language(ones(1,TimeHorizon),2*ones(1,TimeHorizon)), ...
		'X0', Polyhedron('lb',x0','ub',x0') , ...
		'U' , U );

	%%%%%%%%%%%%%%%%%%%
	%% Create Target %%
	%%%%%%%%%%%%%%%%%%%

	eta_target = 2*TimeHorizon*eta_w;
	X_target = Polyhedron('lb',-eta_target,'ub',eta_target);
	

end
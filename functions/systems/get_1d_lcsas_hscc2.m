function [ lcsas , x0 , TimeHorizon , X_target ] = get_1d_lcsas_hscc2(varargin)
	%get_1d_lcsas_hscc2.m
	%Description:
	%	Creates a simple, one-dimensional LCSAS to illustrate the concept of "differentiating only when smart".
	%
	%Usage:
	%	[ lcsas , x0 , TimeHorizon , X_target ] = get_1d_lcsas_hscc2()
	%	[ lcsas , x0 , TimeHorizon , X_target ] = get_1d_lcsas_hscc2()

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	n_x = 1;
	n_w = 1;

	eta_w = 0.5;
	f1 = 0;
	f2 = 0;
	eta_u = 10*(eta_w);

	TimeHorizon = 4;

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	W1 = Polyhedron('lb',-eta_w,'ub',eta_w);

	ad1 = Aff_Dyn(	1 , 1 , f1 , zeros(1,n_x), ...
					W1 , W1 , ...
					1 , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%%
	%% Create Mode 1 %%
	%%%%%%%%%%%%%%%%%%%

	W2 = W1;

	ad2 = Aff_Dyn(	-1 , 1 , f2 , zeros(1,n_x), ...
					W2 , W2 , ...
					1 , zeros(1,1) );

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%
	x0 = 0.25*eta_w;
	U = Polyhedron('lb',-eta_u,'ub',eta_u);

	lcsas = LCSAS( ...
		[ad1,ad2] , ...
		Language(ones(1,TimeHorizon),2*ones(1,TimeHorizon)), ...
		'X0', Polyhedron('lb',x0','ub',x0') , ...
		'U' , U );

	%%%%%%%%%%%%%%%%%%%
	%% Create Target %%
	%%%%%%%%%%%%%%%%%%%

	eta_target = 5*eta_w;
	X_target = Polyhedron('lb',-eta_target,'ub',eta_target);
	

end
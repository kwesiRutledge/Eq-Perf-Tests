function [lcsas_out,P_u,Pw1,Pw2,eta_v,eta_x0] = get_1d_lcsas( varargin )
	%Description:
	%	Create a one-dimensional system to quickly synthesize belief trees.
	%
	%Usage:
	%	[lcsas_out,P_u] = get_1d_lcsas()

	%% Constants %%

	a = 0.9;
	b = 1;

	Pw1 = Polyhedron('lb',0.0,'ub',0.5);
	Pw2 = Polyhedron('lb',-0.5,'ub',0.0);

	eta_v = 0.25;
	Pv = Polyhedron('lb',-eta_v,'ub',eta_v);

	eta_x0 = 0.25;
	P_x0 = Polyhedron('lb',-eta_x0,'ub',eta_x0);

	%% Algorithm %%

	eta_u = 5;
	P_u = Polyhedron('lb',-eta_u,'ub',eta_u);

	ad1 = Aff_Dyn(a,b,zeros(1),eye(1),Pw1,Pv);
	ad2 = Aff_Dyn(a,b,zeros(1),eye(1),Pw2,Pv);

	lcsas_out = LCSAS( [ad1,ad2] , Language([1,1,1],[2,2,2],[1,2,1]) , 'X0' , P_x0 );

end
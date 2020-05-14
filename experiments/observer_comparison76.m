function [results] = observer_comparison76( varargin )
	%observer_comparison76.m
	%Description:
	%	Comparing the method for detecting projection inclusion that I have against sadra's condition.
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	% if nargin >= 1
	% 	load_data_flag = varargin{1};
	% end

	% if nargin >= 3
	% 	c_sq.dim_x = varargin{2};
	% 	c_sq.dim_y = varargin{3};
	% end

	% if nargin >= 4
	% 	verbosity = varargin{4};
	% end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	test_name = 'observer_comparison76';
	save_file_name = [ 'results/' test_name '_results.mat'];

	%% Create System

	L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
	T = 4;

	%Create a simple Language Constrainted Switching System
	A1 = [0,1;0.1,-0.05];
	B1 = [0;1];
	C1 = [1,0];
	
	n_x = size(A1,1);
	n_u = size(B1,2);
	n_y = size(C1,1);

	eta_v = 0.1; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
	eta_v = 0.3;
	Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

	eta_u = 0; eta_x0 = 0.3;
	P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

	f1 = eta_w*[0;1];
	f2 = eta_w*[1;0];
	f3 = -f1;
	f4 = -f2;

	aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
						Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

	lcsas0 = LCSAS( aff_dyn_list , L1 );

	results.lcsas = lcsas0;

	%% Debugging Variables

	verbosity = 0;
	update_freq = 50;

	results.parameters.dim = size(A1,1);

	%%%%%%%%%%%
	%% Tests %%
	%%%%%%%%%%%

	disp(' ')
	disp(['Beginning ' test_name '.' ])
	disp('The objective of this test is to use Sadraddini''s Polytope Inclusion Constraint')
	disp('to construct a BeliefGraph.' )
	disp(' ')

	%% Get System %%

	%% Create Initial Nodes %%



	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Save The Larger Variables to a Data File %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%save([ save_file_name '.mat'])

end
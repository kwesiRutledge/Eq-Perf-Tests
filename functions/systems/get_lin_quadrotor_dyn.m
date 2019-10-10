function [ quad_lcsas ] = get_lin_quadrotor_dyn(varargin)
%Description:
	%	Creates the dynamics for a quadrotor helicopter modeled after the Ascending Technologies' ResearchPilot.
	%	System should be a 10 dimensional.
	%
	%
	%Example Usage:
	%	ad = get_lin_quadrotor_sys()
	
	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	allowable_nargins = [0,1];
	if ~any(allowable_nargins==nargin)
		error(['This function does not support having ' num2str(nargin) ' inputs!'])
	end



	switch nargin
		case 0
			dt = 1;
		case 1
			dt = varargin{1};
		otherwise
			error('It looks like this number of arguments isn''t currently supported. :(')
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Constants (Copied from Dustin Webb's Github) %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	prop_diameter = 0.1397;
	rotor_dist = 0.25;
	w = rotor_dist+prop_diameter; %width (m)
	h = 0.1; %height (m)
	d = rotor_dist+prop_diameter; %depth (m)
	m = 0.26; %mass (kg)

	rotor_to_com = rotor_dist/2;

	g = 9.81;

	%assume a spherical cow in a vacuum
	%or in this case a solid cuboid for the moments of inertia
	%I_h = 1/12*m*(w^2+d^2);
	I_w = 1/12*m*(h^2+d^2);
	I_d = 1/12*m*(h^2+w^2);

	state_dims = 10;
	input_dims = 3;

	n_w = 2;
	eta_w = 0.5;
	Pw1 = Polyhedron('ub',eta_w*ones(1,n_w),'lb',-eta_w*ones(1,n_w));
	Pw2 = Pw1 + (eta_w/2)*ones(n_w,1);

	n_v = state_dims;
	eta_v = 0.3;
	Pv = Polyhedron('ub',eta_v*ones(1,n_v),'lb', -eta_v*ones(1,n_v));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Creating the Dynamics Array %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	A=zeros(state_dims);
	A(1:3,4:6) = eye(3);
	A(4:6,7:8) = [0,g;-g,0;0,0];
	A(7:8,9:10) = eye(2);

	B=zeros(state_dims, input_dims);

	B(4:6,1) = [0,0,1/m]';
	B(9:10,2:3) = diag([rotor_to_com/I_w, rotor_to_com/I_d]);

	c = zeros(state_dims,1);

	B_w = zeros(state_dims,2);
	B_w(3+[1:2],[1:2]) = eye(2);

	% Discretize elements
	sys1 = ss(A,B,eye(state_dims),0);
	dsys1 = c2d(sys1,dt);

	sys2 = ss(A,B_w,eye(state_dims),0);
	dsys2 = c2d(sys2,dt);

	A_d = dsys1.A;
	B_d = dsys1.B;
	c_d = c; %Because c=0
	Bw_d = dsys2.B;

	ad_nominal = Aff_Dyn(	A_d,B_d,c_d,eye(state_dims), ...
							Pw1, Pv, ...
							Bw_d , eye(n_v) );
	ad_strongW = Aff_Dyn(	A_d,B_d,c_d,eye(state_dims), ...
							Pw2, Pv, ...
							Bw_d , eye(n_v) );

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Creating the Mode Language L %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	lang_length = 6;
	strongw_start = 4;
	word_arr = {};
	for word_idx = 1:(lang_length-strongw_start+1)
		word_arr{word_idx} = [ones(1,strongw_start-1+(word_idx-1)),2*ones(1,lang_length-strongw_start+1-(word_idx-1))];
	end
	L = Language(word_arr);

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	quad_lcsas = LCSAS([ad_nominal,ad_strongW],L);


end
function [lk_dyn_cont, lk_dyn_disc] = get_lk_aff_dyn(varargin)
	%get_lk_dyn.m
	%Summary:
	%	Returns the continuous and discretized dynamics of the lane keeping system systems when the
	%	discretization constant is given.
	%
	%Usage:
	%	[cont_sys,disc_sys] = get_lk_aff_dyn()
	%	[cont_sys,disc_sys] = get_lk_aff_dyn(delta_t)
	%
	%Inputs:
	%	None
	%
	%Outputs:
	%	lk_dyn_cont - The Continuous Time system
	%				  Aff_Dyn object
	%
	%	lk_dyn_disc - The Discretized System
	%				  Aff_Dyn object

	%% Constants
	if nargin == 0
		delta_t = 0.1;
	else
		delta_t = varargin{1};
	end

	m = 1650; %kg
	C_af = 133*10^3; %N
	C_ar = 98.8*10^3; %N
	u = 25; %m/s
	a = 1.11; %m
	b = 1.59; %m
	I_z = 2315; %kg m^2

	%% Define Continuous Time Matrices
	n_x = 4;
	A = [	0,	1,							u,	0;
			0,	-(C_af+C_ar)/(m*u),			0,	(b*C_ar-a*C_af)/(m*u) - u;
			0,	0,							0,	1;
			0,	(b*C_ar - a*C_af)/(I_z*u),	0,	-(a^2*C_af+b^2*C_ar)/(I_z*u)];

	B = [	0;
			C_af/m;
			0;
			a*(C_af/I_z)];

	B_w = [0;0;-1;0];

	% Assume that the system can only observe:
	%	1. the lateral deviation from the  of the lane and
	%	2. roughly the difference in heading between the lane and the vehicle's direction of travel
	C = eye(n_x);%[1,0,0,0;0,0,1,0;0,0,0,1];
	n_y = size(C,1);

	%Description of the Noise
	%Process Noise (road curvature) is bounded by 0.05
	eta_w = 0.05;
	%Measurement noise is unclear
	eta_v = 0.001;

	lk_dyn_cont = Aff_Dyn(A,B,zeros(n_x,1),C,eta_w,eta_v, B_w, eye(n_y));

	%% Discretization
	temp_sys1 = ss(A,B,C,0);
	temp_sys2 = ss(A,B_w,C,0);

	temp_dsys1 = c2d(temp_sys1,delta_t);
	temp_dsys2 = c2d(temp_sys2,delta_t);

	lk_dyn_disc = Aff_Dyn(	temp_dsys1.A,temp_dsys1.B,zeros(n_x,1),C,...
							eta_w,eta_v,...
							temp_dsys2.B, eye(n_y));


end
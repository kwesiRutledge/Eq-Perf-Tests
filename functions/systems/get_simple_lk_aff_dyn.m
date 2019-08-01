function [simple_LK] = get_simple_lk_aff_dyn(dt)
	%System Parameters
	A = [ zeros(2,1) [1;-20] ];
	n_x = size(A,1);

	B = [ 0; 1];
	F = zeros(n_x,1);
	
	C = eye(2); %[1 0];
	n_y = size(C,1);

	temp_sys = ss(A,B,C,0);
	dt = 0.1;
	temp_dsys = c2d(temp_sys,dt);

	temp_sys2 = ss(A,[1;0],C,0);
	temp_dsys2 = c2d(temp_sys2,dt);

	eta_w = 0.05; eta_v = 0.1;

	simple_LK = Aff_Dyn(	temp_dsys.A,temp_dsys.B,F,C,...
							eta_w,eta_v, ...
							temp_dsys2.B, eye(n_y) );
end
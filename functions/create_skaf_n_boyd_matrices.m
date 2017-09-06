function [G,H,C_big,x0_mat] = create_skaf_n_boyd_matrices(sys,T)

	%Constants

	%Find G,H, and x0_mat Matrices
	G 		= calc_w_effect_mat(sys.A,T);
	H 		= calc_u_effect_mat(sys.A,sys.B,T);
	x0_mat 	= calc_x0_mat(sys.A,sys.x0,T);

	%Calculate Big C Matrix
	C_at_each_n = {};
	for i = 1:T
		C_at_each_n{i} = sys.C; 
	end

	C_big = [ blkdiag(C_at_each_n{:}) zeros( size(sys.C) * [ T 0 ; 0 1 ] ) ];

end
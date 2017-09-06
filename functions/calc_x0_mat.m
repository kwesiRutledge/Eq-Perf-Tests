function [x0_mat] = calc_x0_mat(A,x0,T)

	%Constants

	%Creating matrix.
	x0_mat = x0;

	for i = 1 : T
		x0_mat = [x0_mat; (A^i)* x0];
	end
end
function [ w , err ] = find_w_that_completes( ad , x_t , u_t , x_tp1 )
	%Description:
	%	Given an Aff_Dyn object, find a value w that satisfies
	%	
	%		x_tp1 = A x + B u + w + f
	%	
	%	where the matrices are given by the affine dynamics object.
	%
	%
	%Usage:
	%	w = ad.find_w_that_completes( x_t , u_t , x_tp1 )
	%	[w,err] = ad.find_w_that_completes( x_t , u_t , x_tp1 )

	%% Constants
	A = ad.A;
	B = ad.B;
	B_w = ad.B_w;
	f = ad.f;

	%% Algorithm
	w = pinv(B_w) * ( x_tp1 - A * x_t - B * u_t - f );

	err = x_tp1 - A * x_t - B * u_t - B_w * w - f;

end
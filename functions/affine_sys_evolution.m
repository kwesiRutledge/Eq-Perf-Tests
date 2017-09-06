function [x_conc] = affine_sys_evolution(sys,ctrlr,T,delta,mu)
	% UNFINISHED
	% This functions should allow the system's evolution to occur

	%Constants
	A = sys.A;
	B = sys.B;
	C = sys.C;
	E = sys.E;

	x0 = sys.x0;
	m = sys.m;
	d = sys.d;

	K = ctrlr.K;
	L = ctrlr.L;

	%Start to apply
	%x = x0;
	x_conc = [x0;zeros(size(x0))];

	for t = 1 : T
		%Compute the next part of x_concatenated which is [x(t);x_hat(t)]
		x_conc(:,t+1) = [ 	A 	B*K ; L*C A+B*K-L*C ] * x_conc(:,t) + [E*delta(t);L*mu(:,t)];
	end

end
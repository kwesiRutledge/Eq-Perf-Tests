classdef Aff_Dyn
	% Aff_Dyn( A, B, f , C )
	% Aff_Dyn( A, B, f , C, eta_w , eta_v )
	% Aff_Dyn( A, B, f , C, eta_w , eta_v, B_w , C_v  )
	% Dscrete time system is of the form
	%
	%	x(k+1) = A x(k) + B u(k) + Bw w(k) + f
	%	y(k)   = C x(k) + C_v v(k)
	%

	properties
		A;
		B;
		B_w;
		f;
		C;
		C_v;
		eta_w;
		eta_v;
		P_w;
		P_v;
		x0;
	end

	methods
		%Constructor
		function ad = Aff_Dyn(A, B, f , C, w_info , v_info, B_w , C_v)

			%Constants
			allowed_nargins = [ 4 , 6 , 8 ];

			%Input processing
			if ~any(nargin == allowed_nargins)
				error([ num2str(nargin) ' number of inputs is not supported.']);
			end

			%Assign basic variables
			ad.A = A;
			ad.B = B;
			ad.f = f;
			ad.C = C;

			n = size(A,1);
			p = size(C,1);
			m = size(B,2);

			ad.x0 = zeros(n,1);

			switch nargin
			case 4
				ad.B_w = eye(n);
				ad.C_v = eye(p);
				ad.eta_w = 1;
				ad.eta_v = 1;
			case 6
				ad.B_w = eye(n);
				ad.C_v = eye(p);
				if isa(w_info,'Polyhedron')
					ad.eta_w = NaN;
					ad.eta_v = NaN;

					ad.P_w = w_info;
					ad.P_v = v_info;
				else
					ad.eta_w = eta_w;
					ad.eta_v = eta_v;
				end
			case 8
				ad.B_w = B_w;
				ad.C_v = C_v;
				if isa(w_info,'Polyhedron')
					ad.eta_w = NaN;
					ad.eta_v = NaN;

					ad.P_w = w_info;
					ad.P_v = v_info;
				else
					ad.eta_w = eta_w;
					ad.eta_v = eta_v;
				end
			end

		end
	end

end
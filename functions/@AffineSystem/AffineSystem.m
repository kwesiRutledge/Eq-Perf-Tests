classdef AffineSystem
	%Definition:
	%	Object which defines discrete time affine systems of the form
	%
	%		x(k+1) = A x(k) + B u(k) + Bw w(k) + f
	%		y(k)   = C x(k) + C_v v(k)
	%
	%Usage:
	%	ad = AffineSystem( A, B, f )
	%	ad = AffineSystem( A, B, f , C )
	%	ad = AffineSystem( A, B, f , eta_w ,  C , eta_v )
    %	ad = AffineSystem( A, B, f , Poly_w , C , Poly_v )
	% 	ad = AffineSystem( A, B, f , eta_w ,  C , eta_v, B_w , C_v  )
    %	ad = AffineSystem( A, B, f , Poly_w , C , Poly_v, B_w , C_v  )
    %
	%
	%FAQ:
	%	Q - If I don't provide the value for B_w or C_v are they assigned to identity, zero, or some other matrix?
	%	A - If B_w or C_v are not provided, then they are assumed to be identity matrices.

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
		function ad = AffineSystem( A, B, f , w_info , C , v_info, B_w , C_v)

			% Constants
			% =========
			allowed_nargins = [ 3 , 4 , 6 , 8 ];
            
			% Input processing
			% ================

			if ~any(nargin == allowed_nargins)
				error([ num2str(nargin) ' number of inputs is not supported.']);
            end
            
            if nargin >= 6
            	%If the disturbance w and v info has been given.
	            if ~((isscalar(w_info) && isscalar(v_info)) || ...
	                 (isa(w_info,'Polyhedron') && isa(v_info,'Polyhedron')) )
	                error('w_info and v_info may only be scalars or ''Polyhedron'' objects.')
	            end

	        end

			% Define Matrices
			% ===============
			
			ad.A = A;
			ad.B = B;
			ad.f = f;

			n = size(A,1);
			m = size(B,2);

			if nargin >= 6
				ad.C = C;
			else
				ad.C = zeros(n);
			end

			p = size(ad.C,1);

			% Define Disturbance Sets
			% =======================

			ad.x0 = zeros(n,1);

			%

			switch nargin
			case 3
				ad.C = zeros(n);
				p = size(ad.C,1);

				ad.B_w = eye(n);
				ad.C_v = eye(p);
				ad.eta_w = 1;
				ad.eta_v = 1;
                
                %Create the infinity norm balls that the above eta_w
                %and eta_v define.
                ad.P_w = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));
                ad.P_v = Polyhedron('lb',-ones(1,p),'ub',ones(1,p));

			case 4
				ad.C = C;
				p = size(C,1);

				ad.B_w = eye(n);
				ad.C_v = eye(p);
				ad.eta_w = 1;
				ad.eta_v = 1;
                
                %Create the infinity norm balls that the above eta_w
                %and eta_v define.
                ad.P_w = Polyhedron('lb',-ones(1,n),'ub',ones(1,n));
                ad.P_v = Polyhedron('lb',-ones(1,p),'ub',ones(1,p));
			case 6
				ad.C = C;
				p = size(C,1);

				ad.B_w = eye(n);
				ad.C_v = eye(p);
				if isa(w_info,'Polyhedron')
					ad.eta_w = NaN;
					ad.eta_v = NaN;

					ad.P_w = w_info;
					ad.P_v = v_info;
				else
					ad.eta_w = w_info;
					ad.eta_v = v_info;
                    
                    %Create the infinity norm balls that the above eta_w
                    %and eta_v define.
                    ad.P_w = w_info*Polyhedron('lb',-ones(1,n),'ub',ones(1,n));
                    ad.P_v = v_info*Polyhedron('lb',-ones(1,p),'ub',ones(1,p));
				end
			case 8
				ad.C = C;
				p = size(C,1);

				ad.B_w = B_w;
				ad.C_v = C_v;
				if isa(w_info,'Polyhedron')
					ad.eta_w = NaN;
					ad.eta_v = NaN;

					ad.P_w = w_info;
					ad.P_v = v_info;
				else
					ad.eta_w = w_info;
					ad.eta_v = v_info;
                    
                    %Create the infinity norm balls that the above eta_w
                    %and eta_v define.
                    ad.P_w = w_info*Polyhedron('lb',-ones(1,size(B_w,2)),'ub',ones(1,size(B_w,2)));
                    ad.P_v = v_info*Polyhedron('lb',-ones(1,size(C_v,2)),'ub',ones(1,size(C_v,2)));
				end
			end

        end
	end

end
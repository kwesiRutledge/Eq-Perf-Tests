function [G,H,C_big,x0_mat] = create_skaf_n_boyd_matrices(varargin)
	% create_skaf_n_boyd_matrices
	%	Description:
	%		Creates MPC matrices that are used in Skaf and Boyd's work.
	%		This is modified to allow for missing observations if necessary.
	%
	%	Usage:
	%		create_skaf_n_boyd_matrices(sys,T)
	%		create_skaf_n_boyd_matrices(sys,T,'missing',m_locs)
	%	
	%	Inputs:
	%		sys - 		A struct containing system matrices and an initial condition
	%					for the desired system.
	%					Required fields: .A,.B,.C,.x0
	%
	%		T - 		Time horizon for the MPC matrices.
	%
	%		text_code - A string that signals for special behavior.
	%					Acceptable values: 'missing'	
	%
	%		m_locs - 	An array containing the indices of the 'missing' data occurences.
	%
	%	Outputs:
	%		G - 

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if (nargin ~= 2) & (nargin ~= 4)
		error(['Inappropriate number of arguments. (Received ' num2str(nargin) ')'])
	end

	sys = varargin{1};
	T 	= varargin{2};

	if nargin > 2
		text_code = varargin{3};
		m_locs = varargin{4};
	end

	%Constants
	if ~exist('m_locs')
		m_locs = [];
	end

	%Find G,H, and x0_mat Matrices
	G 		= calc_w_effect_mat(sys.A,T);
	H 		= calc_u_effect_mat(sys.A,sys.B,T);
	x0_mat 	= calc_x0_mat(sys.A,sys.x0,T);

	%Calculate Big C Matrix
	C_at_each_n = {};
	for i = 1:T
		if any(m_locs == (i-1))
			C_at_each_n{i} = zeros(size(sys.C));
		else
			C_at_each_n{i} = sys.C;
		end 
	end

	C_big = [ blkdiag(C_at_each_n{:}) zeros( size(sys.C) * [ T 0 ; 0 1 ] ) ];

end
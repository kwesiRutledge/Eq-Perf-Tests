function [varargout] = create_skaf_n_boyd_matrices(varargin)
	% create_skaf_n_boyd_matrices
	%	Description:
	%		Creates MPC matrices that are used in Skaf and Boyd's work.
	%		This is modified to allow for missing observations if necessary.
	%
	%	Usage:
	%		[G,H,C_big,x0_mat] = create_skaf_n_boyd_matrices(sys,T)
	%		[G,H,C_big,x0_mat] = create_skaf_n_boyd_matrices(sys,T,'missing',m_locs)
	%		[G,H,C_big,x0_mat,E_bar,G_bar] = create_skaf_n_boyd_matrices(sys,T)
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

	if (nargin ~= 2) & (nargin ~= 4) & (nargin ~=6)
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

	%%%%%%%%%%%%%%%%%%%%%%
	%% Optional Outputs %%
	%%%%%%%%%%%%%%%%%%%%%%

	if nargout >= 5
		% Check to see if E/B_w or G/C_v exist
		% If not just error out.

		if ~isfield(sys,'E') & ~isfield(sys,'B_w')
			error('E or B_w does not exist!')
		end

		if isfield(sys,'B_w')
			E = sys.B_w;
		else
			E = sys.E;
		end

		%Create E_bar
		E_at_each_n = {};
		for i = 1:T
			E_at_each_n{i} = E;
		end
		E_bar = blkdiag(E_at_each_n{:});

	end

	if nargout >= 6

		if ~isfield(sys,'G') & ~isfield(sys,'C_v')
			error('G or C_v does not exist!')
		end

		if isfield(sys,'C_v')
			C_v = sys.C_v;
		else
			C_v = sys.G;
		end

		%Create G_bar
		G_at_each_n = {};
		for i = 1:T
			if any(m_locs == (i-1))
				G_at_each_n{i} = zeros(size(C_v));
			else
				G_at_each_n{i} = C_v;
			end 
		end
		G_bar = blkdiag(G_at_each_n{:});	

	end

	%%%%%%%%%%%%%%%%%%%%
	%% Create Outputs %%
	%%%%%%%%%%%%%%%%%%%%

	varargout{1} = G;
	varargout{2} = H;
	varargout{3} = C_big;
	varargout{4} = x0_mat;

	if nargout >= 5
		varargout{5} = E_bar;
	end

	if nargout >= 6
		varargout{6} = G_bar;
	end
end
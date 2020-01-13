function [ spacecraft_lcsas ] = get_spacecraft_lcsas1(varargin)
%Description:
	%	Creates the dynamics for a quadrotor helicopter modeled after the Ascending Technologies' ResearchPilot.
	%	System should be a 10 dimensional.
	%
	%
	%Usage:
	%	ad = get_spacecraft_lcsas1()
	%	ad = get_spacecraft_lcsas1('dt',0.1)
	%	ad = get_spacecraft_lcsas1('dt',1)
	%	ad = get_spacecraft_lcsas1('disturb_info',P_w,P_v)
	%	ad = get_spacecraft_lcsas1('disturb_info',eta_w,eta_v)
	
	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	arg_idx = 1;
	while arg_idx <= nargin
		switch varargin{arg_idx}
			case 'dt'
				dt = varargin{arg_idx+1};
				arg_idx = arg_idx + 2;
			case 'disturb_info'
				if isa(varargin{arg_idx+1},'Polyhedron')
					P_w = varargin{arg_idx+1};
					P_v = varargin{arg_idx+2};
					arg_idx = arg_idx + 3;
				elseif isscalar(varargin{arg_idx})
					eta_w = varargin{arg_idx+1};
					eta_v = varargin{arg_idx+2};

					arg_idx = arg_idx + 3;

					%Construct Polyhedron with this.
					n_w = 2; n_v = 6;
					P_w = Polyhedron('lb',-eta_w*ones(1,n_w),'ub',eta_w*ones(1,n_w));
					P_v = Polyhedron('lb',-eta_v*ones(1,n_v),'ub',eta_v*ones(1,n_v));

				else
					error('Unexpected input for disturb_info.')
				end
			case 'L'
				L = varargin{arg_idx+1};
				if ~isa(L,'Language')
					error('Expected for L input to be a Language object.')
				end
				arg_idx = arg_idx + 2;
			otherwise
				error('It looks like this number of arguments isn''t currently supported. :(')
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Constants (Copied from Dustin Webb's Github) %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ~exist('dt')
		dt = 1;
	end

	%state
	%
	%	  [    x    ]
	% 	  [ \dot{x} ]
	% x = [   v_x   ]
	%	  [    y    ]
	%	  [ \dot{y} ]
	%	  [   v_y   ]
	n_x = 6;
	n_u = 2;
	n_w = 2;

	A = zeros(n_x);
	A(1,2) = 1;
	A(2,3) = 1;
	A(3,3) = -0.5;
	A(4,5) = 1;
	A(5,6) = 1;
	A(6,6) = -0.5;

	B = zeros(n_x,n_u);
	B(3,:) = [ 1 , 0.1 ];
	B(6,:) = [ 0.1 , 1 ];

	B_w = zeros(n_x,n_w);
	B_w(2,1) = 1; B_w(5,2) = 1;
    
    C = [1,zeros(1,n_x-1); 0,1,zeros(1,n_x-2) ];
    n_v = size(C,1);

	if ~exist('P_w')
		eta_w = 0.3; eta_v = 0.1;
		Pw0 = Polyhedron('lb',-eta_w*ones(1,n_w),'ub',eta_w*ones(1,n_w));
		Pw1 = Pw0 + eta_w*[1;0];
		Pw2 = Pw0 + eta_w*[0;1];
		Pv = Polyhedron('lb',-eta_v*ones(1,n_v),'ub',eta_v*ones(1,n_v));
	end

	% Discretize elements
	sys1 = ss(A,B,eye(n_x),0);
	dsys1 = c2d(sys1,dt);

	sys2 = ss(A,B_w,eye(n_x),0);
	dsys2 = c2d(sys2,dt);

	A_d = dsys1.A;
	B_d = dsys1.B;
	Bw_d = dsys2.B;

	% Create first Aff_Dyn

	ad0 = Aff_Dyn(	A_d,B_d,zeros(n_x,1),C,...
					Pw0,Pv, ...
					Bw_d,eye(n_v));

	ad1 = Aff_Dyn(  A_d,B_d,zeros(n_x,1),C,...
					Pw1,Pv, ...
					Bw_d,eye(n_v));

	ad2 = Aff_Dyn(  A_d,B_d,zeros(n_x,1),C,...
					Pw2,Pv, ...
					Bw_d,eye(n_v));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Creating the Mode Language L %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ~exist('L')
		L = Language([1,2,2,1,1,1],[1,3,3,1,1,1],[1,2,3,1,1,1]);
	end

	%%%%%%%%%%%%%%%%%%
	%% Create LCSAS %%
	%%%%%%%%%%%%%%%%%%

	spacecraft_lcsas = LCSAS([ad0,ad1,ad2],L);

end
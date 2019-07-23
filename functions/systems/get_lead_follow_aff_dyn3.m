function [dyn_cont, dyn_disc,foll_offsets] = get_lead_follow_aff_dyn3(varargin)
	%get_lead_follow_aff_dyn3.m
	%Summary:
	%	Returns the continuous and discretized dynamics of the leader-follower system
	%	proposed by Necmiye in a meeting on June 27.
	%	This version differs in that:
	%	- It designs its states such that all are designed to track/have an origin that is centered on the lead agent (#1).
	%	- The controller is meant to design controls for the followers, not the leader.
	%
	%Usage:
	%	[dyn_cont,dyn_disc] = get_lead_follow_aff_dyn3(xdim,ydim,dt)
	%	[dyn_cont,dyn_disc] = get_lead_follow_aff_dyn3(xdim,ydim,dt, 'disturb_info',eta_w, eta_v)
	%	[dyn_cont,dyn_disc] = get_lead_follow_aff_dyn3(xdim,ydim,dt, 'cube_params',r, K)



	%% Input Processing
	allowable_nargin = [3 6];

	if any(allowable_nargin == nargin)
		foll_cube.dim_x = varargin{1};
		foll_cube.dim_y = varargin{2};
		delta_t = varargin{3};

		if nargin > 3
			scrolling_idx = 4;
			while scrolling_idx < nargin	
				switch varargin{scrolling_idx}
				case 'disturb_info'
					eta_w = varargin{scrolling_idx+1};
					eta_v = varargin{scrolling_idx+2};
					scrolling_idx = scrolling_idx + 3;
				case 'cube_params'
					r = varargin{scrolling_idx+1};
					K = varargin{scrolling_idx+2};
					scrolling_idx = scrolling_idx + 3;
				otherwise
					error(['Input #' num2str(scrolling_idx) 'not recognized by function.'])
				end
			end
		end
		if ~exist('eta_w')
			eta_w = 0.2;
			eta_v = 0.1;
		end
		if ~exist('r')
			r = 2;
			K = 10;
		end


	else
		error(['Improper number of input arguments. Expecting any of these numbers for nargin: ' num2str(allowable_nargin)])
	end

	%% Constants

	%Follower Cube Parameters
	foll_cube.n = foll_cube.dim_x*foll_cube.dim_y;

	%Create desired "offsets". That is, if the leader is considered to be at the origin,
	%then these coordinates are the locations of the followers in the cube.
	foll_cube.offsets = [];
	for cube_idx = 1:foll_cube.n
		row_ind = floor(cube_idx/foll_cube.dim_x)+1;
		col_ind = rem(cube_idx,foll_cube.dim_x);
		if col_ind == 0
			col_ind = foll_cube.dim_x;
		end
		foll_cube.offsets(1,cube_idx) = [ -r*(col_ind) ];
	end
	foll_cube.offsets = [foll_cube.offsets; kron(linspace(r/2,-r/2,foll_cube.dim_y),ones(1,foll_cube.dim_x))];

	%% Create Continuous Time System
	n = (foll_cube.n+1)*2;
	p = n;
	A = zeros(n);

	subB = [];
	for block_idx = 1:foll_cube.dim_y
		for in_row_idx = 1:foll_cube.dim_x
			if in_row_idx == 2
				subB(1+(block_idx-1)*foll_cube.dim_x+in_row_idx,1) = -1;
			elseif in_row_idx > 2
				subB(1+(block_idx-1)*foll_cube.dim_x+in_row_idx,(block_idx-1)*foll_cube.dim_x+in_row_idx) = -1;
			end
			%Always add this.
			subB(1+(block_idx-1)*foll_cube.dim_x+in_row_idx,(block_idx-1)*foll_cube.dim_x+in_row_idx) = 1;
		end
	end

	B = kron(eye(2),subB);

	%Create disturbance matrix by using information about where the rows begin and end.
	B_w = zeros(n,2);
	B_w(1,1) = 1; B_w(1+foll_cube.n+1,2) = 1;
	for row_idx = 1:foll_cube.dim_y
		B_w(1+(row_idx-1)*foll_cube.dim_x+1,1) = -1;
		B_w(1+foll_cube.n+1+(row_idx-1)*foll_cube.dim_x+1,2) = -1;
	end

	%% Create Discrete Time System
	temp_sys1 = ss(A,B,eye(n),0);
	temp_sys3 = ss(A,B_w,eye(n),0);

	temp_dsys1 = c2d(temp_sys1,delta_t);
	temp_dsys3 = c2d(temp_sys3,delta_t);

	disc_aff_dyn = Aff_Dyn(temp_dsys1.A,temp_dsys1.B,zeros(n,1),eye(n),...
							eta_w,eta_v,...
							temp_dsys3.B,eye(n));

	%% Outputs
	dyn_cont = Aff_Dyn(A,B,zeros(n,1),eye(n),eta_w/delta_t,eta_v/delta_t,B_w,eye(p));
	dyn_disc = disc_aff_dyn;

	foll_offsets = foll_cube.offsets;

end
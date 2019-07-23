function [dyn_cont, dyn_disc,foll_offsets] = get_lead_follow_aff_dyn2(varargin)
	%get_lead_follow_aff_dyn2.m
	%Summary:
	%	Returns the continuous and discretized dynamics of the leader-follower system
	%	proposed by Necmiye in a meeting on June 27.
	%	This version focuses on ERROR states instead of the true states which is what is
	%	considered in get_lead_follow_aff_dyn.m.
	%Usage:
	%	[dyn_cont,dyn_disc] = get_lead_follow_aff_dyn(xdim,ydim,dt)
	%	[dyn_cont,dyn_disc] = get_lead_follow_aff_dyn(xdim,ydim,dt, 'disturb_info',eta_w, eta_v)
	%	[dyn_cont,dyn_disc] = get_lead_follow_aff_dyn(xdim,ydim,dt, 'cube_params',r, K)



	%% Input Processing
	allowable_nargin = [3 6];

	if any(allowable_nargin == nargin)
		xdim = varargin{1};
		ydim = varargin{2};
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
		else
			eta_w = 0.2;
			eta_v = 0.1;

			r = 2;
			K = 10;
		end


	else
		error(['Improper number of input arguments. Expecting any of these numbers for nargin: ' num2str(allowable_nargin)])
	end

	%% Constants

	%Follower Cube Parameters
	foll_cube.dim_x = xdim;
	foll_cube.dim_y = ydim;
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
	subA = [];
	subA = [ zeros(1,foll_cube.dim_x+1);
			[K*eye(foll_cube.dim_x), zeros(foll_cube.dim_x,1)] + [zeros(foll_cube.dim_x,1),-K*eye(foll_cube.dim_x) ] ]; 

	for block_ind = 2:foll_cube.dim_y
		for in_row_ind = 1:foll_cube.dim_x
			if in_row_ind == 1
				subA(1+(block_ind-1)*foll_cube.dim_x+in_row_ind,1) = K;
				subA(1+(block_ind-1)*foll_cube.dim_x+in_row_ind,(block_ind-1)*foll_cube.dim_x+in_row_ind+1) = -K;
			else
				subA(1+(block_ind-1)*foll_cube.dim_x+in_row_ind,(block_ind-1)*foll_cube.dim_x+[in_row_ind:in_row_ind+1]) = [K,-K];
			end
		end
	end

	A = kron(eye(2),subA);
	n = size(A,1);

	B = [[1;-ones(foll_cube.n,1);zeros(1+(foll_cube.n),1)],[zeros(1+(foll_cube.n),1);1;-ones(foll_cube.n,1)]];

	B_w = eye(2*(1+foll_cube.n));

	node_hghts = K*fliplr(linspace(-r/2,r/2,foll_cube.dim_y))';
	%F = [0;-K*r*ones(foll_cube.n,1); 0; K*(r/2);0;-K*(r/2);0]; %Specialized for 2 x 2 case
	F = [0;-K*r*ones(foll_cube.n,1); 0;kron(node_hghts,[1;zeros(foll_cube.dim_x-1,1)])];

	%Used desired offsets to update the F matrix
	dF = [];
	for cube_idx = 1:foll_cube.n
		row_ind = floor(cube_idx/foll_cube.dim_x)+1;
		col_ind = rem(cube_idx,foll_cube.dim_x);
		if col_ind == 0
			col_ind = foll_cube.dim_x;
		end

		if col_ind == 1
			dF(1+cube_idx,1) = -K*foll_cube.offsets(1,cube_idx);
			dF(1+foll_cube.n+1+cube_idx,1) = -K*foll_cube.offsets(2,cube_idx);
		else
			dF(1+cube_idx,1) = -K*(foll_cube.offsets(1,cube_idx) - foll_cube.offsets(1,cube_idx-1));
			dF(1+foll_cube.n+1+cube_idx,1) = -K*(foll_cube.offsets(2,cube_idx) - foll_cube.offsets(2,cube_idx-1));
		end
	end

	F = F + dF;

	%% Create Discrete Time System
	temp_sys1 = ss(A,B,eye(n),0);
	temp_sys2 = ss(A,F,eye(n),0);
	temp_sys3 = ss(A,B_w,eye(n),0);

	temp_dsys1 = c2d(temp_sys1,delta_t);
	temp_dsys2 = c2d(temp_sys2,delta_t);
	temp_dsys3 = c2d(temp_sys3,delta_t);

	disc_aff_dyn = Aff_Dyn(temp_dsys1.A,temp_dsys1.B,temp_dsys2.B,eye(n),...
							eta_w,eta_v,...
							temp_dsys3.B,eye(n));

	%% Outputs
	dyn_cont = Aff_Dyn(A,B,F,eye(n),eta_w/delta_t,eta_v/delta_t);
	dyn_disc = disc_aff_dyn;

	foll_offsets = foll_cube.offsets;

end
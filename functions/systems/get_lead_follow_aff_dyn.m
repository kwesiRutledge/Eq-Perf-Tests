function [dyn_cont, dyn_disc] = get_lead_follow_aff_dyn(varargin)
	%get_lead_follow_aff_dyn.m
	%Summary:
	%	Returns the continuous and discretized dynamics of the leader-follower system
	%	proposed by Necmiye in a meeting on June 27.
	%Usage:
	%	[dyn_cont,dyn_disc] = get_lead_follow_aff_dyn(xdim,ydim,dt)
	%	[dyn_cont,dyn_disc] = get_lead_follow_aff_dyn(xdim,ydim,dt, eta_w, eta_v)



	%% Input Processing
	allowable_nargin = [3 5];

	if any(allowable_nargin == nargin)
		xdim = varargin{1};
		ydim = varargin{2};
		delta_t = varargin{3};
		switch nargin
		case 5
			eta_w = varargin{4};
			eta_v = varargin{5};
		otherwise
			eta_w = 0.2;
			eta_v = 0.1;
		end


	else
		error(['Improper number of input arguments. Expecting any of these numbers for nargin: ' num2str(allowable_nargin)])
	end

	%% Constants

	%Follower Cube Parameters
	foll_cube.dim_x = xdim;
	foll_cube.dim_y = ydim;
	foll_cube.n = foll_cube.dim_x*foll_cube.dim_y;

	r = 2;

	K = 10;

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

	B = [[1;zeros(1+2*(foll_cube.n),1)],[zeros(1+(foll_cube.n),1);1;zeros((foll_cube.n),1)]];

	B_w = eye(2*(1+foll_cube.n));

	%Specialized to 2x2 case
	node_hghts = K*fliplr(linspace(-r/2,r/2,foll_cube.dim_y))';
	%F = [0;-K*r*ones(foll_cube.n,1); 0; K*(r/2);0;-K*(r/2);0];
	F = [0;-K*r*ones(foll_cube.n,1); 0;kron(node_hghts,[1;zeros(foll_cube.dim_x-1,1)])];

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

end
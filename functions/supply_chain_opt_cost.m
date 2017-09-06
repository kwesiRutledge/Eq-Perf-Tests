function [ cost ] = supply_chain_opt_cost( varargin )
	%Calculate the cost for this inventory management problem.

	%Process inputs
	is_final_cost = varargin{1};
	x 			  = varargin{2};
	
	if is_final_cost
		%do nothing. No input needed.
	else
		u = varargin{3};
	end

	if nargin > 3
		error('Too many inputs given.')
	end

	%Constants
	s = 8;
	h = 1;
	b = 4;
	g = 4;

	%Calculate
	if is_final_cost
		cost = max( [g*x -b*x] );
	else
		cost = max( [h*x -b*x] ) + s*abs(u);
	end
end
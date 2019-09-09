function [Lambda,constrs] = get_H_polyt_inclusion_constr( obj , H_x , h_x , H_y, h_y )
	%Description:
	%	Provides the constraints and optimization variables needed to define the constraint: \mathbb{X} \subseteq \mathbb{Y}
	%	where \mathbb{X} = \{ x | H_x * x <= h_x } and \mathbb{Y} = \{ y | H_y * y <= h_y }
	%
	%Usage:
	%
	%Inputs:
	%	H_x 	- Hyperplane representation of the x-polytope's normal vectors
	%	h_x 	- Hyperplane representaiton of the X Polytope's offset values
	%	H_y 	- Etc.
	%	h_y 	- Etc.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(H_x,2);

	q_x = size(H_x,1);
	q_y = size(H_y,1);

	if n ~= size(H_y,2)
		error(['The two polytopes appear to be in different dimension. X in dimension ' num2str(n) ', while Y in dimension ' num2str(size(H_y,2))])
	end

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	Lambda = sdpvar(q_y,q_x,'full');

	constrs = [];
	constrs = constrs + [Lambda >= 0]; %Lambda must be nonnegative
	constrs = constrs + [ Lambda * H_x == H_y];
	constrs = constrs + [ Lambda *h_x <= h_y];

end
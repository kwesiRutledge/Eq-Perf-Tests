function [ constraints_out , yalmip_vars_out ] = no_other_belief_sequences_reached( cg , all_paths , active_path_flags , lcsas_in , K_test , k_test , x0 )
	%Description:
	%	This function receives the set of all reasonable belief sequences for the system, the set of active paths, and the input lcsas
	%	and attempts to determine if the controller with gains described by K_test and k_test will never lead to beliefs that are inactive (i.e.
	%	beliefs corresponding to active_path_flags(i) = false .)
	%
	%Usage:
	%	[ constrain_other_beliefs_empty , y_arr ] = cg.no_other_belief_sequences_reached( all_paths , active_path_flags , lcsas_in , K , k )
	%

	%% Input Processing %%

	if size(all_paths,2) ~= length(active_path_flags)
		error(['The number of active_path_flags (' num2str(length(active_path_flags)) ') should be equal to the number of paths (' num2str(size(all_paths,2)) ').'])
	end

	if isempty(lcsas_in.X0)
		error('The X0 field in lcsas_in should be nonempty for this function.')
	end

	if ~isa(lcsas_in.X0,'Polyhedron')
		error(['Expected X0 to be a Polyhedron, but it is of class ' class(lcsas_in.X0) ])
	end

	if isempty(lcsas_in.U)
		error('The U field in lcsas_in should be nonempty for this function.')
	end

	if ~isa(lcsas_in.U,'Polyhedron')
		error(['Expected U to be a Polyhedron, but it is of class ' class(lcsas_in.U) ])
	end

	%% Constants

	eps0 = 10^(-4);

	Px0 = lcsas_in.X0;
	Pu  = lcsas_in.U;

	TimeHorizon = size(all_paths,1);

	inactive_path_indices = find( active_path_flags == false )'; %Need this to be a row vector for some reason.

	%% Algorithm

	y = {}; %Create dual YALMIP Variables
	constraints_out = [];

	for inactive_path_index = inactive_path_indices

		% Extract Path
		current_path = all_paths(:,inactive_path_index);

		%Get the closed_loop_consistent_internal_behavior_set with respect to
		%this path when there is controller defined by (K_test,k_test) applied.

		[ H_PhiI , h_PhiI ] = lcsas_in.get_consistent_internal_behavior_matrices( TimeHorizon , current_path(end) , Pu , Px0 );

		%Get Closed Loop Consistent Internal Behavior Set Matrices (Function of K)
		[ H_cl , h_cl ] = lcsas_in.get_closed_loop_consistent_internal_behavior_set_matrices( ...
			H_PhiI , h_PhiI , ...
			x0 , ...
			K_test , k_test , current_path );

		% Construct Empty Polyhedron Constraint
		y{end+1} = sdpvar(size(H_cl,1),1,'full');
		constraints_out = constraints_out + ...
			[ H_cl'*y{end} == 0] + [ h_cl'*y{end} <= -eps0 ] + [y{end} >= 0];


	end

	% Create outputs
	yalmip_vars_out = y;

end
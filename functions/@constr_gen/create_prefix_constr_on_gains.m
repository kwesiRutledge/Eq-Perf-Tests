function [constraints] = create_prefix_constr_on_gains( obj , lcsas , Q_set , r_set )
	%create_sadraddini_cond.m
	%Description:
	%	This function should use the language in L to constrain the gains defined by
	%	Q_set and r_set.
	%
	%	Note that it is assumed that the word L{i} is associated with gains Q_set{i} and r_set{i}
	%
	%Usage:
	%	[constraints] = cg.create_prefix_constr_on_gains( lcsas , L , Q_set , r_set )
	%Inputs:
	%

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if ~isa(lcsas,'LCSAS')
		error('The second input must be a LCSAS object.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	n = size(lcsas.Dyn(1).A,1);
	m = size(lcsas.Dyn(1).B,2);
	p = size(lcsas.Dyn(1).C,1);
	wd = size(lcsas.Dyn(1).B_w,2); %Assume that the size of the disturbance doesn't change for any of the included dynamics
	vd = size(lcsas.Dyn(1).C_v,2);

	L = lcsas.L;

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	constraints = [];

	%Iterate through every path
	for path_idx1 = 1:length(L.words)
		for path_idx2 = path_idx1+1:length(L.words)
			path1 = L.words{path_idx1};
			path2 = L.words{path_idx2};

			%Truncate one if necessary.
			if length(path1) < length(path2)
				path2 = path2(1:length(path1));
			elseif length(path1) > length(path2)
				path1 = path1(1:length(path2));
			end

			%Check the length of the matching prefixes
			match_vec = ~(path1 == path2);
			
			if isempty(find(match_vec,1))
				%If the words are completely the same, then return the length of either word.
				shared_pref_length = length(path1);
			else
				shared_pref_length = find(match_vec,1)-1;
			end

			%Constrain the matrices appropriately.
			if shared_pref_length > 0
				constraints = constraints +  [Q_set{path_idx1}( [1:shared_pref_length*m] , [1:shared_pref_length*p] ) == Q_set{path_idx2}( [1:shared_pref_length*m] , [1:shared_pref_length*p] )];
				constraints = constraints +  [r_set{path_idx1}( [1:shared_pref_length*m] ) == r_set{path_idx2}([1:shared_pref_length*m]) ];
			end
		end
	end
end
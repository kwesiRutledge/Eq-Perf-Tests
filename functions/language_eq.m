function [l_eq] = language_eq(L1,L2)
% language_eq.m
%	Description:
%		Returns true if the set defined by the cell matrix L1 is the same as the set defined by L2.
%	
%	Inputs:
%		L1 - Cell array of numerical matrices
%		L2 - Cell array of numerical matrices

%% Constants

%% Verify L1 is a subset of L2

L1_vals_found = false(length(L1),1);
for L1_ind = 1:length(L1)
	%Verify that each
	for L2_ind = 1:length(L2)
		if (length(L1{L1_ind}) == length(L2{L2_ind}))
			L1_vals_found(L1_ind) = L1_vals_found(L1_ind) || all(L1{L1_ind} == L2{L2_ind});
		end
	end
end

%% Verify L2 is a subset of L1
L2_vals_found = false(length(L2),1);
for L2_ind = 1:length(L2)
	%Verify that each
	for L1_ind = 1:length(L1)
		if (length(L1{L1_ind}) == length(L2{L2_ind}))
			L2_vals_found(L2_ind) = L2_vals_found(L2_ind) || all(L1{L1_ind} == L2{L2_ind});
		end
	end
end

l_eq = all(L1_vals_found) && all(L2_vals_found);

end
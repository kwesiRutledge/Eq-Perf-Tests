function [ MatchingSetOfPaths ] = FindPathsWithPrefix( SetOfPaths , Prefix )
	%Description:
	%
	%Inputs:
	%	SetOfPaths - A Matrix of Language objects.
	%				 Each column is a path (a sequence of Languages).
	%				 The number of columns is the number of paths.
	%	Prefix - A vector of Language objects.
	%			 Should be smaller than the length of the paths.

	%% Constants %%

	[Tp1, NumPaths] = size(SetOfPaths);
	T = Tp1 - 1;

	PrefixLength = length(Prefix);

	%% Input Processing %%

	if length(Prefix) > Tp1
		error(['The input Prefix to FindPathsWithPrefix() is of length ' num2str(length(Prefix)) ' but the paths from SetOfPaths have length ' num2str(Tp1) '.' ])
	end

	%% Algorithm %%
	MatchingSetOfPaths = [];

	for path_index = 1:NumPaths
		TempPath = SetOfPaths(:,path_index);

		%disp(TempPath(1:PrefixLength) == Prefix)

		if TempPath(1:PrefixLength) == Prefix
			MatchingSetOfPaths = [ MatchingSetOfPaths , TempPath ];
		end

	end


end
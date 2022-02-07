function [tf] = IsProperlyBranching( SetOfPaths , L )
	%Description:
	%	Identifies if the set of paths given as SetOfPaths is properly branching according to the language L.

	%% Constants %%

	[Tp1, NumPaths] = size(SetOfPaths);
	T = Tp1 - 1;

	%% Algorithm %%

	tf = true;

	for path_index = 1:NumPaths
		temp_path = SetOfPaths(:,path_index); %Get A Path

		%At each time, evaluate if there is a branch.
		for t = 1:T
			branchDetected = temp_path(t).cardinality() > temp_path(t+1).cardinality();

			if branchDetected
				LangAt_t = temp_path(t);
				PathsWithSharedPrefix = SetOfPaths.FindPathsWithPrefix( temp_path([1:t]) );

				if isempty(PathsWithSharedPrefix)
					error('What the--')
				end

				%Validate that at the branching time, there is another path which can be used to make up L.
				LangsAt_tp1 = PathsWithSharedPrefix(t+1,:);
                switch length(LangsAt_tp1)
                    case 1
                        if LangsAt_tp1 ~= L
                            tf = false;
                            return;
                        end
                    otherwise
                        if (LangsAt_tp1(1).union(LangsAt_tp1([2:end]))) ~= LangAt_t
                            tf = false;
                            return
                        end
                end

			end

		end


	end



end
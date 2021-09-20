function tf = PathSubsetCoversBehaviorsFromWord( PathsSubset , Word , System )
	%Description:
	%
	%Usage:
	%	tf = Paths.PathSubsetCoversBehaviorsFromWord( Word , System )

	%% Input Processing %%

	if ~isa(System,'LCSAS')
		error(['System input to PathSubsetCoversBehaviorsFromWord() is of class ' class(System) ' but expected class LCSAS.' ])
	end

	%% Constants %%

	[T,NumPaths] = size(PathsSubset); 

	%% Algorithm %%

	% Create InternalBehaviorSets.
	ibs_array = [];
	for path_index = 1:NumPaths
		ibs_array = [ibs_array; InternalBehaviorSet( System , PathsSubset(:,path_index) ) ];
	end

	[tf] = ibs_array.CoversWordBehaviors( Word );

end
function [ possible_subsets_of_paths , subsets_as_binary_flags1 ] = organize_subsets_of_paths( unchecked_knowledge_sequences , subsets_as_binary_flags0 , organization_strategy )
	%Description:
	%	This function organizes the subsets of unchecked_knowledge_sequences according to some rule.
	%	Some valid rules or "organization_strategies" are:
	%	- "None"

	%% Constants %%

	num_unchecked_sequences = size(unchecked_knowledge_sequences,2);

	%% Algorithm %%

	switch organization_strategy
	case {'None','none','AscendingCardinality'}
		%Do nothing.
		%Return everything that was input to the function.
		subsets_as_binary_flags1 = subsets_as_binary_flags0;

		possible_subsets_of_paths = {};
		for subset_index = 1:length(subsets_as_binary_flags1)
			possible_subsets_of_paths{subset_index} = unchecked_knowledge_sequences(:,subsets_as_binary_flags1{subset_index});
		end

	case {'DescendingCardinality'}

		subsets_as_binary_flags1 = fliplr(subsets_as_binary_flags0);


		possible_subsets_of_paths = {};
		for subset_index = 1:length(subsets_as_binary_flags1)
			possible_subsets_of_paths{subset_index} = unchecked_knowledge_sequences(:,subsets_as_binary_flags1{subset_index});
		end

	case {'DescendingCardinality+MostlyIgnoreAllModeSequence'}

		temp_subsets_as_binary_flags = fliplr(subsets_as_binary_flags0);

		top_half = {}; bottom_half = {};
		for subset_index = 1:length(temp_subsets_as_binary_flags)
			temp_subset = temp_subsets_as_binary_flags{subset_index};
			if temp_subset(end) == 1
				bottom_half{end+1} = temp_subset;
			else
				top_half{end+1} = temp_subset;
			end
		end

		subsets_as_binary_flags1 = { top_half{:} , bottom_half{:} };

		possible_subsets_of_paths = {};
		for subset_index = 1:length(subsets_as_binary_flags1)
			possible_subsets_of_paths{subset_index} = unchecked_knowledge_sequences(:,subsets_as_binary_flags1{subset_index});
		end

	case {'DescendingCardinality+PreferAllModeSequence'}

		temp_subsets_as_binary_flags = fliplr(subsets_as_binary_flags0);

		top_half = {}; bottom_half = {};
		for subset_index = 1:length(temp_subsets_as_binary_flags)
			temp_subset = temp_subsets_as_binary_flags{subset_index};
			if temp_subset(end) == 1
				top_half{end+1} = temp_subset;
			else
				bottom_half{end+1} = temp_subset;
			end
		end

		subsets_as_binary_flags1 = { top_half{:} , bottom_half{:} };

		possible_subsets_of_paths = {};
		for subset_index = 1:length(subsets_as_binary_flags1)
			possible_subsets_of_paths{subset_index} = unchecked_knowledge_sequences(:,subsets_as_binary_flags1{subset_index});
		end


	otherwise
		error(['Unexpected organization_strategy given: ' organization_strategy ])

	end

end
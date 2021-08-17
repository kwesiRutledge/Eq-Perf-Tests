function [node_out] = GetNodeWithHighestCardinality( node_array )
	%Description:
	%	Self-described?

	%% Constants

	%% Algorithm
	max_card = -1;
	max_index = -1;

	for node_index = 1:length(node_array)
		node_i = node_array(node_index);
		card_i = node_i.subL.cardinality();

		if card_i > max_card
			max_card = card_i;
			max_index = node_index;
		end
	end

	node_out = node_array(node_index);

end
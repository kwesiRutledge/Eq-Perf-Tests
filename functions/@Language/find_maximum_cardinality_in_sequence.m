function [ max_card , max_card_index ] = find_maximum_cardinality_in_sequence( KnowledgeSequenceIn )
	%Description:
	%	This function searches through the elements of KnowledgeSequenceIn to find the Language with the
	%	with the largest cardinality.
	%
	%Usage:
	%	max_card = KnowledgeSequence.find_maximum_cardinality_in_sequence()
	
	%% Input Processing %%

	if ~isvector(KnowledgeSequenceIn)
		error('Input to the function find_maximum_cardinality_in_sequence() must be a vector.')
	end

	%% Algorithm
	max_card = -1;
	max_card_index = -1;

	for lang_index = 1:length(KnowledgeSequenceIn)

		temp_lang = KnowledgeSequenceIn(lang_index);

		if temp_lang.cardinality() > max_card
			max_card_index = lang_index;
			max_card = temp_lang.cardinality();
		end


	end


end
function [ tf , matching_lang_idcs ] = prefix_in_list_of_languages( prefix_in , list_of_langs )
	%Description:
	%	The function returns true if the prefix given (prefix_in)
	%	exists in one of the languages in list_of_langs.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	if ~isa(list_of_langs(1),'Language')
		error('The input array is not an array of Language.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	tf = false;
	matching_lang_idcs = [];
	for lang_idx = 1:length(list_of_langs)
		temp_lang = list_of_langs(lang_idx);

		for word_idx = 1:temp_lang.cardinality()
			temp_word = temp_lang.words{word_idx};

			if length(temp_word) >= length(prefix_in)
				if all(prefix_in == temp_word([1:length(prefix_in)]))
					tf = true;
					matching_lang_idcs = [matching_lang_idcs,lang_idx];
				end
			end
		end
	end

end
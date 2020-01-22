function [matching_idcs, matching_words] = find_all_words_with_pref( obj, pref )
	%Language.find_all_words_with_pref()
	%Description:
	%	Identifies the words and the indices of any word in the current language that
	%	has the prefix given.
	%
	%Usage:
	%	[idx,matching_word] = L.find_all_words_with_pref(pref)
	%	[idcs,matching_words] = L.find_all_words_with_pref(pref)

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	pref_len = length(pref);

	if pref_len == 0
		error('Given an empty word for prefix!')
	end

	if pref_len > obj.find_longest_length()
		error('Given prefix is longer than all words in Language.')
	end

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	%Match pref with a word from L
	matching_mat = [];
	for word_ind = 1:length(obj.words)
		if length(obj.words{word_ind}) >= pref_len
			matching_mat = [matching_mat; obj.words{word_ind}(1:pref_len)];
		else
			%If the word is not long enough to match, then place some nonsense in the
			%corresponding row of the matching matrix.
			matching_mat = [matching_mat; Inf(1,pref_len)];
		end
	end
	matching_mat = repmat(pref,length(obj.words),1) == matching_mat;	
    
    %Showing which words match the observed prefix
    if pref_len == 1
        matching_locs = matching_mat;
    elseif pref_len > 1
        matching_locs = all(matching_mat')';
    else
        error('pref_len is a nonpositive integer.')
    end

    temp_idcs = [1:length(matching_locs)];
	matching_idcs = temp_idcs(matching_locs);

	matching_words = Language({obj.words{matching_idcs}});

end
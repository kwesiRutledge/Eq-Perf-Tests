function subset_flag = subseteq(obj,L_in)
			%Description:
			%	Returns true if the Language 'L_in' is a subset of the current language.
			%	Otherwise returns false.
			%
			%Usage:
			%	in_lang_flag = L.contains(L_in)
			%	in_lang_flag = L.contains( Language([1,2,1,1,1]) )

			%% Variables

			subset_flag = false;

			num_words = obj.cardinality();

			%% Algorithm

			word_flags = false(1,num_words);

			%Test every word in the object 'obj'.
			for word_idx = 1:num_words
				temp_word = obj.words{word_idx};
				if L_in.contains( temp_word )
					word_flags(word_idx) = true;
				end
			end

			%obj is a subset or equal to L_in if and only if all of the word flags is now true.
			subset_flag = all(word_flags);

		end
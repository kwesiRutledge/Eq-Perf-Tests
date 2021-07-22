function subset_flag = subseteq(obj,L_in)
	%Description:
	%	Returns true if the Language 'L_in' is a subset of the current language.
	%	Otherwise returns false.
	%
	%Usage:
	%	in_lang_flag = L.contains(L_in)
	%	in_lang_flag = L.contains( Language([1,2,1,1,1]) )

	%% Input Processing %%

	subset_type = inputprocessing_subseteq(obj,L_in);

	%% Variables

	subset_flag = false;

	%% Algorithm

	switch subset_type
	case 'scalar'
		num_words = obj.cardinality();
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
	case 'vector'
		subset_flag = true;

		for language_index = 1:length(obj)
			temp_obj_lang = obj(language_index);
			temp_L_in_lang = L_in(language_index);

			temp_flag = temp_obj_lang.subseteq(temp_L_in_lang);
			subset_flag = subset_flag && temp_flag;
		end
	otherwise
		error(['Unexpected subset_type value: ' subset_type ])
	end

end

function [subset_type] = inputprocessing_subseteq(obj,L_in)
	%Description:
	%	

	if ~all( size(obj) == size(L_in) )
		error(['The size of obj is ' num2str(size(obj)) ' while the size of L_in is ' num2str(size(L_in)) '.'  ])
	end

	if isscalar(obj) && isscalar(L_in)
		subset_type = 'scalar';
	elseif isvector(obj) && isvector(L_in)
		subset_type = 'vector';

		if length(obj) ~= length(L_in)
			error(['Length of obj is ' num2str(length(obj)) ', while the length of L_in is ' num2str(length(L_in)) '.' ])
		end

	else
		error(['The sizes of input obj and L_in do not match!'])
	end

end
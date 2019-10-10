classdef Language
	%Description:
	%	This class is meant to contain the methods needed to model languages used in the Ozay group's prefix-based methods.
	%
	%Construction:
	%	L1 = Language([1,2],[3,4,5],[1,4,5])
	%	L2 = Language({[1,2,3],[2,3,3]})
	%
	%Methods:
	%	- contains
	%	- subseteq
	%	- union
	%	- powerset
	%	- is_eq
	%	- all_words_start_with_root
	%	- find_longest_length
	%	- find_shortest_length
	%	- find_a_word_with_pref
	%
	%Member Variables:
	%	words 	- A cell array of numeric arrays.
	%			  In general languages can contain words of different lengths, which is why this is a cell array.

	properties
		words;
	end

	methods
		function L = Language(varargin)
			%Description:
			%	Creates language by receiving each word as a numeric array.
			%
			%Examples:
			%	L1 = Language([1,2],[3,4,5],[1,4,5])

			L.words = {};


			if (nargin == 1) && iscell(varargin{1})
				%Check to see if each element of the cell array is a row vector.
				temp_c_arr = varargin{1};
				for word_idx = 1:length(temp_c_arr)
					if size(temp_c_arr,1) ~= 1
						error('All words in the language should be row vectors.')
					end
				end

				%Save
				L.words = temp_c_arr;
			else
				for word_idx = 1:nargin
					%Check to see if each argument is numeric.
					if ~isnumeric(varargin{word_idx})
						error('All inputs to a language must be numeric.')
					end

					if size(word_idx,1) ~= 1
						error('All words should be row vectors.')
					end

					L.words{word_idx} = varargin{word_idx};
				end
			end
		end

		function in_lang_flag = contains(obj,word_in)
			%Description:
			%	Returns true if the word 'word_in' is in the language's words/
			%	Otherwise returns false.
			%
			%Usage:
			%	in_lang_flag = L.contains([1,2,1,1,1])

			%% Variables

			in_lang_flag = false;

			%% Algorithm

			for lang_idx = 1:length(obj.words)
				if size(obj.words{lang_idx},2) == size(word_in,2)
					if all(obj.words{lang_idx} == word_in)
						in_lang_flag = true;
					end
				end
			end

		end

		function [L_union] = union(obj,langs_in)
			%Description:
			%
			%Inputs:
			%	langs_in 	- An array of Language objects.
			%				  Example: [ Language([1,2],[1,3]) , Language([1,5]) ]

			%% Input Checking

			if ~isa(langs_in(1),'Language')
				error('The inputs to this function should be a cell array of Language class objects.')
			end

			%% Variables

			L_union = obj;

			%% Algorithm

			for lang_idx = 1:length(langs_in)
				temp_lang = langs_in(lang_idx);
				%Iterate through every word in the language
				for word_idx = 1:length(temp_lang.words)
					temp_word = temp_lang.words{word_idx};

					%Add to union if it was not already there.
					if ~L_union.contains(temp_word)
						L_union.words{end+1} = temp_word;
					end
				end
			end

		end

		function powerset_L = powerset(obj)
			%Description:
			%	Computes an array of languages where each language in the array is composed of some subset of obj.

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			powerset_L = [];

			for comb_length = 1:length(obj.words)
				temp_combs = nchoosek([1:length(obj.words)],comb_length);
				for comb_ind = 1:size(temp_combs,1)
					%Create temporary language
					temp_lang = Language();
					temp_lang.words = { obj.words{temp_combs(comb_ind,:)} };
					% Update the powerset
					powerset_L = [powerset_L, temp_lang ];
				end
			end

		end

		function [l_eq_flag] = is_eq(obj,L_in)
			%Description:
			%	Verifies if every word in the current language is in L and vice versa.

			%% Verify L1 is a subset of L2

			obj_vals_found = false(length(obj.words),1);
			for L1_ind = 1:length(obj.words)
				%Verify that each
				for L2_ind = 1:length(L_in.words)
					if (length(obj.words{L1_ind}) == length(L_in.words{L2_ind}))
						obj_vals_found(L1_ind) = obj_vals_found(L1_ind) || all(obj.words{L1_ind} == L_in.words{L2_ind});
					end
				end
			end

			%% Verify L2 is a subset of L1
			L_vals_found = false(length(L_in.words),1);
			for L2_ind = 1:length(L_in.words)
				%Verify that each
				for L1_ind = 1:length(obj.words)
					if (length(obj.words{L1_ind}) == length(L_in.words{L2_ind}))
						L_vals_found(L2_ind) = L_vals_found(L2_ind) || all(obj.words{L1_ind} == L_in.words{L2_ind});
					end
				end
			end

			l_eq_flag = all(obj_vals_found) && all(L_vals_found);

		end

		function [reached_root] = all_words_start_with_root(obj)
			%Description:
			%	This function returns true if all words in the language begin with 1, the first node in a belief graph.

			reached_root = true;
			for word_idx = 1:length(obj.words)
				reached_root = reached_root && (obj.words{word_idx}(1) == 1);
			end
		end

		function longest_T = find_longest_length(obj)
			%Description:
			%	Searches through all elements of the words for this node and determines
			%	how long the longest word is.

			longest_T = -1;

			for L_idx = 1:length(obj.words)
				longest_T = max(length(obj.words{L_idx}),longest_T);
			end
		end

		function shortest_T = find_shortest_length(obj)
			%Description:
			%	Searches through all elements of the words for this node and determines
			%	how long the longest word is.
			%
			%Usage:
			%	shortest_T = Language.find_shortest_length()

			shortest_T = Inf;

			for L_idx = 1:length(obj.words)
				shortest_T = min(length(obj.words{L_idx}),shortest_T);
			end
		end

		function card = cardinality(obj)
			%Description:
			%	Returns the cardinality of the set.
			%
			%Assumption:
			%	Assumes every word in the language is unique.

			card = length(obj.words);

		end

	end

end

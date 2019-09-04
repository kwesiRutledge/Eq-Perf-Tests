classdef Language
	%Description:
	%	This class is meant to contain the methods needed to model languages used in the Ozay group's prefix-based methods.
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

			L.words = {}

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

		function in_lang_flag = contains(obj,word_in)
			%Description:
			%	Returns true if the word 'word_in' is in the language's words/
			%	Otherwise returns false.

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

	end

end

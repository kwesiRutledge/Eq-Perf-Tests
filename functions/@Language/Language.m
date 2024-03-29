classdef Language
	%Description:
	%	This class is meant to contain the methods needed to model languages used in the Ozay group's prefix-based methods.
	%
	%Construction:
	%	L1 = Language([1,2],[3,4,5],[1,4,5])
	%	L2 = Language({[1,2,3],[2,3,3]})
	%
	%Methods:
	%	- Language
	%	- contains
	%	- contains_prefix
	%	- subseteq
	%	- union
	%	- powerset
	%	- is_eq
	%	- eq (==)
	%	- ne (~=)
	%	- all_words_start_with_root
	%	- find_longest_length
	%	- find_shortest_length
	%	- cardinality
	%	- get_all_symbols_at_idx
	%	- langtostr
	%	- lang2str
	%	- find_a_word_with_pref
	%	- find_all_words_with_pref
	%	- set_minus
	%	- create_belief_sequences_of_length
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
			%	L2 = Language({[1,2],[3,4,5],[1,4,5]})

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

		function [in_lang_flag,pos_in_lang] = contains(obj,word_in)
			%Description:
			%	Returns true if the word 'word_in' is in the language's words/
			%	Otherwise returns false.
			%	The output pos_in_lang is the first position where word_in
			%	occurs in the Language L.
			%
			%Usage:
			%	in_lang_flag = L.contains([1,2,1,1,1])
			%	[in_lang_flag,pos_in_lang] = L.contains([1,2,1,1,1])

			%% Input Checking

			if isvector(obj) && ~isscalar(obj)
				%If there is not a single Language given as an input, but multiple
				%Call the function for each language in obj.
				[in_lang_flag,pos_in_lang] = obj.contains_array_input(word_in);
				return;
			end

			%% Variables

			in_lang_flag = false;
			pos_in_lang = -1;

			%% Algorithm

			for lang_idx = 1:length(obj.words)
				if size(obj.words{lang_idx},2) == size(word_in,2)
					if all(obj.words{lang_idx} == word_in)
						in_lang_flag = true;
						pos_in_lang = lang_idx;
						return;
					end
				end
			end

		end

		function [in_lang_flag,pos_in_lang] = contains_array_input(obj,word_in)
			%Description:
			%	Returns an array describing whether or not the word 'word_in' is 
			%	each of the languages in obj (obj is an array).
			%
			%Usage:
			%	in_lang_flag = L_arr.contains([1,2,1,1,1])

			%% Variables
			
			in_lang_flag = false(size(obj));
			pos_in_lang = -1*ones(size(obj));

			%% Algorithm

			for obj_index = 1:length(obj)
				obj_i = obj(obj_index);
				[ in_lang_flag(obj_index) , pos_in_lang(obj_index) ] = obj_i.contains(word_in);
			end

		end

		function [in_lang_flag,pos_in_lang] = contains_prefix(obj,prefix_in)
			%Description:
			%	Returns true if the word 'prefix_in' is in the language's words/
			%	Otherwise returns false.
			%
			%Usage:
			%	in_lang_flag = L.contains_prefix([1,2,1,1,1])
			%	[in_lang_flag, pos_in_lang] = L.contains_prefix([1,2,1,1,1])

			%% Input Processing

			%% Variables

			in_lang_flag = false;
			pos_in_lang = [];

			%% Algorithm for Vector obj

			% if isvector(obj)
			% 	%If this is a vector of languages, then search through each language.
			% 	in_lang_flag = false;
			% 	for obj_index = 1:length(obj)
			% 		temp_flag = obj(obj_index).contains_prefix(prefix_in);
			% 		if temp_flag
			% 			in_lang_flag = true;
			% 			pos_in_lang = [pos_in_lang,obj_index];
			% 		end
			% 	end
			% end

			%% Scalar Obj Algorithm

			for lang_idx = 1:obj.cardinality()
				if length(obj.words{lang_idx}) >= length(prefix_in)
					if all(obj.words{lang_idx}([1:length(prefix_in)]) == prefix_in)
						in_lang_flag = true;
						pos_in_lang = [ pos_in_lang , lang_idx ] ;
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

		function [powerset_L,powerset_word_indices] = powerset(obj)
			%Description:
			%	Computes an array of languages where each language in the array is composed of some subset of obj.
			%	Does not contain the emptyset which should normally be in a powerset definition.
			%
			%Outputs:
			%	powerset_L = An array of language objects where each language object is a subset
			%				 of the language object obj.
			%	powerset_word_indices = A cell array of indices corresponding to the word indexes
			%							that make up powerset_L. See example for explanation.
			%
			%Example:
			%	Consider the language L = Language([1,2,3,4],[5,6,5,8],[9,10,11,12])
			%	The results of [powerset_L,powerset_word_indices] = L.powerset() should be the following
			%		powerset_L = [ Language([1,2,3,4]) , Language([5,6,7,8]) , Language([9,10,11,12]) , ...
			%					   Language([1,2,3,4],[5,6,7,8]) , Language([1,2,3,4],[9,10,11,12]) , ...
			%					   Language([5,6,7,8],[9,10,11,12]) , Language([1,2,3,4],[5,6,5,8],[9,10,11,12]) ]
			%		powerset_word_indices = { [1] , [2] , [3] , [1,2] , [1,3] , [2,3] , [1,2,3] }
			%	Note that powerset_word_indices is very compact, but loses the information of WHAT the
			%	words in the language L are. Use it wisely.

			%%%%%%%%%%%%%%%
			%% Constants %%
			%%%%%%%%%%%%%%%

			%%%%%%%%%%%%%%%
			%% Algorithm %%
			%%%%%%%%%%%%%%%

			powerset_L = [];
			powerset_word_indices = {};

			for comb_length = 1:length(obj.words)
				temp_combs = nchoosek([1:length(obj.words)],comb_length);
				for comb_ind = 1:size(temp_combs,1)
					%Create temporary language
					temp_lang = Language();
					temp_lang.words = { obj.words{temp_combs(comb_ind,:)} };
					% Update the powerset
					powerset_L = [powerset_L, temp_lang ];
					powerset_word_indices{end+1} = temp_combs(comb_ind,:);
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

		function [tf] = eq( obj, comparison_obj )
			%Description:
			%	Defines the == operator.
			%
			%Usage:
			%	tf = L1.eq(L2)
			%	tf = ( L1 == L2 )

			%% Input Processing %%

			sizeOfObj = size(obj);
			sizeOfComparison = size(comparison_obj);

			if ~all(all( sizeOfObj == sizeOfComparison ))
				error('The dimensions of the two objects is not the same.')
			end

			%% Algorithm

			if isscalar(obj)
				tf = obj.is_eq(comparison_obj);
			elseif isvector(obj)
				tf = true;
				for obj_index = 1:length(obj)
					tf = tf && obj(obj_index).is_eq(comparison_obj(obj_index));
				end
			elseif ismatrix(obj)
				tf = true;
				for obj_index1 = 1:sizeOfObj(1)
					for obj_index2 = 1:sizeOfObj(2)
						obj_at_12 = obj(obj_index1,obj_index2);
						comparison_at_12 = obj(obj_index1,obj_index2);
						tf = tf && obj_at_12.is_eq( comparison_at_12 );
					end
				end
			else
				error('The equality function only supports scalar (1x1 Language objects), vector or matrix (MxN Language object matrices) comparisons.')
			end

		end

		function tf = ne( obj, comparison_obj )
			%Description:
			%	Defines the ~= operator.
			%
			%Usage:
			%	tf = L1.ne(L2)
			%	tf = ( L1 ~= L2 )

			tf = ~( obj == comparison_obj );

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

		function symb_arr = get_all_symbols_at_idx(obj,idx)
			%Description:
			%	For each word, w, in the language 'obj', this function retrieves the symbol at idx 'idx', w[idx]
			%	and returns the list of all such symbols.
			%

			%% Algorithm %%

			symb_arr = [];
			for word_idx = 1:length(obj.words)
				symb_arr = [symb_arr,obj.words{word_idx}];
			end
			symb_arr = unique(symb_arr);

		end

		function char_arr = langtostr( obj )
			%Description:
			%	Converts a langauge of numeric arrays into a set of strings.

			%% Constants %%

			n_words = length(obj.words);

			%% Algorithm %%

			char_arr = '{ ';
			for word_idx = 1:n_words
				char_arr = [char_arr,num2str(obj.words{word_idx})];
				if word_idx ~= n_words
					char_arr = [char_arr,' , '];
				end
			end
			char_arr = [char_arr,' }'];
		end

		function char_arr = lang2str( obj )
			char_arr = obj.langtostr();
		end

		function lang_out = set_minus( obj , word_in )
			%Description: 
			%	Computes the set minues of the current language (obj) by the input word.

			% Input Processing

			if iscell(word_in)
				error('set_minus is not yet defined for inputs of type ''cell''.')
			end

			% Algorithm
			word_collection = {};

			for word_idx = 1:obj.cardinality()
				temp_word = obj.words{word_idx};
				if all(temp_word == word_in)
					continue;
				else
					word_collection{length(word_collection)+1} = temp_word;
				end
			end

			lang_out = Language(word_collection);

		end

		function [ seq_matrix , LK ] = create_belief_sequences_of_length(obj,TimeHorizon)
			%create_belief_sequences_of_length
			%Description:
			%	Defines all possible belief sequences that are based on the language obj
			%	To be a valid belief sequence, the sequence must satisfy
			%	- the Language at time t, must be a subset or equal to the belief at time t-1
			%
			%Usage:
			%	seq_matrix = L.create_belief_sequences_of_length(T)
			%
			%Outputs:
			%	seq_matrix - A matrix of "belief sequences".
			%				 This is a matrix of Language objects, where each column represents a possible, but
			%				 unchecked belief sequence. There should always be TimeHorizon number of rows.
			%	LK - A cell array containing each of the "seq_matrix" objects created on the way to reaching the
			%		 desired length of TimeHorizon

			%% Constants

			%% Algorithm 

			LK = {};
			t0 = 0;
			LK{t0+1} = obj;
			for t = 1:TimeHorizon-1
				%Initialize next level of LK
				LK{t+1} = [];
				% Collect all languages available on the previous level.
				for belief_index = 1:size(LK{t},2)
					temp_knowl_seq = LK{t}(:,belief_index);
					temp_last_lang = temp_knowl_seq(end);

					last_lang_powerset = temp_last_lang.powerset();
					for powerset_idx = 1:length(last_lang_powerset)
						% Add one of the elements of the powerset to the end of the knowledge sequence.
						LK{t+1} = [ LK{t+1} , [ temp_knowl_seq ; last_lang_powerset(powerset_idx) ] ];
					end

				end

			end

			seq_matrix = LK{end};

		end

		function [ L_Trim ] = trim_by( obj , amount_to_remove )
			%Description:
			%	Trims each word of the language by the amount amount_to_remove
			%	and returns the trimmed language.
			%	Returns an error if amount_to_remove is too large.

			%% Input Checking

			L_in = obj;

			if amount_to_remove == 0
				L_Trim = L_in;
				return;
			end

			if amount_to_remove < 0
				error(['The amount_to_remove input to trim_by() is negative (' num2str(amount_to_remove) '). Please make sure that it is always positive!' ])
			end

			if L_in.cardinality() == 0
				L_Trim = L_in;
				warning('The input language to trim_by() contains no words.')
				return;
			end

			%% Algorithm

			L_Trim = Language();
			for word_index = 1:L_in.cardinality()

				temp_word = L_in.words{word_index};

				if length(temp_word) <= amount_to_remove
					error(['The amount to remove (' num2str(amount_to_remove) ') is too large for L_in.words{' num2str(word_index) '} which is of length ' num2str(length(L_in.words{word_index})) '.' ])
				end

				L_Trim = L_Trim.union( Language( L_in.words{word_index}([1:end-amount_to_remove]) ) );

			end

		end

		function subset_flag = le(obj,L_in)
			%Description:
			%	
			%Usage:
			%	subset_flag = le(L1,L2)
			%	subset_flag = L1 <= L2;
			%

			subset_flag = obj.subseteq(L_in);
		end

		function subset_flag = ge(obj,L_in)
			%Description:
			%
			%Usage:
			%	subset_flag = ge(L1,L2)
			%	subset_flag = L1 >= L2;
			%

			subset_flag = L_in.subseteq(obj);
		end

		function [ star_sequence , star_index ] = find_minimal_covering_in( sequence0 , candidate_sequences )
			%Description:
			%	Determines which of the knowledge sequences (columns) in candidate_sequences covers
			%	sequence0.

			% Constants
			seq0_len = length(sequence0);
			[cand_len,num_candidates] = size(candidate_sequences);

			% Input Processing
			if seq0_len ~= cand_len
				error(['The length of the input sequence is ' num2str(seq0_len) ' and it cannot be covered by sequences of length ' num2str(cand_len) '.' ])
			end

			% Algorithm
			covering_sequences = []; covering_sequence_indices = [];

			for sequence_index = 1:num_candidates
				sequence_i = candidate_sequences(:,sequence_index);

				if sequence_i == sequence0
					continue; %skip the exact same one.
				end

				disp(sequence_i >= sequence0)

				if sequence_i >= sequence0
					covering_sequences = [ covering_sequences , sequence_i ];
					covering_sequence_indices = [ covering_sequence_indices , sequence_index ];
				end
			end

			star_sequence = covering_sequences;
			star_index = covering_sequence_indices;
			return

			%Find the element in the sequence which has the smallest prefix
			num_covering = size(covering_sequences,2);
			switch num_covering
				case 0
					star_sequence = [];
					star_index = -1;
				case 1
					star_sequence = covering_sequences;
					star_index = covering_sequence_indices;
				otherwise
					%Search for the minimal covering sequence from covering_sequences
					for prefix_length = [ 1 : seq0_len - 1 ]
						for covering_index = 1:num_covering
							covering_i = covering_sequences(:,covering_index);

							%When you find the first elment such that
							% sequence_i(prefix_length) ~= sequence0(prefix_length)
							%then return
							if covering_i(prefix_length) ~= sequence0(prefix_length)
								star_sequence = covering_i;
								star_index = covering_sequence_indices(covering_index);
								return;
							end
						end
					end

			end

		end


	end

end

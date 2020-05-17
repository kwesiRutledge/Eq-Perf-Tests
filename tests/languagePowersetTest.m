%languagePowersetTest.m
%Description:
%	This test is meant to test the Language.powerset() function.

%% Powerset Test #1
%  Test the number of choices

L1 = Language([1],[2],[3],[4],[5]);
powerset_cardinality = 0;
for k = 1:L1.cardinality()
	powerset_cardinality = powerset_cardinality + nchoosek(L1.cardinality(),k);
end

L1_powerset = L1.powerset();

assert(powerset_cardinality == length(L1_powerset))

%% Powerset Test #2
%	Test the content of powerset

L2 = Language([1],[2],[3],[4]);
L2_powerset = L2.powerset();

for powerset_idx = 1:length(L2_powerset)
	%Check that this contains all of the items that I want
	temp_lang = L2_powerset(powerset_idx);

	assert( temp_lang == Language([1]) || ...
			temp_lang == Language([2]) || ...
			temp_lang == Language([3]) || ...
			temp_lang == Language([4]) || ...
			temp_lang == Language([1],[2]) || ...
			temp_lang == Language([1],[3]) || ...
			temp_lang == Language([1],[4]) || ...
			temp_lang == Language([2],[3]) || ...
			temp_lang == Language([2],[4]) || ...
			temp_lang == Language([3],[4]) || ...
			temp_lang == Language([1],[2],[3]) || ...
			temp_lang == Language([1],[2],[4]) || ...
			temp_lang == Language([1],[3],[4]) || ...
			temp_lang == Language([2],[3],[4]) || ...
			temp_lang == L2 )
end

%% Powerset Test #3
%	Test the indices returned by of powerset

L3 = Language([5],[6],[7],[8]);
[~ , L3_powerset_idcs] = L3.powerset();

for powerset_idx = 1:length(L3_powerset_idcs)
	%Check that this contains all of the items that I want
	temp_idcs = L3_powerset_idcs{powerset_idx};

	contains_valid_single_digit = false;
	contains_valid_two_digits = false;
	contains_valid_three_digits = false;
	contains_valid_four_digits = false;
	switch length(temp_idcs)
		case 1
			contains_valid_single_digit =	(temp_idcs == [1]) || (temp_idcs == [2]) || ...
											(temp_idcs == [3]) || (temp_idcs == [4]);
		case 2
			contains_valid_two_digits = all(temp_idcs == [1,2]) || all(temp_idcs == [1,3]) || ...
										all(temp_idcs == [1,4]) || all(temp_idcs == [2,3]) || ...
										all(temp_idcs == [2,4]) || all(temp_idcs == [3,4]);
		case 3
			contains_valid_three_digits = 	all(temp_idcs == [1,2,3]) || all(temp_idcs == [1,2,4]) || ...
											all(temp_idcs == [1,3,4]) || all(temp_idcs == [2,3,4]);
		case 4
			contains_valid_four_digits = all(temp_idcs == [1:4]);
		otherwise
			disp('You should not reach here.')
			assert(false)
	end
	

	assert( contains_valid_single_digit || contains_valid_two_digits || contains_valid_three_digits || ...
			contains_valid_four_digits );
end
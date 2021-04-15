function tests = languageTest
	%disp(localfunctions)
	tests = functiontests(localfunctions);

function test_eq1(testCase)
	%languageTest_eq1.m
	%Description:
	%	Tests the equality operator on the Language object
	%	when two of the same language are compared.

	l1 = Language([1,2,3]);
	l2 = Language([1,2,3]);

	assert(l1 == l2)

function test_eq2(testCase)
	%test_eq2.m
	%Description:
	%	Tests the equality operator on the Language object
	%	when two different languages are compared.

	l1 = Language([1,2,3]);
	l2 = Language([1,2,2]);

	assert(~(l1 == l2))

function test_eq3(testCase)
	%test_eq3
	%Description:
	%	Tests the equality operator on the Language object
	%	when two Language arrays are compared which are not of the same size.

	l1 = Language([1,2,3]);
	l2 = Language([1,2,2]);

	seq1 = [ l1 ; l1 ; l2 ];
	seq2 = [ l1 ; l1 ];

	try
		seq1 == seq2
		assert(false)
	catch e
		assert(strcmp(e.message,'The dimensions of the two objects is not the same.'))
	end

function test_eq4(testCase)
	%test_eq4
	%Description:
	%	Tests the equality operator on the Language object
	%	when two Language matrices are compared which are not of the same size.

	l1 = Language([1,2,3]);
	l2 = Language([1,2,2]);

	mat1 = [ l1 , l2 ; l1 , l2 ; l1, l1 ];
	mat2 = [ l1 , l2 ; l1 , l2 ];

	try
		mat1 == mat2
		assert(false)
	catch e
		assert(strcmp(e.message,'The dimensions of the two objects is not the same.'))
	end

function test_eq5(testCase)
	%test_eq5
	%Description:
	%	Tests the equality operator on the Language object
	%	when two Language matrices are compared which are the same.

	l1 = Language([1,2,3]);
	l2 = Language([1,2,2]);

	mat1 = [ l1 , l2 ; l1 , l2 ; l1, l1 ];
	mat2 = [ l1 , l2 ; l1 , l2 ; l1, l1 ];

	assert(mat1 == mat2)

function test_eq6(testCase)
	%test_eq6
	%Description:
	%	Tests the equality operator on the Language object
	%	when two Language matrices are compared which are not the same.

	l1 = Language([1,2,3]);
	l2 = Language([1,2,2]);

	mat1 = [ l1 , l2 ; l1 , l2 ; l1, l1 ];
	mat2 = [ l1 , l2 ; l1 , l2 ; l1, l2 ];

	assert(mat1 == mat2)

function test_eq7(testCase)
	%test_eq7
	%Description:
	%	Tests the equality operator on the Language object
	%	when two Language vectors are compared which are not the same.

	l1 = Language([1,2,3]);
	l2 = Language([1,2,2]);

	vec1 = [ l1 ; l2 ; l1 ];
	vec2 = [ l1 ; l2 ; l2 ];

	assert(~(vec1 == vec2))

function test_eq8(testCase)
	%test_eq8
	%Description:
	%	Tests the equality operator on the Language object
	%	when two Language vectors are compared which are the same.

	l1 = Language([1,2,3]);
	l2 = Language([1,2,2]);

	vec1 = [ l1 ; l2 ; l1 ];
	vec2 = [ l1 ; l2 ; l1 ];

	assert(vec1 == vec2)

function test_contains_prefix1(testCase)
	%test_contains_prefix1
	%DescriptionL
	%	Tests the contains_prefix() function for a simple prefix contained
	%	by a given language.

	% Constants
	l1 = Language([1,2,3]);

	assert(l1.contains_prefix([1,2]))


function test_contains_prefix2(testCase)
	%test_contains_prefix2
	%DescriptionL
	%	Tests the contains_prefix() function for a simple prefix contained
	%	by a given language with multiple words.

	% Constants
	l1 = Language([3,4,4],[1,2,3],[3,4,5]);

	assert(l1.contains_prefix([1,2]))

function test_contains_prefix3(testCase)
	%test_contains_prefix3
	%DescriptionL
	%	Tests the contains_prefix() function for a simple prefix contained
	%	by a given language with multiple words where prefix is contained
	%	in two words in the language.

	% Constants
	l1 = Language([3,4,4],[1,2,3],[3,4,5]);

	[ contains_flag , contains_indicies ] = l1.contains_prefix([3,4]);

	assert(all( contains_indicies == [1,3] ) )

function test_contains_prefix4(testCase)
	%test_contains_prefix4
	%DescriptionL
	%	Tests the contains_prefix() function for a simple prefix that is
	%	NOT contained in a given language.
	
	% Constants
	l1 = Language([3,4,4],[1,2,3],[3,4,5]);

	assert( ~l1.contains_prefix([3,1]) )

function test_create_belief_sequences_of_length1(testCase)
	%test_create_belief_sequences_of_length1
	%Description:
	%	Tests the create_unchecked_belief_sequences_of_length() function for a simple
	%	Language which contains a single word.
	
	% Constants
	l1 = Language([3,4,4]);

	% Algorithm
	b1 = l1.create_belief_sequences_of_length(4);
	num_sequences = size(b1,2);

	assert( (num_sequences == 1) && all( b1 == [l1 ; l1 ; l1 ; l1] ) )

function test_create_unchecked_belief_sequences_of_length2(testCase)
	%test_create_unchecked_belief_sequences_of_length2
	%Description:
	%	Tests the create_unchecked_belief_sequences_of_length() function for a simple
	%	Language which contains a two words.
	
	% Constants
	l1 = Language([3,4,4],[1,2,1]);

	% Algorithm
	b1 = l1.create_belief_sequences_of_length(2);

	%b1(2,2)
	num_sequences = size(b1,2);

	assert( (num_sequences == 3) && ...
			( all( b1(:,1) == [ l1 ; Language([3,4,4]) ] ) ) && ...
			( all( b1(:,2) == [ l1 ; Language([1,2,1]) ] ) ) && ...
			( all( b1(:,3) == [ l1 ; l1 ] ) ) )

function test_get_feasible_combinations_of_beliefs1(testCase)
	%test_get_feasible_combinations_of_beliefs1
	%Description:
	%	Tests the get_feasible_combinations_of_beliefs() function error handling

	% Constants
	l1 = Language([3,4,4],[1,2,1]);
	oneD_sys = get_1d_lcsas();

	% Algorithm
	b1 = l1.create_belief_sequences_of_length(2);
	oneD_sys.X0 = [];

	try
		[ a , b , c ] = oneD_sys.get_feasible_combinations_of_beliefs( b1 );
	catch e 
		disp(e.message)
		assert(strcmp(e.message,['get_feasible_combinations_of_beliefs() requires nonempty choice of X0 in lcsas object.']))
	end

function test_get_feasible_combinations_of_beliefs2(testCase)
	%test_get_feasible_combinations_of_beliefs2
	%Description:
	%	Tests the get_feasible_combinations_of_beliefs() function error handling
	%	Input X0 is not a polyhedron.

	% Constants
	l1 = Language([3,4,4],[1,2,1]);
	oneD_sys = get_1d_lcsas();

	% Algorithm
	b1 = l1.create_belief_sequences_of_length(2);

	oneD_sys.X0 = 2;

	try
		[ a , b , c ] = oneD_sys.get_feasible_combinations_of_beliefs( b1 );
	catch e 
		disp(e.message)
		assert(strcmp(e.message,'get_feasible_combinations_of_beliefs() requires X0 to be a Polyhedron.'))
	end

function test_get_feasible_combinations_of_beliefs3(testCase)
	%test_get_feasible_combinations_of_beliefs3
	%Description:
	%	Tests the get_feasible_combinations_of_beliefs() function error handling
	%	Input U is empty

	% Constants
	l1 = Language([3,4,4],[1,2,1]);
	oneD_sys = get_1d_lcsas();

	% Algorithm
	b1 = l1.create_belief_sequences_of_length(2);

	try
		[ a , b , c ] = oneD_sys.get_feasible_combinations_of_beliefs( b1 );
	catch e 
		disp(e.message)
		assert(strcmp(e.message,'get_feasible_combinations_of_beliefs() requires nonempty choice of U in lcsas object.'))
	end

function test_get_feasible_combinations_of_beliefs4(testCase)
	%test_get_feasible_combinations_of_beliefs3
	%Description:
	%	Tests the get_feasible_combinations_of_beliefs() function on the test LCSAS
	%	from get_simple

	% Constants
	[ oppRotLCSAS , TimeHorizon , Pu , Pw , x0 , Px0 , X_Target ] = get_opposing_rotation_lcsas();
	l1 = oppRotLCSAS.L;

	% Algorithm
	b1 = l1.create_belief_sequences_of_length(2);

	[ possible_LK_realizations , possible_choices , choices_as_binary_flags ] = oppRotLCSAS.get_feasible_combinations_of_beliefs( b1 );

	disp(possible_LK_realizations{1})
	disp(possible_LK_realizations{1}(end).words{1})
	disp(length(possible_LK_realizations))

	assert(length(possible_LK_realizations) == 1)

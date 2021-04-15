function tests = test_knowledge_synth1()
	% disp(localfunctions)
	tests = functiontests(localfunctions);
end

function include_relevant_libraries()
	%Description:
	%	Attempts to add the relevant libraries/toolboxes to the path.

	%% Constants

	is_on_personal_mac = strcmp(getenv('USER'),'kwesirutledge');
	is_on_great_lakes = (strcmp(getenv('USER'),'krutledg') && isunix) ;

	%% Algorithm

	%Include Yalmip
	cd ..
	include_fcns2({'mosek','gurobi','MPT3'})
	cd lcsas

	%Add Local Functions to Path
	addpath(genpath('../../functions'));

	%Add System Repository to path
	system_repo_loc = '~/Documents/Michigan/Research/systemRepository';
	addpath(genpath(system_repo_loc))
end

function test_no_disturbance_synth1(testCase)
	%Description:
	%	Creates a system with no disturbance (Linear Inverted Pendulum)

	% Include Libraries
	lipm1 = FixedHeightPlanarLIPM();
	A1 = lipm1.A();
	B1 = lipm1.B();

	lipm2 = lipm1;
	lipm2.zbar_cm = 1.1; %Before this was 1.
	A2 = lipm2.A();
	B2 = lipm2.B();

	%

end

function test_get_cl_ibs1(testCase)
	%Description:
	%	Test for get_closed_loop_consistent_internal_behavior_set_matrices()
	%	for language containing 1 mode.

	% Initializing Test
	include_relevant_libraries()

	% Constants

	[oneD_sys,P_u,Pw1,Pw2,eta_v,eta_x0] = get_1d_lcsas();
	L = oneD_sys.L;
	T = length(L.words{1});
	temp_last_lang = Language(L.words{2}); %Creating simple language with one word.

	[ n_x , n_u , n_y , n_w , n_v ] = oneD_sys.Dimensions();

	x0 = zeros(n_x,1);

	K = zeros( T*n_x , T*n_w ); %randi([-10, 10], T*n_x , T*n_w );
	k = randi([-10, 10], T*n_x , 1 );

	% Create MPC Matrices
	[Sw,Su,~,J,f_bar] = oneD_sys.get_mpc_matrices('word',temp_last_lang.words{1});

	[ H_PhiI , h_PhiI ] = oneD_sys.get_consistent_internal_behavior_matrices( T , temp_last_lang , P_u , oneD_sys.X0 )

	%% Algorithm

	[~,word_id] = oneD_sys.L.contains(temp_last_lang.words{1});

	nonK_prefactor = [ 	Sw ;
		zeros(n_u*T,n_w*T) ;
		eye(n_w*T) ;
		zeros(n_x,n_w*T) ];

	K_prefactor = [ Su ;
					eye(n_u*T);
					zeros(n_w*T,n_u*T);
					zeros(n_x,n_u*T) ];

	H_p = H_PhiI * ( nonK_prefactor + K_prefactor*K);

	h_independent_factor = h_PhiI - ...
		H_PhiI * [ J*x0 + Sw*f_bar ;
				zeros(n_u*T,1) ;
				zeros(n_w*T,1) ; 
				x0 ];

	h_dependent_factor = - ...
		H_PhiI * [ Su ;
				eye(n_u*T) ;
				zeros(n_w*T,n_u*T) ; 
				zeros(n_x,n_u*T) ];

	h_p = h_independent_factor + h_dependent_factor*k;

	%%  Compare with Function Result

	[ H_cl , h_cl ] = oneD_sys.get_closed_loop_consistent_internal_behavior_set_matrices( H_PhiI , h_PhiI , x0 , K , k , repmat(temp_last_lang,T,1) );

	%% Comparison

	assert(all(all( H_cl == H_p )) && all(h_cl==h_p) )

end

function test_get_cl_ibs2(testCase)
	%Description:
	%	Test for get_closed_loop_consistent_internal_behavior_set_matrices()
	%	for language containing two words.

	% Initializing Test
	include_relevant_libraries()

	% Constants

	[oneD_sys,P_u,Pw1,Pw2,eta_v,eta_x0] = get_1d_lcsas();
	L = oneD_sys.L;
	T = length(L.words{1});
	temp_last_lang = Language(L.words{2}); %Creating simple language with one word.

	[ n_x , n_u , n_y , n_w , n_v ] = oneD_sys.Dimensions();

	x0 = zeros(n_x,1);

	K = zeros( T*n_x , T*n_w ); %randi([-10, 10], T*n_x , T*n_w );
	k = randi([-10, 10], T*n_x , 1 );

	% Create MPC Matrices
	[Sw,Su,~,J,f_bar] = oneD_sys.get_mpc_matrices('All Words');

	[ H_PhiI , h_PhiI ] = oneD_sys.get_consistent_internal_behavior_matrices( T , temp_last_lang , P_u , oneD_sys.X0 )

	%% Algorithm %%

	H_p = []; h_p = [];
	tll_card = temp_last_lang.cardinality();

	for tll_index = 1:temp_last_lang.cardinality()
 		[~,word_id] = oneD_sys.L.contains(temp_last_lang.words{tll_index});

 		nonK_prefactor = ...
 			[ zeros(n_x*(T+1),(tll_index-1)*n_w*T), Sw{word_id}, zeros(n_x*(T+1),(tll_card-tll_index)*n_w*T) ;
				zeros(n_u*T,tll_card*n_w*T) ;
				eye(tll_card*n_w*T) ;
				zeros(tll_card*n_x,tll_card*n_w*T) ];

		K_prefactor = [	Su{word_id} ;
						eye(n_u*T);
						zeros(tll_card*n_w*T,n_u*T);
						zeros(tll_card*n_x,n_u*T) ] ;

		H_p = [ H_p ; H_PhiI * ( nonK_prefactor + K_prefactor * K ) ];
		
		h_independent_factor = h_PhiI - ...
			H_PhiI * [ J{word_id}*x0 + Sw{word_id}*f_bar{word_id} ;
					zeros(n_u*T,1) ;
					zeros(tll_card*n_w*T,1) ; 
					repmat(x0,tll_card,1) ];

		h_dependent_factor = - ...
			H_PhiI * [ Su{word_id} ;
					eye(n_u*T) ;
					zeros(tll_card*n_w*T,n_u*T) ; 
					zeros(tll_card*n_x,n_u*T) ];

		h_p = [ h_p ; h_independent_factor + h_dependent_factor*k ];

	end

	%%  Compare with Function Result

	[ H_cl , h_cl ] = oneD_sys.get_closed_loop_consistent_internal_behavior_set_matrices( H_PhiI , h_PhiI , x0 , K , k , repmat(temp_last_lang,T,1) );

	%% Comparison

	assert(all(all( H_cl == H_p )) && all(h_cl==h_p) )

end
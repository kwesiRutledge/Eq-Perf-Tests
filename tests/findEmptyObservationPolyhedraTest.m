% findEmptyObservationPolyhedraTest.m
%Description:
%	Testing the ability of find_empty_observation_polyhedra() to correctly identify which items of a Polyhedron
%	list/array are empty or not.

%% Test 1: Give Two Sets, Both Full

L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);

%Create a simple Language Constrainted Switching System
A1 = [0,1;0.1,-0.05];
B1 = [0;1];
C1 = [1,0];

n_x = size(A1,1);
n_u = size(B1,2);
n_y = size(C1,1);

eta_v = 0.1; eta_w = 0.2;
Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
eta_v = 0.3;
Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

eta_u = 0; eta_x0 = 0.3;
P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

f1 = eta_w*[0;1];
f2 = eta_w*[1;0];
f3 = -f1;
f4 = -f2;

aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

bg0 = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true );

n = 3;
P1 = Polyhedron( 'lb' , -ones(1,n) , 'ub' , ones(1,n) );
P2 = 2*P1;

P_list = [P1,P2];

[ empty_set_idcs , empty_set_flags ] = bg0.find_empty_observation_polyhedra( P_list );

assert( all(empty_set_flags == false(2,1)) )

%% Test 2: Give Two Sets, One Full and One Empty

L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);

%Create a simple Language Constrainted Switching System
A1 = [0,1;0.1,-0.05];
B1 = [0;1];
C1 = [1,0];

n_x = size(A1,1);
n_u = size(B1,2);
n_y = size(C1,1);

eta_v = 0.1; eta_w = 0.2;
Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
eta_v = 0.3;
Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

eta_u = 0; eta_x0 = 0.3;
P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

f1 = eta_w*[0;1];
f2 = eta_w*[1;0];
f3 = -f1;
f4 = -f2;

aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

bg0 = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true );

n = 3;
P1 = Polyhedron( 'lb' , -ones(1,n) , 'ub' , ones(1,n) );
P2 = Polyhedron('A',[1;-1],'b',[2;-3]);

P_list = [P1,P2];

[ empty_set_idcs , empty_set_flags ] = bg0.find_empty_observation_polyhedra( P_list );

assert( all(empty_set_flags == [false;true]) )

%% Test 3: Give Two Sets, Both Empty

L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);

%Create a simple Language Constrainted Switching System
A1 = [0,1;0.1,-0.05];
B1 = [0;1];
C1 = [1,0];

n_x = size(A1,1);
n_u = size(B1,2);
n_y = size(C1,1);

eta_v = 0.1; eta_w = 0.2;
Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
eta_v = 0.3;
Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

eta_u = 0; eta_x0 = 0.3;
P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

f1 = eta_w*[0;1];
f2 = eta_w*[1;0];
f3 = -f1;
f4 = -f2;

aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

bg0 = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true );

n = 3;
P1 = Polyhedron( 'A' , [1;-1] , 'b' , [1;-4] );
P2 = Polyhedron('A',[1;-1],'b',[2;-3]);

P_list = [P1,P2];

[ empty_set_idcs , empty_set_flags ] = bg0.find_empty_observation_polyhedra( P_list );

assert( all(empty_set_flags == [true;true]) )

%  =================================
%% Test 4: Sets based on Language L1

L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);

%Create a simple Language Constrainted Switching System
A1 = [0,1;0.1,-0.05];
B1 = [0;1];
C1 = [1,0];

n_x = size(A1,1);
n_u = size(B1,2);
n_y = size(C1,1);

eta_v = 0.1; eta_w = 0.2;
Pv1 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));
Pw1 = Polyhedron('lb',-eta_w*ones(1,n_x) ,'ub',eta_w*ones(1,n_x));
eta_v = 0.3;
Pv2 = Polyhedron('lb',-eta_v*ones(1,n_y) ,'ub',eta_v*ones(1,n_y));

eta_u = 0; eta_x0 = 0.3;
P_u = Polyhedron('lb',-eta_u*ones(1,n_u) ,'ub',eta_u*ones(1,n_u));
P_x0 = Polyhedron('lb',-eta_x0*ones(1,n_x),'ub',eta_x0*ones(1,n_x));

f1 = eta_w*[0;1];
f2 = eta_w*[1;0];
f3 = -f1;
f4 = -f2;

aff_dyn_list = [	Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1), ...
					Aff_Dyn(A1,B1,f1,C1,Pw1,Pv2) ];

lcsas0 = LCSAS( aff_dyn_list , L1 , 'X0' , P_x0 );

bg0 = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true );

[L1_powerset, word_idx_powerset] = L1.powerset();

Y_Set_of = [];
for word_idx = 1:L1.cardinality();
	%For each word, create the output set of the system.
	temp_sigma = L1.words{word_idx};
	first_symb = temp_sigma(1);

	%Compute Sets
	Y_Set_of = [Y_Set_of, lcsas0.Dyn(first_symb).C * P_x0 + lcsas0.Dyn(first_symb).P_v];
end

for powerset_idx = (L1.cardinality()+1):length(word_idx_powerset)
	powerset_elt = word_idx_powerset{powerset_idx};
	powerset_elt_L = L1_powerset(powerset_idx);

	%For each word, create the output set of the system.
	temp_sigma = powerset_elt_L.words{1};
	[ ~ , L_idx_of_sigma ] = L1.contains(temp_sigma);

	%Compute Sets
	Y_Set_of = [Y_Set_of, Y_Set_of(L_idx_of_sigma)];
	for sigma_idx = 2:powerset_elt_L.cardinality()
		%For each word, create the output set of the system.
		temp_sigma = powerset_elt_L.words{sigma_idx};
		[ ~ , L_idx_of_sigma ] = L1.contains(temp_sigma);

		Y_Set_of(end) = Y_Set_of(end).intersect( Y_Set_of(L_idx_of_sigma) );
	end

end	

[ empty_set_idcs , empty_set_flags ] = bg0.find_empty_observation_polyhedra( Y_Set_of );

assert( all(empty_set_flags == [false(length(Y_Set_of),1)]) )
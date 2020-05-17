% emptyBeliefGraphTest1.m
%Create Empty BG Test

L1 = Language([1,2,1,2],[3,4,3,4],[5,1,1,1]);
T = 4;

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
lcsas0.X0 = P_x0;

empty_bg = BeliefGraph( lcsas0 , P_u , P_x0 , 'return_empty' , true );

assert( isempty(empty_bg.N) & isempty(empty_bg.E) )
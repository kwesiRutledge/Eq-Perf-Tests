% reduced order observer design for ACC

clear all;
close all;

con = constants_normal; 

kappa = con.f1_bar/con.mass;
ekt = exp(-kappa*con.dT);

A = [ ekt 				0 		0 ; 
      (ekt-1)/kappa 	1 		con.dT;
      0 				0 		1]; 
B = (1/(con.f1_bar^2))*[ con.f1_bar*(1-ekt) ;
                         con.mass*(1-ekt) - con.dT*con.f1_bar ;
                         0 ];
B_cond_number = max(abs(B));
B = B/B_cond_number;

E = [0 ; con.dT^2/2; con.dT];
K = [ (con.f0_bar/con.f1_bar)*(ekt-1) ;
          (con.f0_bar/(con.f1_bar^2) )*(con.dT*con.f1_bar + con.mass*(ekt-1)) ;
          0 ];
      
% Ad = [Aaa, Aau;
%       Aua, Auu];
Aaa = A(1:2,1:2); 
Aau = A(1:2,3); %becomes C
Aua = A(3,1:2); 
Auu = A(3,3); %becomes A

Ba = B(1:2,1);
Bu = B(3,1);
Fa = K(1:2,1);
Fu = K(3,1);

Ea = E(1:2,1);
Eu = E(3,1);

% disturbance model

Hd = [1.;-1.]
hd = [con.al_max;-con.al_min]

W = Polyhedron(Hd,hd);

p = 'inf';

mu = 0.05*(con.vl_max-con.vl_min); 
% we are only estimating lead vehicle speed, allow 5% error in estimation

L = ro_observer(Auu, Aau, Eu, Ea, mu, W, p)


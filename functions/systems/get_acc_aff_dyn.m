function [dyn_disc] = get_acc_aff_dyn()
	%Constants
	m = 1370; %kg
	k0 = 7.58; %kg
	k1 = 9.9407;
	Ts = 0.5;
	eta_w = 0.1;
	eta_v = 0.05;

	kappa = k1/m;

	%Create ACC System
	acc0.A = [  exp(-kappa*Ts),             0,  0 ;
	            (exp(-kappa*Ts)-1)/kappa,   1,  Ts;
	            0,                          0,  1];
	acc0.B = (1/(k1.^2))*[ (1-exp(-kappa*Ts))*k1 ; m*(1-exp(-kappa*Ts))-k1*Ts ; 0 ];
	acc0.C = [1, 0, 0; 0, 1, 0];
	acc0.E = [0; (Ts^2)/2; Ts ];
	acc0.f = [-(k0/k1)*(1-exp(-kappa*Ts)) ; -(k0/(k1^2))*(m*(1-exp(-kappa*Ts)) - k1*Ts) ; 0 ];

	acc0.eta_w = eta_w;
	acc0.eta_v = eta_v;

	dyn_disc = Aff_Dyn(acc0.A,acc0.B,acc0.f,acc0.C,acc0.eta_w,acc0.eta_v,acc0.E,eye(size(acc0.C,1)));
end
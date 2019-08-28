function [results] = observer_comparison54( varargin )
	%Description:
	%	This script is meant to test the formation control example that Necmiye recommended for the Belief Graph Control Methods.
	%	It first tests the feasibility of recovery when the mode is perfectly observed and then is meant to show that when the mode
	%	is not observed, our previous solution cannot be directly applied.

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%
	cube_dim.x = 3;
	cube_dim.y = 2;

	deltat = 0.1;

	eta_w = 0.25;
	eta_v = 0.5;
	eta_u = 3;

	[ad1,map] = get_consensus_dyn(cube_dim.x,cube_dim.y,deltat, ...
									'disturb_params',eta_w,eta_v);

	results.params.cube_dim = cube_dim;
	results.params.dt = deltat;
	results.params.ad1 = ad1;
	results.params.des_pos = map;

	% Get the Consensus System for when a strong wind is blowing to the north
	n_w = size(ad1.B_w,2);
	unit_cube = Polyhedron('lb',-ones(1,n_w),'ub',ones(1,n_w));

	P_w1 = eta_w * unit_cube;

	P_w2 = P_w1 + [0;eta_w];

	[ad2,~] = get_consensus_dyn(cube_dim.x,cube_dim.y,deltat, ...
									'disturb_params',P_w2,eta_v*Polyhedron('lb',-ones(1,size(ad1.C_v,2)),'ub',ones(1,size(ad1.C_v,2))));

	results.params.PW1 = P_w1;
	results.params.PW2 = P_w2;
	results.params.ad2 = ad2; 

	% Create the 'all-encapsulating' system
	P_w3 = Polyhedron('lb',-eta_w*ones(1,n_w),'ub',eta_w*ones(1,n_w)+[0,eta_w]);
	[ad3,~] = get_consensus_dyn(cube_dim.x,cube_dim.y,deltat, ...
									'disturb_params',P_w3,eta_v*Polyhedron('lb',-ones(1,size(ad1.C_v,2)),'ub',ones(1,size(ad1.C_v,2))));

	results.params.PW3 = P_w3;
	results.params.ad3 = ad3;

	ad_arr = [ad1,ad2,ad3];

	% Create Input Set
	Pu = eta_u*Polyhedron('lb',-ones(1,size(ad_arr(1).B,2)),'ub',ones(1,size(ad_arr(1).B,2)));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Attempt Recovery on similar words %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	M1 = 1.5;

	len_nominal = 3;
	len_switch = 3;

	%Get All Words
	L2 = {};
	for word_idx =1:length(len_switch)
		temp_word = [ones(1,len_nominal+word_idx-1),2*ones(1,len_switch-word_idx+1)];
		L2{word_idx} = temp_word;
	end
	% word1 = [ones(1,len_nominal),2*ones(1,len_switch)];
	% word2 = [ones(1,len_nominal),1,2*ones(1,len_switch-1)];
	% word3 = [ones(1,len_nominal),1,1,2*ones(1,len_switch-2)];

	wordC = [ones(1,len_nominal),3*ones(1,len_switch-1),2]; %The "L-Star" like word

	%Synthesize for The individual words
	%L2 = {word1,word2,word3};
	[optim_out1,contr1] = ad_arr.rec_synthesis('Equalized','prefix','Minimize M2',M1,L2,'Pu',Pu,'System Type','Switched');

	[optim_out2,contr2] = ad_arr.rec_synthesis('Equalized','prefix','Minimize M2',M1,{wordC},'Pu',Pu,'System Type','Switched');

	results.exp1.L = L2;
	results.exp1.optim_out_ideal = optim_out1;
	results.exp1.contr_ideal = contr1;
	results.exp1.optim_out_conservative = optim_out2;
	results.exp1.contr_conservative = contr2;

end
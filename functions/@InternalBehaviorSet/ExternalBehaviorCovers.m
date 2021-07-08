function [tf] = ExternalBehaviorCovers(ibs1,ibs2)
	%Description:
	%
	%Usage:
	%	tf = ibs1.ExternalBehaviorCovers(ibs2)

	%% Input Checking

	%% Constants
	System = ibs1.System;
	cg = constr_gen(0);

	ops0 = sdpsettings('verbose',0);

	%% Algorithm

	Rt_1 = ibs1.SelectExternalBehavior();
	Rt_2 = ibs2.SelectExternalBehavior();

	external_beh_dim = size(Rt_1,1);

	[dual_vars, constr] = cg.create_sadraddini_AH_inclusion_constr( ...
							zeros(external_beh_dim,1) , Rt_2 , [ibs2.A;ibs2.Ae;-ibs2.Ae] , [ibs2.b;ibs2.be;-ibs2.be] , ...
							zeros(external_beh_dim,1) , Rt_1 , [ibs1.A;ibs1.Ae;-ibs1.Ae] , [ibs1.b;ibs1.be;-ibs1.be] );

	diagnostics = optimize(constr,[],ops0);

	tf = (diagnostics.problem == 0);

end
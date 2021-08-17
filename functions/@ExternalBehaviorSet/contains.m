function [tf] = contains(ebs,external_behavior)
	%Description:
	%	Determines if an external behavior is part of this target ExternalBehaviorSet.
	%
	%Usage:
	%	tf = ebs.contains( external_behavior )

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	Settings = ebs.Settings;
	System = ebs.System;
	[ n_x , n_u , n_y , n_w , n_v ] = System.Dimensions();

	ibs = ebs.ParentInternalBehaviorSet;

	ib_dim = ibs.Dim;
	eb_dim = ebs.Dim;

	verbosity = 1;
	options = sdpsettings('verbose',verbosity);

	first_hypotheses = ebs.KnowledgeSequence(1);
	final_hypotheses = ebs.KnowledgeSequence(end);
	num_final_hypotheses = final_hypotheses.cardinality();

	%%%%%%%%%%%%%%%
	%% Algorithm %%
	%%%%%%%%%%%%%%%

	switch Settings.fb_type
	case 'state'

		ibsAsPolyhedron = ibs.ToPolyhedron();
		x0 = external_behavior([1:n_x]);

		if ebs.t == 0
			tf = ibsAsPolyhedron.contains(x0);
			return
		end

		%Construct w's that would complete the 
		w_hypothesis = {};
		for hypothesis_index = 1:first_hypotheses.cardinality()
			w_hypothesis{hypothesis_index} = [];
		end

		for tau = 0:ebs.t-1
			x_tau = external_behavior([tau*n_x+1:(tau+1)*n_x]);
			x_taup1 = external_behavior([(tau+1)*n_x+1:(tau+2)*n_x]);
			u_tau = external_behavior(n_x*(ebs.t+1)+[tau*n_u+1:(tau+1)*n_u]);

			for hypothesis_index = 1:first_hypotheses.cardinality()
				hypothesis = first_hypotheses.words{hypothesis_index};
				hypothesisDynamicsAtTau = System.Dyn( hypothesis(tau+1) );

				w_tau = hypothesisDynamicsAtTau.find_w_that_completes( x_tau , u_tau , x_taup1 );
				w_hypothesis{hypothesis_index} = [ w_hypothesis{hypothesis_index} ; w_tau ];
			end

		end

		% Construct internal_behavior
		internal_behavior = external_behavior;
		for hypothesis_index = 1:first_hypotheses.cardinality()
			internal_behavior = [ internal_behavior ; w_hypothesis{hypothesis_index} ];
		end
		internal_behavior = [internal_behavior;x0];

		
		tf = ibsAsPolyhedron.contains(internal_behavior);

		% % Construct Optimization Variables
		% internal_behavior = sdpvar(ib_dim,1,'full');

		% % Construct Constraints
		% ib_match_constraint = [ [eye(eb_dim) , zeros(eb_dim,ib_dim-eb_dim)]*internal_behavior == external_behavior ];

		% valid_ib_constraint = 	[ ibs.A * internal_behavior <= ibs.b ] + ...
		% 						[ ibs.Ae * internal_behavior == ibs.be ];

		% %Run optimization
		% diagnostics = optimize(	ib_match_constraint+valid_ib_constraint, ...
		% 						[], options);

		% tf = (diagnostics.problem == 0);
	otherwise
		error(['Unexpected fb_type for ExternalBehaviorSet using contains(): ' Settings.fb_type ])
	end

end
function [results] = observer_comparison53( varargin )
%	Description:
%		This example will test 2 new classes for constructing the BeliefGRAPH (no longer a belief tree).

	%%%%%%%%%%%%%%%%%%%%%%
	%% Input Processing %%
	%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	dim = 2;

	L1 = {[1,1,1],[1,2,1],[3,2,1]};
	L2 = {[1,1,2],[1,2,2],[1,2,4]};

	%Create a simple Language Constrainted Switching System
	A1 = eye(dim);
	B1 = eye(dim);
	C1 = eye(dim);
	f1 = [0;1];

	eta_v = 0; eta_w = 0.2;
	Pv1 = Polyhedron('lb',-eta_v*ones(1,dim) ,'ub',eta_v*ones(1,dim));
	Pw1 = Polyhedron('lb',-eta_w*ones(1,dim) ,'ub',eta_w*ones(1,dim));

	ad1 = Aff_Dyn(A1,B1,f1,C1,Pw1,Pv1);

	f2 = [1;0];
	f3 = -f1;
	f4 = -f2;

	lcss1 = [	ad1,...
				Aff_Dyn(A1,B1,f2,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f3,C1,Pw1,Pv1),...
				Aff_Dyn(A1,B1,f4,C1,Pw1,Pv1)];

	%%%%%%%%%%%%%%%%%
	%% Experient 1 %%
	%%%%%%%%%%%%%%%%%

	disp('Experiment 1: Create a Dummy Belief Node')

	bn0 = BeliefNode({L1{[1:2]}},1);

	results.exp1.L = L1;
	results.exp1.bn = bn0;

	disp('Done!')
	disp(' ')

	%%%%%%%%%%%%%%%%%%
	%% Experiment 2 %%
	%%%%%%%%%%%%%%%%%%

	disp('Experiment 2: Implementing Power Set Operation for our Belief Node''s subL property')

	bn1 = BeliefNode(L2,1);

	c_level = [bn1];

	%Get All Combinations of the node's subset
	node_p_set = {};
	node_ind = 1;
	for comb_length = 1:length(c_level(node_ind).subL)
		temp_combs = nchoosek([1:length(c_level(node_ind).subL)],comb_length);
		for comb_ind = 1:size(temp_combs,1)
			node_p_set{end+1} = temp_combs(comb_ind,:);
		end
	end

	temp2 = bn1.idx_powerset_of_subL();

	results.exp2.L = L2;
	results.exp2.node_p_set = node_p_set;
	results.exp2.function_pset = temp2;

	disp('Comparing Languages.')
	disp(['node_p_set == bn1.powerset_of_subL(): ' num2str( language_eq(node_p_set,temp2) )])

	disp('Done')
	disp(' ')

	%%%%%%%%%%%%%%%%%%
	%% Experiment 3 %%
	%%%%%%%%%%%%%%%%%%

	disp('Experiment 3: Testing the Post Operator');
	
	%Define Sets
	eta_u = 0.25; eta_x0 = 0.5;
	P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
	P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));

	bn2 = BeliefNode(L2,0);
	temp_post = bn2.post(lcss1,P_u,P_x0);

	results.exp3.L = L2;
	results.exp3.bn = bn2;
	results.exp3.post = temp_post;

	disp('Possible Belief Nodes at the Next Step:')
	temp_post

	disp('Done')
	disp(' ')

	%%%%%%%%%%%%%%%%%%
	%% Experiment 4 %%
	%%%%%%%%%%%%%%%%%%

	clear temp_post

	disp('Experiment 4: Verifying that the Post Operator can create a proper array of BeliefNode objects.');
	
	%Define Sets
	% eta_u = 0.25; eta_x0 = 0.5;
	% P_u = Polyhedron('lb',-eta_u*ones(1,dim) ,'ub',eta_u*ones(1,dim));
	% P_x0 = Polyhedron('lb',-eta_x0*ones(1,dim),'ub',eta_x0*ones(1,dim));

	bn3 = BeliefNode(L2,1);
	temp_post = bn3.post(lcss1,P_u,P_x0);

	results.exp4.L = L2;
	results.exp4.bn = bn2;
	results.exp4.post = temp_post;

	%%%%%%%%%%%%%%%%%%
	%% Experiment 5 %%
	%%%%%%%%%%%%%%%%%%

	disp('Experiment 5: Testing the Belief Node Equality function')

	bn4 = BeliefNode(L2,1);

	disp(['Comparing nodes:'])
	disp(['bn4 == bn2 : ' num2str(bn4.is_eq(bn2))])
	disp(['bn4 == bn3 : ' num2str(bn4.is_eq(bn3))])

	disp('Done!')
	disp(' ')

	%%%%%%%%%%%%%%%%%%
	%% Experiment 6 %%
	%%%%%%%%%%%%%%%%%%

	disp('Experiment 6: Testing the first iteration of the BeliefGraph construction')

	% bg = BeliefGraph();
	% bg = bg.construct(lcss1,L2, P_u, P_x0)
	L3 = Language();
	L3.words = L2;

	bg = BeliefGraph(lcss1,L3, P_u, P_x0)

	disp('BeliefGraph successfully created!')
	
	results.exp6.bg = bg;

	disp('Done!')
	disp(' ')


	%%%%%%%%%%%%%%%%%%
	%% Experiment 7 %%
	%%%%%%%%%%%%%%%%%%

	disp('Experiment 7: Testing BeliefGraph visualization for debugging purposes')

	bg.plot()

	disp('Are the belief nodes the same as the graph visualizes?')
	disp(['bg.get_leaf_node_idxs() == [6,7,8,9] ? ' num2str(all(bg.get_leaf_node_idxs() == [6,7,8,9])) ])
	disp(' ')
	disp(['bg.pre(8) == [4,5]?' num2str(all(bg.pre(8) == [4,5]))])
	disp(' ')
	disp(['bg.all_words_start_with_root(L2) ? ' num2str(L3.all_words_start_with_root())])
	temp_lang = Language([2,2,3,2,2]);
	disp(['bg.all_words_start_with_root({[2,2,3,2,2]}) ? ' num2str(temp_lang.all_words_start_with_root())])
	disp(' ')
	
	temp_lang = bg.prepend_any_valid_node([8]);
	disp(['bg.prepend_any_valid_node([8]) == {[4,8],[5,8]} [' num2str(temp_lang.words{1}) '] , [' num2str(temp_lang.words{2}) ']' ]) 
	temp_lang = Language([8],[9]);
	temp_lang2 = bg.prepend_any_valid_node(temp_lang);
	disp(['bg.prepend_any_valid_node({[8],[9]}) == {[4,8],[5,8],[5,9]} [' num2str(temp_lang2.words{1}) '] , [' num2str(temp_lang2.words{2}) '], [' num2str(temp_lang2.words{3}) ']' ])
	disp(' ')

	temp_lang1 = Language([1,2],[4,5]);
	temp_lang2 = Language([4,5],[3,4],[1:3],[1,2]);
	temp_lang = temp_lang1.union(temp_lang2);
	disp(['temp_lang1.union(temp_lang2) = {' num2str(temp_lang.words{1}) ',' num2str(temp_lang.words{2}) ',' num2str(temp_lang.words{3}) ', ...'])
	disp(' ')
	[belief_lang,leaf_node_idxs] = bg.get_belief_language();
	disp(' ')
	disp('Done!')
	disp(' ')

	results.exp7.belief_lang = belief_lang;

end
function sample_out = sample_once_from( polytope_in )
	%Description:
	%	Samples a point randomly from the polytope provided 


	%% Constants %%
	num_verts = size(polytope_in.V,1);
	mu0 = 1; %Mean of 1 for all variables is important.

	%% Algorithm %%
    
    if num_verts == 1
        %If there is only one vertex given, then return it.
        %It is the only thing in this set.
        sample_out = polytope_in.V';
        return
    end

	convex_comb = exprnd(mu0,num_verts,1);
	convex_comb = convex_comb./repmat(sum(convex_comb),num_verts,1);

	sample_out = polytope_in.V'*convex_comb;

end
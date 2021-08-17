function sample_matrix_out = sample_n_times_from( polytope_in , n )
	%Description:
	%	Samples the polytope polytope_in and returns the matrix of samples

	%% Algorithm %%
	sample_matrix_out = [];

	for sample_index = 1 : n
		sample_matrix_out = [ sample_matrix_out , sample_once_from( polytope_in ) ];
	end

end
function [ initial_nodes ] = get_initial_beliefnodes( varargin )
	%Desctiption:
	%	This function uses the definition of the current belief nodes to
	%
	%Usage:
	%	[ initial_nodes ] = bg.get_initial_beliefnodes()
	%	[ initial_nodes ] = bg.get_initial_beliefnodes(P_x0)

	%% Input Processing %%
	bg = varargin{1};

	%Check to see if the initial state set is embedded in the belief graph's LCSAS member.
	lcsas_in = bg.lcsas;
	if isfield(lcsas,'X0')
		X0 = lcsas.X0;
	else
		if nargin ~= 2
			error('Two arguments are needed if LCSAS does not have an initial state set defined.')
		end
		X0 = varargin{2};
	end

	%% Algorithm %%

	Y_Set_of = {};
	for word_idx = 1:bg.ModeLanguage.cardinality()
		%For each word, create the output set of the system.
		temp_sigma = bg.ModeLanguage.words{word_idx};
		first_symb = temp_sigma(1);

		%Compute Sets
		Y_Set_of{word_idx} = lcsas_in.Dyn(first_symb).C * X0 + lcsas_in.Dyn(first_symb).P_v;
		
	end

end
function deconstruct_and_save_to(cbc, filenameIn)
	%deconstruct_and_save_to
	%Description:
	%	Save each of the elements of the ConsistentBeliefController to the directory given as filenameIn.
	%

	%% Constants
	System = cbc.System;
	KnowledgeSequences = cbc.KnowledgeSequences;

	L = System.L;

	T = length(L.words{1});

	date_converted = datestr(now,'ddmmmyyyy-HHMM');

	%% Input Processing
	if ~exist('filenameIn')
		filenameIn = ['deconstructed_controller_' date_converted '.mat'];
	end



	%% Algorithm


	% Save Dynamics Info

	lcsas_A_matrices = {};
	lcsas_B_matrices = {};
	lcsas_K_matrices = {};
	for dynamicsIndex = 1 : length(System.Dyn)
		lcsas_A_matrices{dynamicsIndex} = System.Dyn(dynamicsIndex).A;
		lcsas_B_matrices{dynamicsIndex} = System.Dyn(dynamicsIndex).B;
		lcsas_K_matrices{dynamicsIndex} = System.Dyn(dynamicsIndex).f;

		lcsas_W_A_matrices{dynamicsIndex}  = System.Dyn(dynamicsIndex).P_w.A;
		lcsas_W_b_matrices{dynamicsIndex}  = System.Dyn(dynamicsIndex).P_w.b;
		lcsas_W_Ae_matrices{dynamicsIndex} = System.Dyn(dynamicsIndex).P_w.Ae;
		lcsas_W_be_matrices{dynamicsIndex} = System.Dyn(dynamicsIndex).P_w.be;

	end

	x0 = System.X0.V';

	% Save Input Constraints for System
	U_A = System.U.A;
	U_b = System.U.b;

	% Save Language Data for System
	LAsCellArray = L.words;

	% Save Knowledge Sequence Info in a matrix containing matrices.
	for knowl_seq_index = 1:size(KnowledgeSequences,2)
		for tauPOne = 1:T
			% Save the information about indices
			tempLanguage = KnowledgeSequences(tauPOne,knowl_seq_index);
			tempLanguageAsMatrix = [];
			for wordIndex = 1:tempLanguage.cardinality()
				tempLanguageAsMatrix = [ tempLanguageAsMatrix ; tempLanguage.words{wordIndex} ];
			end
			knowl_matrix{tauPOne,knowl_seq_index} = tempLanguageAsMatrix;
		end
	end

	% Save Internal Behavior Set Objects
	for knowl_seq_index = 1:size(KnowledgeSequences,2)
		for tauPOne = 1:T
			% Save internal behavior set objects.
			ibs_matrix_A{tauPOne,knowl_seq_index}  = cbc.ConsistencySets{tauPOne,knowl_seq_index}.ParentInternalBehaviorSet.A;
			ibs_matrix_b{tauPOne,knowl_seq_index}  = cbc.ConsistencySets{tauPOne,knowl_seq_index}.ParentInternalBehaviorSet.b;

			ibs_matrix_Ae{tauPOne,knowl_seq_index} = cbc.ConsistencySets{tauPOne,knowl_seq_index}.ParentInternalBehaviorSet.Ae;
			ibs_matrix_be{tauPOne,knowl_seq_index} = cbc.ConsistencySets{tauPOne,knowl_seq_index}.ParentInternalBehaviorSet.be;


			% knowl_matrix{tauPOne,knowl_seq_index} = tempLanguageAsMatrix;
		end
	end


	% Save Gain Information
	controller_K_matrices = cbc.K_set;
	controller_k_matrices = cbc.k_set;

	% Save the file
	save(filenameIn, ...
		'date_converted', ...
		'lcsas_A_matrices','lcsas_B_matrices','lcsas_K_matrices', ...
		'lcsas_W_A_matrices','lcsas_W_b_matrices','lcsas_W_Ae_matrices','lcsas_W_be_matrices', ...
		'x0', 'U_A' , 'U_b' , 'LAsCellArray', ...
		'knowl_matrix', ...
		'ibs_matrix_A', 'ibs_matrix_b', 'ibs_matrix_Ae', 'ibs_matrix_be', ...
		'controller_K_matrices', 'controller_k_matrices');

end
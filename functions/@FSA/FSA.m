classdef FSA
	%Description:
	%	The model of a finite-state automaton according to the Wang et. al. paper from 2007.
	
	%Member Variables
	%	X 		= Finite State Space of the Automaton
	%			Data Structure: |X| x 1 Vector
	%
	%	E 		= Finite Event Space of the Automaton
	%			  Data Structure: |E| x 1 Vector
	%	
	%	f 		= Partial Transition Function f(x,e) = y or (x,e,y) \in f
	%		  	  Data Structure: |f| x 3 Matrix
	%
	%	Gamma 	= Active Event Function
	%			  Describes whether state-action pair is feasible by whether or not the pair exists
	%				in the following relation. (pair would be a row in the matrix)
	%			  Data Structure: |G| x 2 Vector
	%
	%	x0 		= Initial State
	%			  Data Structure: 1 x 1 Vector
	%

	properties
		X;
		E;
		f;
		Gamma;
		x0;
		ind_func;
	end

	methods
		%Constructor
		function fsa = FSA( X , E , f , Gamma , x0 )
			fsa.X = X;
			fsa.E = E;
			fsa.f = f;
			fsa.Gamma = Gamma;
			fsa.x0 = x0;
		end
		%Set of transitions TR(G)
		function [ sot ] = TR(obj)
			sot = obj.f(:,[1 2]);
		end 
		%Index Function Definition
		% function [ obj ] = def_script_I(obj,script_I)
		% 	obj.ind_func = script_I;
		% end
		function [ I_val ] = eval_I(obj,x0,e)
			I_val = obj.ind_func( findrows(obj.ind_func(:,[1 2]),[x0 e]) , 3 );
		end
		% Information Mapping, or Projection function
		function [ theta_s ] = info_map(obj,s)
			%Description:
			%	Recursion is written in the paper.

			if isempty(s)
				theta_s = [];
			else
				if s(end) == 1
					theta_s = [ obj.info_map(s([1:end-1])), s(end) ];
				else
					theta_s = [ obj.info_map(s([1:end-1])) ];
				end
			end
		end
		%%%%%%%%%%%%%%%%%%%%%%
		%% HELPER FUNCTIONS %%
		%%%%%%%%%%%%%%%%%%%%%%
		function [cp] = cart_prod(obj, set1 , set2)
			%Description:
			%	For every state in set1 (a column vector) and set2, create a pair (row in the result cp).
			cp = [];
			for item1 = set1'
				for item2 = set2'
					cp = [cp; item1,item2];
				end
			end
		end
	end


end
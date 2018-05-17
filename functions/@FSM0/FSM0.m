classdef FSM0
	%Description:
	%	This class is a Finite State Machine as defined by De Santis et. al. in a relevant 2017 Automatica paper.
	%	We use their FSM for our "automata" because they've developed algorithms that can discuss the executions
	%	of such objects in a flexible and intuitive way.

	%Member Variables:
	%	X - 	A set of potential states that the fsm can be in.
	%			Data Type: |X| x 1, Array of integers.
	%
	%	X0 - 	The set of potential initial states of the fsm.
	%			Data Type: |X0| x 1, Array of integers.
	%
	%	Y - 	The set of possible outputs.
	%			Data Type: |Y| x 1, Array of integers.
	%
	%	H -		Relation between given states and output values.
	%			Data Type: |X| x 2 Matrix. First entry of each row is value from X, second entry is output from Y
	%
	%	Delta - Transition Relation Between States
	%			Data Type: 2 x |Delta| matrix. First entry is originating state (from X), second entry is target state (from X).
	%

	properties
		X;
		X0;
		Y;
		H;
		Delta;
		Pi;
		Theta;
	end

	methods
		%Constructor
		function fsm = FSM0( X , X0 , Y , H , Delta )
			fsm.X = X;
			fsm.X0 = X0;
			fsm.Y = Y;
			fsm.H = H;
			fsm.Delta = Delta;

			fsm.Theta = fsm.create_Theta;
			fsm.Pi = fsm.create_Pi;

		end
		function y_out = H_of(obj,x)
			%Description:
			%	For state x, produces the value y that corresponds to it according to H.
			y_out = obj.H( find(obj.H(:,1)==x) , 2 );
		end 
		function [fsm_Pi] = create_Pi( obj )
			%Description:
			%	Create Pi the set of state-pairs for which, in each pair, both states share the same output value.
			fsm_Pi = [];
			for x1 = obj.X'
				for x2 = obj.X'
					%Add to Pi only if the outputs are identical.
					if obj.H_of(x1) == obj.H_of(x2)
						fsm_Pi = [fsm_Pi; [x1,x2] ];
					end
				end
			end
			obj.Pi = fsm_Pi;
		end
		function [fsm_Theta] = create_Theta( obj )
			%Description:
			%	Create Theta, the set of state-pairs for which, in each pair, the states are the same.
			fsm_Theta = [];
			for x1 = obj.X' %'
				fsm_Theta = [fsm_Theta; x1,x1 ];
			end
			obj.Theta = fsm_Theta;
		end

		function [s_states] = succ( obj , x0 )
			%Description:
			%	Find the set of states that can possibly follow the current state x0.
			%
			%Assumption:
			%	Assumption is that x0 is a SINGLE state.
			%
			s_states = obj.Delta( find(obj.Delta(:,1)==x0) , 2 );
		end

		function [p_states] = pre(obj, x0)
			%Description:
			%	Find the set of states that can possibly precede the current state x0.
			%
			%Assumption:
			%	Assumption is that x0 is a SINGLE state.
			%
			p_states = obj.Delta( find(obj.Delta(:,2)==x0) , 1 );
		end

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

		function [B_kp1] = B_recur(obj,B_k)
			%Description:
			%	Attempting to implement the recursion described by Definition 10 and Lemma 11 in De Santis et. al.
			%	Recursion:
			%		0 - Select the target set Sigma, a subset of Pi.
			%		1 - B1(Sigma) = Sigma
			%		2 - B_{k+1}(Sigma) = { (i,j) in B_k(Sigma) | (pre(i) x pre(j)) \cap B_k(Sigma) \neq \emptyset }
			
			B_kp1 = [];
			for ind = 1:size(B_k,1)
				%Iterate through every pair in B_k
				% 1. Compute the cartesian products of the two entries in tuple ij
				temp_cp = obj.cart_prod(obj.pre( B_k(ind,1) ), obj.pre( B_k(ind,2)) );
				% 2. Intersect the cartesian product with B_k(\Sigma)
				temp_isect = intersect(temp_cp,B_k,'rows');
				% 3. If intersection is nonempty, then add the appropriate row to B_{k+1}(Sigma)
				if ~isempty(temp_isect)
					B_kp1 = [B_kp1; B_k(ind,:)];
				end
			end
		end

	end

end
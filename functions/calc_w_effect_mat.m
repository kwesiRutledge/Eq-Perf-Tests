function [G] = calc_w_effect_mat(varargin)
	%calc_w_effect_mat
	%	Description:
	%		Creates the mpc matrix associated with the effects of disturbances
	%		over the course of an important time horizon.
	%		It has been modified to also handle the switched system in our most recent work.
	%	
	%	Usage:
	%		H = calc_w_effect_mat(A,T)
	%		H = calc_w_effect_mat(sys_arr,sigma)
	%

	%Depending on whether or not the first input is an 'Aff_Dyn' object, do different things.

	if isa(varargin{1},'Aff_Dyn')
		%If the first input is an array of affine dynamics objects,
		%then let's do the complicated version of this script.

		%Input Processing
		sys_arr = varargin{1};
		sig = varargin{2};

		%Constants
		n = size(sys_arr(1).A,1);
		T = length(sig);

		G = zeros(n*(T+1),n*T);

		nonzero_part = []; %Nonzero part of the row.

		for i = 1:T
			if i == 1 
				nonzero_part = [ eye(n) ];
            else
                %Note the integer indexing the word sigma starts from 2 as
                %defined in Skaf et. al.'s original work.
				nonzero_part = [  sys_arr(sig(i)).A*nonzero_part, eye(n) ];
			end

			%Update G Matrix
			G([i*n+1:(i+1)*n],[1:n*i]) = nonzero_part;

		end

	else
		%If the first input is a numeric matrix for affine dynamics,
		%then let's do the simple version of this function.

		%Input Processing
		A = varargin{1};
		T = varargin{2};

		%Constants
		G = zeros(size(A,1)*(T+1),size(A,2)*T);

		for i = 1 : T

			temp = [];
			for k = 1:i
				temp = [temp A^(i-k)];
			end
			
			G([i*size(A,1)+1:(i+1)*size(A,1)],[1:size(A,2)*i]) = temp;
		
		end
	end

end

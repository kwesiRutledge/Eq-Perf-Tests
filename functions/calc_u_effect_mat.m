function [H] = calc_u_effect_mat(A,B,T)

	%Constants
	H = zeros(size(A,1)*(T+1),size(B,2)*T);

	for i = 1: T

		temp = [];
		for k = 1:i
			temp = [temp A^(i-k)*B];
		end
		
		H([i*size(A,1)+1:(i+1)*size(A,1)],[1:size(B,2)*i]) = temp;
	
	end

end
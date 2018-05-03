classdef FHAE_sb
%Description:
%	This class is the state based controller that was developed in my most recent write up. (#31)
	properties
		L;
		F_set;
		u0_set;
	end

	methods
		%Constructor
		function contr = FHAE_sb( L , F_set , u0_set )
			%Simply copy everything over
			contr.L = L;
			contr.F_set = F_set;
			contr.u0_set = u0_set;
		end

		%Apply the correct control
		function u = apply_control( obj , observed_w , y_vec )
			disp('This function doesn''t work yet!')

			ow_len = length(observed_w);

			%Match observed_w with a word from L
			matching_mat = obj.L(:,1:ow_len);
			matching_mat = repmat(observed_w,size(obj.L,1),1) == matching_mat;

			matching_locs = all(matching_mat')';

			first_match_ind = find(matching_locs,1)

			% if ( size(y_vec,1) ~= F_set() )

			u = 2;
		end
	end

end

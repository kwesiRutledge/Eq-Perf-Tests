function [ opt_out, fb ] = eq_rec_design_t( varargin )
%Description:
%	Searches for a feasible, time-based feedback that satisfies the parameters given in the
%	Equalized Recovery Problem:
%		(M1,M2,T)
%
%	There are 2 different problems that we consider:
%	'Min_M2' , 'Feasible Set'. Which describe the purpose of our optimization.
%
%Usage:
%	eq_rec_design_tf( ad , 'Feasible Set' , M! , M2 , T )	
%	eq_rec_design_tf( ad , 'Min_M2' , M1 , T )


% Constants
n = size(ad.A,1);
m = size(ad.B,2);
p = size(ad.C,1);
dd = size(ad.B_d,1);
dm = size(ad.C_m,1);

% Manage Inputs
switch str_in
case 'Feasible Set'

case 'Min_M2'

otherwise
	error(['Unrecognized String: ' str_in ] )

end


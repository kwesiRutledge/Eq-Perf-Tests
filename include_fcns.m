function [] = include_fcns(varargin)
% Description:
%	Includes different libraries, if the user did not properly include the right toolboxes.
%
% Inputs:
%	The inputs will be the desired toolboxes to look for. If nargin ==0 return an error.


if nargin == 0
	error('No toolboxes/libraries given to include.')
end

for arg_ind = 1:nargin
	%Depending on the Toolbox, do different things.
	switch varargin{arg_ind}
		case 'MPT3'
			try
				Polyhedron();
			catch
				warning('MPT3 Toolbox does not appear in path! Attempting one method to fix this...')
				try
					%Add All MPT Toolbox folders AND subfolders to path.
					addpath(genpath([ '../../' 'tbxmanager/toolboxes/mpt/3.1.8/all/mpt3-3_1_8/mpt']) )
					Polyhedron()
				catch
					error('MPT3 not added to path.')
				end

			end
		case 'YALMIP'
			try
				a = sdpvar(1,1,'full');
			catch
				warning('YALMIP does not appear in path! Attempting one method to fix this...');
				try
					addpath(genpath([ '../../' 'YALMIP-master' ]));
					b = sdpvar(1,1,'full');
				catch
					error('YALMIP not added to path.')
				end
			end

		otherwise
			error(['The toolbox name that you provided ' varargin{arg_ind} ' is not a valid library/toolbox. ' ])
	end
end

disp(['All Libraries successfully added to path.'])

end
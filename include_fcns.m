function [] = include_fcns(varargin)
% Description:
%	Includes different libraries, if the user did not properly include the right toolboxes.
%
% Inputs:
%	The inputs will be the desired toolboxes to look for. If nargin ==0 return an error.


if nargin == 0
	warning('No toolboxes/libraries given to include.')
end

if strcmp(getenv('USER'),'kwesirutledge') %Suggests the laptop is in use

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
                        cd([ '../../' 'toolboxes/tbxmanager/toolboxes/mpt/3.1.8/all/mpt3-3_1_8/mpt'])
                        mpt_init
                        Polyhedron();
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
                        addpath(genpath([ '../../' 'toolboxes/YALMIP-master' ]));
                        b = sdpvar(1,1,'full');
                    catch
                        error('YALMIP not added to path.')
                    end
                end
            case 'tbxmanager'
                try
                    addpath(genpath(['../../' 'toolboxes/tbxmanager']) )
                    Polyhedron();
                catch
                    error('tbxmanager was not added to path.')
                end
            case 'mosek'
                try
                    mosekdiag
                catch
                    try
                        addpath(genpath('~/External_Libraries/mosek/9.2/toolbox/r2015a/'))
                        mosekdiag
                    catch
                        error('Mosek was not properly added to path.')
                    end
                end
            case 'gurobi'
                if (exist('gurobi.m') == 2 )
                    %Gurobi is already on the path continue
                    disp('gurobi.m exists!')
                    1;
                else
                    addpath(genpath('/Library/gurobi811/mac64/matlab/'))
                    %Try to find it again
                    if exist('gurobi.m') ~= 2
                        %If the gurobi function still cannot be found,
                        %then produce an error.
                        error('gurobi is still not on path!')
                    end
                end

            otherwise
                error(['The toolbox name that you provided ' varargin{arg_ind} ' is not a valid library/toolbox. ' ])
        end
    end
    disp(['All Libraries successfully added to path.'])
else
    %Depending on the Toolbox, do different things.
    for arg_ind = 1:nargin
    switch varargin{arg_ind}
        case 'MPT3'
            try
                Polyhedron();
            catch
                warning('MPT3 Toolbox does not appear in path!')

            end
        case 'YALMIP'
            try
                a = sdpvar(1,1,'full');
            catch
                warning('YALMIP does not appear in path!');
                try
                    addpath(genpath([ '../../' 'YALMIP-master' ]));
                    b = sdpvar(1,1,'full');
                catch
                    error('YALMIP not added to path.')
                end
            end
    end
    end
    
    warning([ 'Unrecognized user: ' getenv('USER') '. Not adding libraries to path.' ] )
end


end
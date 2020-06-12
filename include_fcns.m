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
elseif (strcmp(getenv('USER'),'krutledg') && isunix) 
    %If the user variable is that of Kwesi's uniqname and this is a Unix machine,
    %then the program is running on the Great Lakes cluster.
    %Use directory paths/structure for that machine instead of Kwesi's own structures.

    path_to_top_of_gl = '~';

    for arg_ind = 1:nargin
        switch varargin{arg_ind}
            case 'MPT3'
                %% Include MPT3
                try
                    disp('- Testing MPT3')
                    Polyhedron();
                    disp('  + MPT3 is already included.')
                catch
                    disp('  + MPT3 Toolbox does not appear in path! Attempting one method to fix this...')
                    try
                        %Add All MPT Toolbox folders AND subfolders to path.
                        addpath(genpath([ path_to_top_of_gl '/matlab_toolboxes/tbxmanager/']) )
                        Polyhedron()
                        disp('  + Successfully added MPT3 to path.')
                    catch
                        error('MPT3 not added to path.')
                    end
                end
            case 'YALMIP'
                %Include YALMIP
                try
                    disp('- Testing YALMIP')
                    sdpvar(1,1,'full');
                    disp('  + YALMIP is already included in the path.')
                catch    
                    disp('  + YALMIP is not currently included on the path! Attempting a fix... ')
                    try
                        %Add 
                        addpath(genpath([ path_to_top_of_gl '/matlab_toolboxes/YALMIP-Master']) )
                        sdpvar(1,1,'full');
                        disp('  + Successfully added YALMIP to path.')
                    end
                end
            case 'tbxmanager'
                %Include YALMIP
                try
                    disp('- Testing tbxmanager')
                    tbxmanager show installed
                    disp('  + tbxmanager is already included in the path.')
                catch    
                    disp('  + tbxmanager is not currently included on the path! Attempting a fix... ')
                    try
                        %Add 
                        addpath(genpath([ path_to_top_of_gl '/matlab_toolboxes/tbxmanager/']) )
                        tbxmanager show installed
                        disp('  + Successfully added tbxmanager to path.')
                    end
                end
            otherwise
                error(['Unrecognized input to include_fcns (' varargin{arg_ind} ').' ] )
        end
    
    end

    disp('Completed adding libraries to the path.')

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
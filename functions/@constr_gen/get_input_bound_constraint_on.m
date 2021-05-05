function [ constraints_out ] = get_input_bound_constraint_on( varargin )
	%get_input_bound_constraint_on.m
	%Description:
	%	The 
	%
	%Usage:
	%	[ constraints_out ] = cg.get_input_bound_constraint_on( lcsas , K , k , 'fb_type' , 'state-disturbance' )
	%	[ constraints_out ] = cg.get_input_bound_constraint_on( lcsas , K , k , 'fb_type' , 'state-disturbance' , 'ActiveGains' , active_gain_vector )
	%	[ constraints_out ] = cg.get_input_bound_constraint_on( lcsas , F , f , 'fb_type' , 'output' )

	%% Input Processing

	[ cg , lcsas , LinearGainArray , OpenLoopGainArray , gibc_settings ] = ip_get_input_bound_constraint_on(varargin{:});

	num_gains = length(LinearGainArray);
	TimeHorizon = length(lcsas.L.words{1});

	%% Algorithm

	% Create PuT
	PuT = 1;
	for t = 1:TimeHorizon
		PuT = PuT * lcsas.U;
	end

	% Create Constraint
	constraints_out = [];

	switch gibc_settings.fb_type
		case 'state-disturbance'

			K = LinearGainArray;
			k = OpenLoopGainArray;

			for gain_index = 1:num_gains
				if gibc_settings.NoLinearGainsProvided
					%Create constraints
					if gibc_settings.ActiveGains(gain_index)
						constraints_out = constraints_out + [ PuT.A * k{gain_index} <= PuT.b ];
					end
				else
					error(['Haven''t written this part of the get_input_bound_constraint_on() function!'])
				end

			end
		otherwise
			error(['Unexpected feedback type should never reach this part of get_input_bound_constraint_on(). What went wrong?'])
	end

	

end

function [ cg , lcsas , LinearGainArray , OpenLoopGainArray , gibc_settings ] = ip_get_input_bound_constraint_on(varargin)
	%Description:
	%	Processes the inputs given to get_input_bound_constraint_on().
	%
	%Usage:
	%	[ cg , lcsas , LinearGainArray , OpenLoopGainArray , gibc_settings ] = ip_get_input_bound_constraint_on(varargin{:})

	%% Constants

	allowable_fb_types = {'state-disturbance'};

	%% Algorithm

	% Check that there is the minimum number of necessary inputs
	if nargin < 4
		error(['get_input_bound_constraint_on() requires at least 4 arguments. Received ' num2str(nargin) '.'])
	end

	cg 					= varargin{1};
	lcsas 				= varargin{2};
	LinearGainArray		= varargin{3};
	OpenLoopGainArray 	= varargin{4};

	% Checking the Necessary Arguments
	if ~isa(lcsas,'LCSAS')
		error(['Expected the second argument of get_input_bound_constraint_on() to be of class LCSAS. Received an object of class ' class(lcsas) '.' ])
	end

	lcsas.check('U');

	NoLinearGainsProvided = true;
	for lin_gain_index = 1:length(LinearGainArray)
		if isa(LinearGainArray{lin_gain_index},'sdpvar')
			NoLinearGainsProvided = false;
			break;
		end
	end

	% Create Default Settings Structure
	num_gains = length(LinearGainArray);

	gibc_settings = struct( ...
		'NoLinearGainsProvided',NoLinearGainsProvided, ...
		'fb_type','state-disturbance', ...
		'ActiveGains',true(1,num_gains));

	% Checking for Optional Arguments
	arg_index = 5;
	while arg_index <= nargin
		switch varargin{arg_index}
			case 'fb_type'
				gibc_settings.fb_type = varargin{arg_index+1};
				arg_index = arg_index + 2;
			case 'ActiveGains'
				gibc_settings.ActiveGains = varargin{arg_index+1};
				arg_index = arg_index + 2;
			otherwise
				error(['Unexpected input to get_input_bound_constraint_on(): ' varargin{arg_index} '.' ])
		end
	end


	% Checking The Settings
	if ~any(strcmp(gibc_settings.fb_type,allowable_fb_types))
		error(['Unexpected fb_type value provided(' gibc_settings.fb_type '). Please see the variable allowable_fb_types for more info on allowable types.'])
	end

	if num_gains ~= length(gibc_settings.ActiveGains)
		error(['There number of gains (' num2str(num_gains) ') does not match the ActiveGains setting''s vector length (' num2str(length(gibc_settings.ActiveGains)) ').'  ])
	end




end
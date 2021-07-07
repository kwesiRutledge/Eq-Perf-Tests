function [ lcsas0 , KnowledgeSequence , ibs_settings ] = ip_InternalBehaviorSet( varargin )
	%Description:
	%	Input processing function for internal_behavior_set().
	%	Verifies that there:
	%	- are at least four inputs,
	%	- LCSAS input checks,

	%% Check Argument Numbers
	expected_number_of_arguments = 3;

	if nargin < expected_number_of_arguments
		error(['Expect at least four arguments.'])
	end

	%% First Required Input: lcsas
	ibs = varargin{1};
	lcsas0 = varargin{2};
	lcsas0.check('X0','U'); %Verify that X0 and U exist and are both Polyhedron objects.

	%% Second Required Input: KnowledgeSequence
	KnowledgeSequence = varargin{3};
	if ~isa(KnowledgeSequence,'Language')
		error(['KnowledgeSequence  is not of type ''Language''. Instead it is of type ' class(KnowledgeSequence) '.'])
	end

	if ~isvector(KnowledgeSequence)
		error(['KnowledgeSequence must be a vector or scalar of Language objects.'])
	end

	%% Define Defaults
	ibs_settings = struct('fb_type','state','reduce_flag',false,'debug',0,'ReturnEarly',false, ...
							'A', [] , 'b' , [] , 'Ae' , [] , 'be' , [] , ...
							'Dim', -1 );

	all_fb_types = {'state','output'};
	all_fields = {'fb_type','reduce_flag','debug','ReturnEarly','A','b','Ae','be','Dim'};
	expected_ebs_fields = {'fb_type','reduce_flag'};

	%% Check the additional inputs
	argument_index = expected_number_of_arguments+1;
	while argument_index <= nargin
		switch varargin{argument_index}
			case 'fb_type'
				ibs_settings.fb_type = varargin{argument_index+1};
				if ~any(strcmp(ibs_settings.fb_type,all_fb_types))
					error(['Unexpected fb_type value: ' ibs_settings.fb_type ])
				end
				argument_index = argument_index + 2;
			case 'reduce_flag'
				ibs_settings.reduce_flag = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'ReturnEarly'
				ibs_settings.ReturnEarly = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'A'
				ibs_settings.A = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'b'
				ibs_settings.b = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'Ae'
				ibs_settings.Ae = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'be'
				ibs_settings.be = varargin{argument_index+1};
				argument_index = argument_index + 2;
			case 'ibs_settings_struct'
				ibs_settings = varargin{argument_index+1};
				argument_index = argument_index + 2;

				fields_available = isfield(ibs_settings,all_fields);
				if ~all( fields_available )
					first_missing_field = all_fields{ find(fields_available,1) };
					error(['The field ' first_missing_field ' is not available in ibs_settings. Please create it!' ])
				end
			case 'ebs_settings_struct'
				ebs_settings = varargin{argument_index+1};
				argument_index = argument_index + 2;

				fields_available = isfield(ebs_settings,expected_ebs_fields);
				if ~all( fields_available )
					first_missing_field = all_fields{ find(fields_available,1) };
					error(['The field ' first_missing_field ' is not available in ibs_settings. Please create it!' ])
				end

				ibs_settings.fb_type = ebs_settings.fb_type;
				ibs_settings.reduce_flag = ebs_settings.reduce_flag;

			otherwise
				error(['Unexpected input to internal_behavior_set: ' varargin{argument_index} ])
		end
	end

	% Check inputs when the ReturnEarly flag is raised.
	if ibs_settings.ReturnEarly
		if isempty(ibs_settings.A) && isempty(ibs_settings.Ae)
			error(['The InternalBehaviorSet constructed using ReturnEarly must contain a form of H constraint.'])
		end

		tempDimArray = [ size(ibs_settings.A,2) ; size(ibs_settings.Ae,2) ];
		if all(tempDimArray > 0) && ( tempDimArray(1) ~= tempDimArray(2) )
			error(['The dimension of the set according to A is ' num2str(size(ibs_settings.A,2)) ' while the dimension of the set according to Ae is ' num2str(size(ibs_settings.Ae,2)) '. Please adjust your input matrices.'  ])
		end 

		ibs_settings.Dim = max(tempDimArray);

	end
end
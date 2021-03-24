function [ PwT ] = find_hypothesis_generating_disturbances( varargin )
	%Description:
	%	Finds the set of all disturbances which 	
	%Usage:
	%	[ PwT ] = lcsas0.find_hypothesis_generating_disturbances( T , L , Px0 , Pu )


	%% Input Processing

	[ lcsas_in , T , L , Px0 , Pu , fb_type ] = find_hypothesis_generating_disturbances_input_processing(varargin{:});

	%% Constants

	card_L = L.cardinality();
	[ n_x , n_u , n_y , n_w , n_v ] = lcsas_in.Dimensions();

	%% Algorithm

	[ ~ , PhiI ] = lcsas_in.consistent_set( T , L , Pu , Px0 );

	PwT = [];
	for word_idx = 1:L.cardinality()
		PwT = [ PwT ; PhiI.projection(n_x*(T+1) + n_u*(T) +[1:n_w*(T)]) ];
	end

end

function [ lcsas_out , T , L , Px0 , Pu , fb_type ] = find_hypothesis_generating_disturbances_input_processing(varargin)
	%Description:
	%	Manages the inputs appropriately.

	%%%%%%%%%%%%%%%
	%% Constants %%
	%%%%%%%%%%%%%%%

	%% Algorithm

	if nargin < 5
		error(['The function find_hypothesis_generating_disturbances() requires at least five input arguments. Received ' num2str(nargin) '.'])
	end

	% Assigning required values

	lcsas_out = varargin{1};
	T = varargin{2};
	L = varargin{3};
	Px0 = varargin{4};
	Pu = varargin{5};

	if ~isa(lcsas_out,'LCSAS')
		error(['The first input is not of type LCSAS. Received type ' class(lcsas_out) '.' ] )
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Assign Default Values for All Other Things %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	fb_type = 'state'; %State feedback option.

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Process Any Extra Inputs %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if nargin > 5
		argidx = 6;
		while argidx <= nargin
			switch varargin{argidx}
				case 'Feedback Type'
					fb_type = varargin{argidx+1};
					argidx = argidx + 2;
				otherwise
					error(['Unexpected input to find_hypothesis_generating_disturbances: ' varargin{argidx} ])
			end
		end
	end
end
function [ observable_ibs_flags ] = FindObservableIBS(varargin)

	%% Input Processing
	[ ibs_array , empty_flags ] = ip_FindObservableIBS(varargin{:});

	%% Constants
	num_ibs = length(ibs_array);

	%% Algorithm
	if isempty(empty_flags)
		empty_flags = ibs_array.IsEmpty(); %Find which ibs are empty.
	end

	observable_ibs_flags = ~empty_flags;

	for ibs_index1 = 1:num_ibs
		
		if empty_flags(ibs_index1) %If the ibs is empty, then ignore it in our considerations.
			continue; 
		end

		ibs1 = ibs_array(ibs_index1);

		for ibs_index2 = ibs_index1+1:num_ibs
			%Collect ibs.
			ibs2 = ibs_array(ibs_index2);

			if ibs2.ExternalBehaviorCovers(ibs1)
				observable_ibs_flags(ibs_index1) = false;
			end

		end
	end

end

function [ ibs_array , empty_flags ] = ip_FindObservableIBS(varargin)
	%Description:
	%	Handles the inputs to FindObservableIBS

	%% Algorithm

	ibs_array = varargin{1};

	% Initialize Variables
	num_ibs = length(ibs_array);

	empty_flags = [];

	varargin_idx = 2;
	while varargin_idx <= nargin
		switch varargin{varargin_idx}
		case 'empty_flags'
			empty_flags = varargin{varargin_idx+1};
			varargin_idx = varargin_idx + 2;



		otherwise
			error(['Unexpected input to FindObservableIBS: ' varargin{varargin_idx} ])
		end

	end


	% Check Inputs
	if ~isempty(empty_flags)
		if length(empty_flags) ~= num_ibs
			error(['The input empty flags has length ' num2str(length(empty_flags)) ', but there are ' num2str(num_ibs) ' InternalBehaviorSet objects in the input ibs_array.'  ])
		end
	end

end
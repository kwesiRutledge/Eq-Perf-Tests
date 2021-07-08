function [ observable_ebs_flags ] = FindObservableEBS(varargin)
	%Description:
	%
	%Usage:
	%	observable_ebs_flags = ebs_array.FindObservableEBS()
	%	observable_ebs_flags = ebs_array.FindObservableEBS('empty_flags',empty_flags_in)


	%% Input Processing
	[ ebs_array , empty_flags ] = ip_FindObservableEBS(varargin{:});

	%% Constants
	num_ebs = length(ebs_array);

	%% Algorithm
	if isempty(empty_flags)
		empty_flags = ebs_array.IsEmpty(); %Find which ibs are empty.
	end

	observable_ebs_flags = ~empty_flags;

	for ebs_index1 = 1:num_ebs
		
		if empty_flags(ebs_index1) %If the ibs is empty, then ignore it in our considerations.
			continue; 
		end

		ebs1 = ebs_array(ebs_index1);

		for ebs_index2 = ebs_index1+1:num_ebs
			%Collect ibs.
			ebs2 = ebs_array(ebs_index2);

			if ebs2 >= ebs1
				observable_ebs_flags(ebs_index1) = false;
			end

		end
	end

end

function [ ebs_array , empty_flags ] = ip_FindObservableEBS(varargin)
	%Description:
	%	Handles the inputs to FindObservableIBS

	%% Algorithm

	ebs_array = varargin{1};

	% Initialize Variables
	num_ebs = length(ebs_array);

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
		if length(empty_flags) ~= num_ebs
			error(['The input empty flags has length ' num2str(length(empty_flags)) ', but there are ' num2str(num_ebs) ' InternalBehaviorSet objects in the input ebs_array.'  ])
		end
	end

end
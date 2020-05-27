#! /bin/bash
#
#load_modules_gl.sh
#Description:
#	This shell script loads the following modules:
#	- matlab/R2019b
#	- gcc
#	- gurobi

 #try
echo "Loading MATLAB..."
module load matlab || { #catch
	echo "  (There is an issue attempting to load matlab.)";
	echo "$module_warning";
	return;
}
echo "  + MATLAB loaded."


echo "Loading gcc..."
#try
module load gcc || { #catch
	echo "  (There is an issue attempting to load gcc.)";
	echo "$module_warning";
	return;
}
echo "  + gcc loaded."

echo "Loading gurobi..."
module load gurobi || { #catch
	echo "  (There is an issue attempting to load gurobi."
	echo "$module_warning";
	return;
}
echo "  + gurobi loaded."
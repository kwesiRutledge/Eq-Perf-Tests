#! /bin/bash
# A simple shell script to:
#   - Login to the great lakes cluster

echo " "
echo "login_gl.sh"
echo " "
echo "Warning: This shell script is designed for MAC OS X."
echo "Results not guaranteed for other operating systems."
echo " "

echo "Opening UM VPN..."
networksetup -connectpppoeservice "UMVPN"
echo "  + VPN opened"
echo " "

ssh krutledg@greatlakes.arc-ts.umich.edu


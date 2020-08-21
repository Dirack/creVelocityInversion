#!/bin/bash
#
# ascFormat.sh (Shell Script)
# 
# Purpose: Convert sfipick picked points
# to Madagascar readable ascii format.
# Important! $1=1 cut and format the first collumn of the file (t0 time)
# If $1=2 cut and format the second collumn of the file (m0 CMP)
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programer: Rodolfo A C Neves (Dirack) 20/08/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

INPUT=$(cat "-" | tr '\t' ' ' | cut -d' ' -f"$1")
NUMBER_OF_POINTS=$(echo "$INPUT" | wc -l)
FORMATED_INPUT=$(echo "$INPUT" | tr '\n' ' ')

if [ "$1" == "1" ]
then
	FILE="t0s"
elif [ "$1" == "2" ]
then
	FILE="m0s"
else
	echo "Error: $(basename $0): \$1 should be equal to 1 or 2!"
	exit 1
fi

echo "${FORMATED_INPUT}n1=$NUMBER_OF_POINTS d1=1 o1=0 data_format=ascii_float in=$2"

INPUT=$(cat "-" | tr '\t' ' ' | cut -d' ' -f"$2")
NUMBER_OF_POINTS=$(echo "$INPUT" | wc -l)
FORMATED_INPUT=$(echo "$INPUT" | tr '\n' ' ')
FILE="$1"

echo "${FORMATED_INPUT}n1=$NUMBER_OF_POINTS d1=1 o1=0 data_format=ascii_float in=$1"


output_directory=$1
output_report=$2
read_1=$3
read_2=$4
fastp_parameters=$6
pairedend=$7

if [ "$pairedend" = true ]; then

    read_1_out=${output_directory}$(basename $3)
    read_2_out=${output_directory}$(basename $4)
    html_report=${output_report}_fastp.html
    json_report=${output_report}_fastp.json

    command="fastp -i $read_1 -I $read_2 \
    -o $read_1_out -O $read_2_out \
    "$fastp_parameters" -h $html_report \
    -j $json_report"
    $command
else

    read_out=${output_directory}$(basename $3)
    html_report=${output_report}_fastp.html
    json_report=${output_report}_fastp.json

    command="fastp -i $read_1 \
    -o $read_out \
    "$fastp_parameters" -h $html_report \
    -j $json_report"
    $command

fi
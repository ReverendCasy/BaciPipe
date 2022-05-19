#!/usr/bin/zsh

eval set -- `getopt --options r:l:s:i:c:p:m:o: --longoptions read_table:,longread_dir:,shortread_dir:,polish_rounds:,config_file:,proteins:,prodigal_model:,out_dir: -- $@`

while true; do
    case ${1} in
        -r|--read_table)
            read_table=${2}
            shift 2 ;;
        -l|--longread_dir)
            longread_dir=${2}
            shift 2 ;;
        -s|--shortread_dir)
            shortread_dir=${2}
            shift 2 ;;
        -i|--polish_rounds)
            polish_rounds=${2}
            shift 2 ;;
        -c|--config_file)
            config_file=${2}
            shift 2 ;;
        -p|--proteins)
            proteins=${2}
            shift 2 ;;
        -m|--prodigal_model)
            prodigal_model=${2}
            shift 2 ;;
        -o|--out_dir)
            out_dir=${2}
            shift 2 ;;
        --)
            break
            shift ;;
        *)
            exit 1 ;;
    esac
done

polish_rounds=${polish_rounds:-1}

strains=( $(cut -f1 ${read_table} | awk 'NR > 1' | sort | uniq) )
declare -A long_reads=()
declare -A forwards=()
declare -A reverses=()

for strain in ${strains[@]}; do
    long_barcode=$(awk -v strain=${strain} -F'\t' '($1==strain && tolower($4)=="nanopore"){print $3}' ${read_table})
    echo ${long_barcode}
    if ((${long_barcode} < 10)); then
        long_barcode="0${long_barcode}"
    fi

    longread_path=$(find ${longread_dir} -type f -name "barcode${long_barcode}_*")
    long_reads[${strain}]=${longread_path}

    illumina_barcode=$(awk -v strain=${strain} -F'\t' '($1==strain && tolower($4)=="illumina"){print $3}' ${read_table})
    echo ${illumina_barcode}
    if ((${illimina_barcode} <  10)); then
        illumina_barcode="0${illumina_barcode}"
    fi

    forward_path=$(find ${shortread_dir} -type f -name "${illumina_barcode}_f*")
    reverse_path=$(find ${shortread_dir} -type f -name "${illumina_barcode}_r*")
    forwards[${strain}]=${forward_path}
    reverses[${strain}]=${reverse_path}

done

echo "strains:" > ${config_file}
for strain in ${strains[@]}; do
    echo "    ${strain}: ${long_reads[${strain}]}" >> ${config_file}
done

echo "forward:" >> ${config_file}
for strain in ${strains[@]}; do
    echo "    ${strain}: ${forwards[${strain}]}" >> ${config_file}
done

echo "reverse:" >> ${config_file}
for strain in ${strains[@]}; do
    echo "    ${strain}: ${reverses[${strain}]}" >> ${config_file}
done

echo "polish_rounds: ${polish_rounds}" >> ${config_file}

if [[ ! -z ${proteins} ]]; then
    echo "reference_proteins: ${proteins}" >> ${config_file}
fi

if [[ ! -z ${prodigal_model} ]]; then
    echo "prodigal_model: ${prodigal_model}" >> ${config_file}
fi

echo "\nout_dir: ${out_dir:-pipeline_out}" >> ${config_file}

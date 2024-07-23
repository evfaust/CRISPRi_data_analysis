#!/bin/bash
############################################################
# Help                                                     #
############################################################
Help()
{
    # Display Help
    echo "Process raw paired-end sequence files (ending with _R1.fastq.gz or _R2.fastq.gz) from a CRISPRi screen
(unzip, merge, convert to fasta, align, count reads, and combine samples into single dataframe all_counts.tsv)."
    echo
    echo "Example with minimum inputs: ./process_data.sh -r /path/to/raw/seqs -d /path/to/library_database -a /path/for/out_alignments -c /path/for/out_counts -w wash_control_name"
    echo "Syntax: scriptTemplate [-h|r|o|d|i|m|t]"
    echo "options:"
    echo "h     Print this help."
    echo "r     Enter the folder path storing raw paired-end .fastq.gz (zipped) sequence files. Ensure paired files have the same prefix and end with '_R1.fastq.gz' or _R2.fastq.gz'. 
Merged .fasta sequence files will be stored here."
    echo "d     Enter the database .fasta file as alignment reference. Must contain the reverse complements of guides."
    echo "o     Enter the desired output folder for storing all future output files, including alignment .tsv files. Script will create folder if it does not exist."
    echo "i     Enter the percent identity for alignment as a fraction. (Default is 0.9.)"
    echo "m     Enter the minimum sequence length for vsearch alignment. (Default is 1.)"
    echo "t     Enter the target coverage for vsearch alignment. (Default is 1.)"
    echo
    echo "Note: it is also recommended to look at read quality using softwares such as fastqc. Expect fails in ###. 
    Q score 30 is considered very good, and N ratio should not exceed 5% (lack of confidence at a certain position - indicates machine malfunction)"
    echo
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

# Set default variables
Align_Percent=0.9 # default align percent
Min_Seq_Length=1 # default minimum sequence length
Target_Coverage=1 # default target coverage

############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while getopts "hr:o:d:i:m:t:c:w:" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        r) # Enter the raw sequence folder
            Raw_Seq_Folder=$OPTARG;;
        o) # Enter the desired output folder for storing alignment .tsv files and future outputs
            Output_Folder=$OPTARG;;
        d) # Enter the database file
            Database_File=$OPTARG;;
        i) # Enter the percent identity for alignment
            Align_Percent=$OPTARG;;
        m) # Enter the minimum sequence length
            Min_Seq_Length=$OPTARG;;
        t) # Enter the target coverage
            Target_Coverage=$OPTARG;;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

# Unzip any files in the Raw_Seq_Folder ending in .fastq.gz
if ls ${Raw_Seq_Folder}/*.fastq.gz 1> /dev/null 2>&1; then
    gunzip ${Raw_Seq_Folder}/*.fastq.gz
fi

# Create unique_samples array containing the prefixes of the paired fastq files
unique_samples=($(ls $Raw_Seq_Folder | grep -E 'R[1-2].fastq$' | cut -d'_' -f1 | sort -u))

# Create Output_Folder directory if it does not already exist
if [ ! -d "$Output_Folder" ]; then
    mkdir -p "$Output_Folder"
fi

# For each unique sample, merge its paired fastq files, convert to fasta, align to the library database, and calculate counts per guide
for sample in ${unique_samples[@]}; 
do
    # merge and convert to fasta
    vsearch --fastq_mergepairs ${Raw_Seq_Folder}/${sample}_R1.fastq --reverse ${Raw_Seq_Folder}/${sample}_R2.fastq --fastaout ${Raw_Seq_Folder}/${sample}_merged.fasta
    echo "hi"
    # align to provided database guide library file
    vsearch --db $Database_File --strand both --id $Align_Percent --top_hits_only --target_cov $Target_Coverage --userfields query+target+id+qrow+trow --minseqlength $Min_Seq_Length --usearch_global ${Raw_Seq_Folder}/${sample}_merged.fasta --userout ${Output_Folder}/${sample}_aligned.tsv
    echo "hi2"
    # calculate counts per guide
    #python3 counter.py -f ${Output_Folder}/${sample}_aligned.tsv -d $Database_File -o ${Output_Counts_Folder}/${sample}_counts.tsv -i 100.0 -w $Wash_Control
done


# Combine counts .tsv files into one .tsv file
#python3 combine_count_files.py -c $Output_Counts_Folder







# TODO
# store statistics/diagnostics of merging and aligning somewhere
# options for strandedness, top_hits_only, etc?
# options for filtering by average Q value or trimming non-overlapping sections?
# user fields options? also put in help for output
# output counts into the same .tsv file column at a time


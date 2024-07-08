#!/bin/bash
############################################################
# Help                                                     #
############################################################
Help()
{
    # Display Help
    echo "Preprocess raw paired-end sequence files from a CRISPRi screen
(merge, convert to fasta, and align)."
    echo
    echo "Syntax: scriptTemplate [-h|r|o|d|i|m|t]"
    echo "options:"
    echo "h     Print this help."
    echo "r     Enter the folder path storing raw paired-end .fastq.gz (zipped) sequence files. Ensure paired files have the same prefix and end with '_R1.fastq.gz' or _R2.fastq.gz'. 
Merged .fasta sequence files will be stored here."
    echo "d     Enter the database .fasta file as alignment reference. Must contain the reverse complements of guides."
    echo "o     Enter the desired output folder for storing alignment .tsv files. Script will create folder if it does not exist."
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
while getopts "hr:o:d:i:m:t:" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        r) # Enter the raw sequence folder
            Raw_Seq_Folder=$OPTARG;;
        o) # Enter the desired output folder for storing alignment .tsv files
            Output_Alignments_Folder=$OPTARG;;
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

# Create Output_Alignments_Folder directory if it does not already exist
if [ ! -d "$Output_Alignments_Folder" ]; then
    mkdir -p "$Output_Alignments_Folder"
fi

# For each unique sample, merge its paired fastq files, convert to fasta, and align to the library database
for sample in ${unique_samples[@]}; 
do
    echo "hi"
    # merge and convert to fasta
    vsearch --fastq_mergepairs ${Raw_Seq_Folder}/${sample}_R1.fastq --reverse ${Raw_Seq_Folder}/${sample}_R2.fastq --fastaout ${Raw_Seq_Folder}/${sample}_merged.fasta
    # align to provided database guide library file
    vsearch --db $Database_File --strand both --id $Align_Percent --top_hits_only --target_cov $Target_Coverage --userfields query+target+id+qrow+trow --minseqlength $Min_Seq_Length --usearch_global ${Raw_Seq_Folder}/${sample}_merged.fasta --userout ${Output_Alignments_Folder}/${sample}.tsv
done


# TODO
# store statistics/diagnostics of merging and aligning somewhere
# options for strandedness, top_hits_only, etc?
# options for filtering by average Q value or trimming non-overlapping sections?
# user fields options? also put in help for output


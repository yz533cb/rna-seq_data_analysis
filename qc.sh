# The quality of raw reads was determined using FastQC and MultiQC packages
for f in "directory/raw_data"/*
do
    if [[ $f =~ * ]]; then
        for data in "$f"/*
        do
            if [[ $data =~ *gz ]]; then
# the file name of raw data ends with “.gz”
                fastqc -o directory/QC "$data" 
            fi
        done
    fi
done

# Visualise the QC results using MultiQC package written with Python
Python3 -m multiqc

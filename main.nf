#!/usr/bin/env nextflow

// Version 0.0

nextflow.enable.dsl = 1

params.inputFile = "$baseDir/Disteche_RNA3-041-042_nova.txt"
params.outputDir = "$baseDir/results"
params.run = "RNA3-041-042_nova"
params.ref_fasta ="/net/bbi/vol1/data/genomes_stage/human/human_star/Homo_sapiens.GRCh38.dna.toplevel.fa.finished"
params.vcf = "/net/bbi/vol1/data/sciRNAseq/benchmarking/resources/human_ref_vcf/common_variants_grch38.vcf"
params.threads = 8
params.k_range = Channel.from(21,22,23,24,25)  // Range of k-values to run SoupOrCell process in parallel 

// Read in a .txt file that contains paths to sample directories 
path_list = Channel.fromPath(params.inputFile)
            .splitText()    // split at each line 
            .map {it.trim()}    // trim any trailing new lines

// Make copy of sample paths for bam and barcode files 
path_list.into{path_list_copy01; path_list_copy02}


/***

Process: Modify bam files to have appropriate tags to run SoupOrCell 

Inputs:  Sample .bam files 

Published: None

Downstream: mergeBams

Notes: 

***/

process processSample {
    cache 'lenient'
    executor 'sge'
    memory '16 GB'
    queue = "shendure-long.q"

    
    // publishDir  path: "${params.outputDir}/", pattern: "*.bam", mode: 'copy'

    input:
        val sample from path_list_copy01

    output:
        file("*.bam") into processed_bams

    script:
    def baseName = new File(sample).name
    key = baseName
    """
    module load python/3.12.1
    python /net/bbi/vol1/data/sciRNAseq/genetic_demultiplexing/reformat_sci_python3_test.py \
        -i "${sample}" \
        -r "${params.run}" \
        -o "${baseName}_socReady"
    """
}

/***

Process: Create a merged bam file from modified bams in previous process. 

Inputs:  - *.bams - all modified bams from each sample 

Published: Merged bam and index bam files in .bam and .bai format. 

Downstream: runSoupOrCell

Notes: Merged bams are saved to a temp folder that will be used for downstream soc analysis 

***/

// save_merged_bam_temp = {params.outputDir + "/" + it - ~/merged.bam/ - ~/merged/ + "/" + it - ~/merged/}
// save_merged_bai_temp = {params.outputDir + "/" + it - ~/merged.bam.bai/ - ~/merged/ + "/" + it - ~/merged/}

process mergeBams {
    cache 'lenient'
    executor 'sge'
    memory '64 GB'
    queue = "shendure-long.q"

    publishDir path: "${params.outputDir}/", pattern: "merged", mode: 'copy'

    input: 
        file input_bam from processed_bams.collect()

    output: 
        file("merged") into bam_dir 

    script:

    """
    # bash watch for errors
    set -ueo pipefail

    module load samtools/1.19  

    echo "testing"
    mkdir merged
    samtools merge -o "merged/${params.run}_merged.bam" ${input_bam} 
    samtools sort "merged/${params.run}_merged.bam" > "merged/${params.run}_merged_sorted.bam"
    samtools index "merged/${params.run}_merged_sorted.bam"
    """
}

/***

Process: Get sample barcodes and filter for => 100 umis. 
         Append run name to barcode as found in modified bam tags. 

Input: sample paths 

Output: Filtered barcodes tsv with => 100 umis and run name appended to barcode 

Downstream: runSoupOrCell

Published: None

Notes:
***/

process getBarcodes{
    cache 'lenient'
    executor 'sge'
    memory  '16 GB'
    queue = "shendure-long.q"

    input:
        val sample_path from path_list_copy02

    output: 
        file("*.tsv") into barcodes

    script:
        def baseName = new File(sample_path).name
    """
    # bash watch for errors
    set -ueo pipefail

    id=\$(md5sum "${sample_path}/umis_per_cell_barcode.txt" | cut -d' ' -f1)
    awk -v run="${params.run}" '\$2 >= 100 {print \$1"_"run}' "${sample_path}/umis_per_cell_barcode.txt" \
    > "\${id}_${baseName}_barcode.tsv"
    """
}

/***

Process: Merge all sample barcodes with >= 100 umis.

Input: Filtered barcodes with >= 100 umis from each sample.

Output: Merged barcodes .tsv with from each sample. 

Published: Merged barcodes in .tsv format 

Notes: 

***/

process mergeBarcodes{
    cache 'lenient'
    executor 'sge'
    memory  '16 GB'
    queue = "shendure-long.q"

    publishDir path: "${params.outputDir}/merged/", pattern: "*.tsv", mode: 'copy'

    input:
        file input_files from barcodes.collect()

    output: 
        file("*.tsv") into merged_barcodes

    script:
    """
    # bash watch for errors
    set -ueo pipefail
   
    cat ${input_files} > ${params.run}_merged_barcodes.tsv
    """
}

/***

Process: Run SoupOrCell on merged bams files 

Input: - .bam = Merged, sorted, and then indexed bam files. 
       - .tsv = 

Output: Merged barcodes .tsv with from each sample. 

Published: Merged barcodes in .tsv format 

Notes: 

***/

process runSoupOrCell {
    cache 'lenient'
    executor 'sge'
    queue = "shendure-long.q"
    memory  '16 GB'
    cpus = '12'
    penv = 'serial'

    input: 
        file input_bam from bam_dir
        file input_barcode from merged_barcodes
        val kval from params.k_range

    output: 
        file("soc_output_k*") into soc_output
        val kval into kval_out
        
    script: 

    """
    # bash watch for errors
    set -ueo pipefail
    
    mkdir -p soc_output_k"$kval"

    singularity exec \
        --bind /net/bbi/vol1/data/  \
        /net/bbi/vol1/nobackup/apptainer/Demuxafy.sif \
        souporcell_pipeline.py -i ${input_bam}/${params.run}_merged_sorted.bam -b ${input_barcode} -f ${params.ref_fasta} -t ${params.threads} -o soc_output_k${kval} -k ${kval} \
        --skip_remap SKIP_REMAP \
        --common_variants ${params.vcf} \
    """

}

/**

Process: Run SoupOrCell summary for the number of droplets classified as doublets, amiguous and assgiend to each cluster. 

**/

process socSummary {
    cache 'lenient'
    executor 'sge'
    queue = "shendure-long.q"
    memory  '32 GB'

    publishDir path: "${params.outputDir}/", pattern: "soc_output_k*", mode: 'copy'

    input: 
        val kval from kval_out
        file soc_dir from soc_output

    output:
        file soc_dir into soc_summary

    script:    
    
    """
    # bash watch for errors
    set -ueo pipefail

    singularity exec \
        --bind /net/bbi/vol1/data/ \
        /net/bbi/vol1/nobackup/apptainer/Demuxafy.sif \
        bash souporcell_summary.sh \
        ${soc_dir}/clusters.tsv > ${soc_dir}/souporcell_summary.tsv
    """
}
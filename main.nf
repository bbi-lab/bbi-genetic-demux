#!/usr/bin/env nextflow

// Version 0.3

nextflow.enable.dsl = 1

params.inputFile = ""
params.outputDir = ""
params.singularity_bindpath=""
params.singularity_path= "" 
params.output_runName = ""
params.skip_remap = false // if true, need to provide a common variant file in params.vcf
params.ref_fasta =""
params.vcf = "" // leave blank to skip common variant and set params.skip_remap to false.
params.threads = 8


// Range of k-values to run SoupOrCell process in parallel in list format 
// eg. ['2', '3']
params.k_range =  ''
kval_range = Channel.fromList(params.k_range).view()

// Make copy of k-value range for socSummary and runSoupOrCell
kval_range.into{kval_range_copy01; kval_range_copy02}

// Read in a .csv file that contains paths to sample directories and run names
run_info= Channel.fromPath(params.inputFile)
            .splitCsv(header: true, sep: ',')
            .map { row -> tuple( row.run_name, row.sample_path)}
            .view()

// Make copy of sample paths for bam and barcode files 
run_info.into{run_info_copy01; run_info_copy02}


/***

Process: Modify bam files to have appropriate tags to run SoupOrCell 

Inputs:  Sample .bam files 

Published: None

Downstream: mergeBams

Notes: 

***/

process processSample {
    cache 'lenient'
    memory '16 GB'

    input:
        tuple val(run_name), val(sample_path) from run_info_copy01

    output:
        file("*.bam") into processed_bams

    script:
        def sampleName = new File(sample_path).name
    """
    module load python/3.12.1
    reformat-sci-bam.py \
        -i "${sample_path}" \
        -r "${run_name}" \
        -o "${run_name}_${sampleName}_socReady"
    """
}

/***

Process: Create a merged bam file from modified bams in previous process. 

Inputs:  - *.bams - all modified bams from each sample 

Published: Merged bam and index bam files in .bam and .bai format. 

Downstream: runSoupOrCell

Notes: Merged bams are saved to a temp folder that will be used for downstream soc analysis 

***/

process mergeBams {
    cache 'lenient'
    memory '64 GB'

    publishDir path: "${params.outputDir}/merged/", pattern: "*sorted*", mode: 'copy'

    input: 
        file input_bam from processed_bams.collect()

    output: 
        file("*sorted*") into bam_dir 

    script:

    """
    # bash watch for errors
    set -ueo pipefail

    module load samtools/1.19  

    samtools merge -o "${params.output_runName}_merged.bam" ${input_bam} 
    samtools sort "${params.output_runName}_merged.bam" > "${params.output_runName}_merged_sorted.bam"
    samtools index "${params.output_runName}_merged_sorted.bam"
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
    memory  '16 GB'

    input:
        tuple val(run_name), val(sample_path) from run_info_copy02

    output: 
        file("*.tsv") into barcodes

    script:
        def sampleName = new File(sample_path).name
    """
    # bash watch for errors
    set -ueo pipefail

    awk -v run="${run_name}" '\$2 >= 100 {print \$1"_"run}' "${sample_path}/umis_per_cell_barcode.txt" \
    > "${run_name}_${sampleName}_barcode.tsv"
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
    memory  '16 GB'

    publishDir path: "${params.outputDir}/merged/", pattern: "*.tsv", mode: 'copy'

    input:
        file input_files from barcodes.collect()

    output: 
        file("*.tsv") into merged_barcodes

    script:
    """
    # bash watch for errors
    set -ueo pipefail
   
    cat ${input_files} > ${params.output_runName}_merged_barcodes.tsv
    """
}

/***

Process: Run SoupOrCell on merged bams files 

Input: - .bam = Merged, sorted, and then indexed bam files. 
       - .tsv = Merged and filtered barcodes

Output: - alt.mtx 
        - ambient_rna.txt
        - clusters_genotypes.vcf
        - clustering.done
        - clusters_tmp.tsv
        - clusters.tsv
        - common_variants_covered_tmp.vcf
        - common_variants_covered.vcf
        - consensus.done
        - depth_merged.bed
        - doublets.err
        - ref.mtx
        - troublet.done
        - variants.done
        - vatrix.done

Published:  - alt.mtx 
            - ambient_rna.txt
            - clusters_genotypes.vcf
            - clustering.done
            - clusters_tmp.tsv
            - clusters.tsv
            - common_variants_covered_tmp.vcf
            - common_variants_covered.vcf
            - consensus.done
            - depth_merged.bed
            - doublets.err
            - ref.mtx
            - troublet.done
            - variants.done
            - vatrix.done
Notes: 

***/


process runSoupOrCell {
    cache 'lenient'
    memory  '16 GB'
    cpus = '12'
    penv = 'serial'

    input: 
        set file(input_bam), file (input_bai) from bam_dir
        file input_barcode from merged_barcodes
        val kval from kval_range_copy01

    output: 
        file("soc_output_k*") into soc_output
        val kval into kval_out
        
    script: 
    """
    # bash watch for errors
    set -ueo pipefail
    
    mkdir -p soc_output_k"$kval"

    singularity exec \
        --bind ${params.singularity_bindpath}  \
        ${params.singularity_path} \
        souporcell_pipeline.py -i ${input_bam} -b ${input_barcode} -f ${params.ref_fasta} -t ${params.threads} -o soc_output_k${kval} -k ${kval} \
        --skip_remap ${params.skip_remap} \
        --common_variants ${params.vcf} \
    """

}

/**

Process: Run SoupOrCell summary for the number of droplets classified as doublets, amiguous and assigned to each cluster. 

Input: - clusters.tsv 

Output: souporcell_summary.tsv 

Published: souporcell_summary.tsv 

Notes: 

**/

process socSummary {
    cache 'lenient'
    memory  '32 GB'

    publishDir path: "${params.outputDir}/", pattern: "soc_output_k*", mode: 'copy'

    input: 
        val kval from kval_range_copy02
        file soc_dir from soc_output

    output:
        file soc_dir into soc_summary

    script:    
    """
    # bash watch for errors
    set -ueo pipefail
    
    souporcell_summary.sh ${soc_dir}/clusters.tsv > ${soc_dir}/souporcell_summary.tsv
    """
}

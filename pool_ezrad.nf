#!/usr/bin/env nextflow

NXF_ANSI_LOG=false
//script parameters
params.reads	 =  "/home/drewx/Documents/pool-ezRAD/SimulationData"
params.reads_ext = "fastq.gz"
params.pattern 	 =  "*R{1,2}.fastq.gz"

//params.pattern 	=  "*_{1,2}.fastq"
params.output   =  "$PWD/dDocent"
params.hcpu	=  4
params.mcpu	=  4
params.barcode  = "barcodes"
params.min_cov  = 4

       
barcodes    	=  file(params.barcode)
reads_pattern 	=  params.reads + "/*" + params.pattern


Channel.fromFilePairs(reads_pattern)
       .ifEmpty{ exit 1, "params.reads empty no reads found" }
       .into{rad_tags; uniq_rad; raw_reads_bwa}     

//raw_reads.map{ it  -> [ it[1][0], it[1][1]] }
//         .set{ raw_reads_FastQC }




// process process_radtags{
//     echo true
//     publishDir  params.output + "/RadTags"
    
//     input:
//         set val(sample),  file(reads) from rad_tags
// 	file barcodes

//     output:
// 	file "process_radtags.log"  into rag_log
// 	file "sample*"  into demultiplexed_reads


// script:
//   (read1, read2) = reads

// """

//     process_radtags \
//         -1 ${read1} \
//         -2 ${read2} \
//         -b barcodes \
//         -e ecoRI \
//         --renz_2 mspI \
//         -r \
//         -i gzfastq

//     rm *rem*

// """


// }




process uniq_RAD{

    publishDir "$PWD/uniqRAD/RawReads"
    publishDir "$PWD/uniqRAD/seqcount_data", pattern: "*.uniq_seq"
    echo true
    
    input:
       set val(sample), file(reads) from uniq_rad

    output:
        file "${sample}.{forward,reverse}"
	file "${sample}.uniq_seq"  into (seqcounts1, seqcounts2, seqcounts3) 

script:
    (foward,reverse) = reads
    PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'


"""
   zcat ${foward} | mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | mawk '!/>/'   > ${sample}.forward
   zcat ${reverse} |mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | mawk '!/>/'  > ${sample}.reverse

   paste -d '-' ${sample}.forward  ${sample}.reverse | mawk '!/NNN/' |  sed 's/-/NNNNNNNNNN/' | perl -e  '$PERLT' > ${sample}.uniq_seq
   
"""
	 
}




// process sample_coverage{

//     echo true
//     publishDir "$PWD/uniqRAD/coveraga"
//     input:
//        file uniqseqs from  seqcounts1

//    output:
// 	file  "${uniqseqs.baseName}.coverage"


// shell:
// '''
//    for i in {2..20};
//    do 
//       echo \$i >> pfile
//    done

//    cat pfile | parallel \
//    --no-notice \
//    -j ${params.mcpu} \
//    "echo -n {}xxx && mawk -v x={} '\\$1 >= x' !{uniqseqs} | wc -l " \
//    | mawk  '{gsub("xxx","\t",\$0); print;}' \
//    | sort -g > !{uniqseqs.baseName}.coverage

// '''


// }




// process combined_coverage{
    
//     publishDir  "$PWD/uniqRAD/coverage"

//     input:
//         file uniqseqs from  seqcounts2.collect()
	
//     output:
// 	file  "combined.coverage"
	

// shell:
// '''
//    for i in {2..20};
//    do 
//    echo \$i >> pfile
//    done
  
//    cat pfile | parallel \
//    --no-notice \
//    -j ${params.mcpu} \
//    "echo -n {}xxx && mawk -v x={} '\\$1 >= x' !{uniqseqs} | wc -l " \
//    | mawk  '{gsub("xxx","\t",\$0); print;}' \
//    | sort -g > combined.coverage

// '''
// // //TO DO
// // //Simplify this shell pipe
// }




process Coverage_filtering{

    echo true
    publishDir "$PWD/uniqRAD/seqcount_data"
    input:
        file uniqseqs from  seqcounts3

    output:
        file "${uniqseqs.baseName}.uniqCperindv" 	

script:
   PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
   //while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";} 

"""    
   parallel --no-notice  \
   -j ${params.mcpu} \
   mawk -v x=${params.min_cov} \\''\$1 >= x'\\' ::: ${uniqseqs} \
   | cut -f2 \
   | perl -e '${PERLT}' > ${uniqseqs.baseName}.uniqCperindv
   wc -l ${uniqseqs.baseName}.uniqCperindv


"""


}

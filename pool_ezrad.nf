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
params.min_sample  = 4
params.clu_perc = 0.8
       
barcodes    	=  file(params.barcode)
reads_pattern 	=  params.reads + "/*" + params.pattern

Channel.fromFilePairs(reads_pattern)
       .ifEmpty{ exit 1, "params.reads empty no reads found" }
       .into{rad_tags; uniq_rad1; uniq_rad2; uniq_rad3; raw_reads_bwa}     

//raw_reads.map{ it  -> [ it[1][0], it[1][1]] }
//         .set{ raw_reads_FastQC }




process process_radtags{
    echo true
    publishDir  params.output + "/RadTags"
    
    input:
        set val(sample),  file(reads) from rad_tags
	file barcodes

    output:
	file "process_radtags.log"  into rag_log
	file "sample*"  into demultiplexed_reads


script:
  (read1, read2) = reads

"""

    process_radtags \
        -1 ${read1} \
        -2 ${read2} \
        -b barcodes \
        -e ecoRI \
        --renz_2 mspI \
        -r \
        -i gzfastq

    rm *rem*

"""


}



process uniq_RAD{

    publishDir params.output + "RawReads"
    publishDir params.output + "uniqRAD/seqcount_data", pattern: "*.uniq_seq"
    
    input:
       set val(sample), file(reads) from uniq_rad1

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




process sample_coverage{

    
    publishDir params.output + "/coveraga"
    
    input:
       file uniqseqs from  seqcounts1
   output:
	file  "${uniqseqs.baseName}.coverage"
	
shell:
'''
   for i in {2..20};
   do 
      echo \$i >> pfile
   done
   cat pfile | parallel \
   --no-notice \
   -j ${params.mcpu} \
   "echo -n {}xxx && mawk -v x={} '\\$1 >= x' !{uniqseqs} | wc -l " \
   | mawk  '{gsub("xxx","\t",\$0); print;}' \
   | sort -g > !{uniqseqs.baseName}.coverage
'''

}



process combined_coverage{
    
    publishDir  params.output + "/uniqRAD/coverage"
    
    input:
        file uniqseqs from  seqcounts2.collect()
    output:
       file  "combined.coverage"
	
shell:
'''

    for i in {2..20};
    do 
    echo \$i >> pfile
    done
    cat pfile | parallel \
    --no-notice \
    -j ${params.mcpu} \
    "echo -n {}xxx && mawk -v x={} '\\$1 >= x' !{uniqseqs} | wc -l " \
    | mawk  '{gsub("xxx","\t",\$0); print;}' \
    | sort -g > combined.coverage

'''
//TO DO
//Simplify this shell pipe

}



uniq_rad2.map{ it  -> [it[1][0], it[1][1]]}
                      .flatten()
                      .set{gunzipped}


process gunzip_fastq{
    
    echo true
    publishDir params.output + "/uniqRAD/seqcount_data"
    input:
        file fastq_gz from gunzipped
	
    output:
        file "${fastq_gz.baseName}*" into fastq_raw
	
script:
"""
  
  gunzip -vk --force  ${fastq_gz}
  
   
"""
}

       

process Coverage_filtering{
    
    echo true
    publishDir params.output + "/uniqRAD/seqcount_data"
	
    input:
        file reads from fastq_raw.collect()
	
    output:
        file "${sample}.uniqCperindv"
	
script:

     sample=reads[1].baseName
     PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'

	
"""

    parallel --no-notice  \
       -j ${params.mcpu} \
       mawk -v x=${params.min_cov} \\''\$1 >= x'\\' ::: ${reads} \
       | cut -f2 \
       | perl -e '${PERLT}' > ${sample}.uniqCperindv

    wc -l ${sample}.uniqCperindv
   
"""

	
}



process seq_filter{

    echo true
    publishDir params.output + "/uniqRAD/coverage"
    
    input:
       file uniqseqs from  uniq_rad3.collect()

   output:
      file  "*uniqCperindv"
      file  "count_k${params.min_cov}_c${params.min_sample}.seqs"
      file  "uniq_k${params.min_cov}_c${params.min_sample}.seqs" into uniq_filtrd


script:
   PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'


shell:
'''

   for ((i = 2; i <= 10; i++));
   do
     echo $i >> ufile
   done

   parallel --no-notice -j !{params.hcpu}  mawk -v x=!{params.min_cov} \\''$1 >= x'\\' :::  *  \
       | cut -f2 \
       | perl -e  '!{PERLT}' \
       > uniqCperindv

   cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\\$1 >= x' uniqCperindv | wc -l" \
       | mawk '{gsub("xxx","\t",$0); print;}' \
       | sort -g > uniqseq.peri.data

   wc -l uniqCperindv | tee - count.uniqCperindv

   mawk -v x=!{params.min_sample} '\$1 >= x' uniqCperindv > uniq_k!{params.min_cov}_c!{params.min_sample}.seqs

   wc -l uniq_k!{params.min_cov}_c!{params.min_sample}.seqs | tee  - count_k!{params.min_cov}_c!{params.min_sample}.seqs

'''


}



process get_contigs{

    echo true
    publishDir params.output + "contigs"
    
    input:
        file uniq_seqs from uniq_filtrd

    output:
	 file "totaluniqseq" into totaluniqseq
	 file "uniq.fasta"
	 file "uniq_Fwd.fasta" into uniq_FWD_reads

shell:
'''

    cut -f 2 !{uniq_seqs}  > totaluniqseq
    
    mawk '{c= c + 1; print ">Contig_" c "\\n" $1}' totaluniqseq > uniq.fasta

    sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f 1 > uniq_Fwd.fasta

'''

}



process cd_hit_FWD{


     publishDir params.output + "contigs"
     
     input:
         file uniq_FWD from uniq_FWD_reads

    output:
        file "${uniq_FWD.baseName}*" into cd_hit_clusters
 	//file "sort_contig.cluster_ids" into test


"""

    cd-hit-est \
       -i ${uniq_FWD} \
       -o ${uniq_FWD.baseName}_cdhit \
       -c ${params.clu_perc} \
       -T ${params.hcpu} \
       -M 0 \
       -g 1
      
"""

}



process  merge_contigs{

    echo true
    publishDir params.output + "clusters"
    
    input:
         file Fwd_clusters from cd_hit_clusters
	 file totaluniqseq

    output:
        file "sort_contig.cluster_ids" 
	file "contig_cluster.totaluniqseq"
        file "Rclusters" into rainbow_clusters

shell:
"""

   mawk '{if (\$1 ~ /Cl/) clus = clus + 1; else  print \$3 "\t" clus}'  ${Fwd_clusters}  \
      | sed 's/[>Contig_,...]//g' \
      | sort -g -k 1  > sort_contig.cluster_ids
   
   paste sort_contig.cluster_ids !{totaluniqseq}  > contig_cluster.totaluniqseq

   sort -k 2,2 -g contig_cluster.totaluniqseq \
      | sed -e 's/NNNNNNNNNN/\\t/g' > Rclusters 

"""

//sort_contig.cluster_ids is tsv file
//col1=contig_ID
//col2=Cluster_ID
//NB: Clustered become 1-indexed
//Rclusters 
//Read_ID	Cluster_ID	Forward_Read	Reverse_Read
}



process rainbow_div{

    echo true
    publishDir params.output + "clusters"
    input:
        file rainbow_clusters


    output:
         file "rainbow_div.out" into rainbow_div
  	    
"""
   
   rainbow \
       div \
       -i ${rainbow_clusters} \
       -f 0.5 \
       -K 10 \
       -o rainbow_div.out

"""
// Output File Format: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>[\t<pre_cluster_id:int>]
//   -k <int>    K_allele, min variants to create a new group [2]
//   -K <int>    K_allele, divide regardless of frequency when num of variants exceed this value [50]
//   The sets the number of reads to apply the freq/k filters	
}



process rainbow_merge{

    echo true
    publishDir params.output + "clusters"
    
    input:
        file   rainbow_div
	

    output:
         file "rainbow_assembly.out" into rainbow_assembly
	 
shell:
'''
   
    rainbow \
        merge \
        -i !{rainbow_div} \
        -r 2\
        -a \
        -o rainbow_assembly.rbasm

    cat rainbow_assembly.rbasm <(echo "E") |sed 's/[0-9]*:[0-9]*://g' |\
    mawk '{
    if (NR == 1) e=$2;
    else if (\$1 ~/E/ && lenp > len1){
        c=c+1; 
        print ">dDocent_Contig_" e "\\n" seq2 "NNNNNNNNNN" seq1; 
        seq1=0; 
        seq2=0;
        lenp=0;
        e=$2;
        fclus=0;
        len1=0;
        freqp=0;
        lenf=0
    }
    else if ($1 ~/E/ && lenp <= len1) {
        c=c+1; 
        print ">dDocent_Contig_" e "\\n" seq1; 
        seq1=0; 
        seq2=0;
        lenp=0;
        e=$2;
        fclus=0;
        len1=0;
        freqp=0;
        lenf=0
    }
    else if ($1 ~/C/) clus=$2;
    else if ($1 ~/L/) len=$2;
    else if ($1 ~/S/) seq=$2;
    else if ($1 ~/N/) freq=$2;
    else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
    else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
    else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {
        seq2 = seq; 
        lenp = len; 
        freqp=freq
    }
    }' > rainbow_contigs.fasta

    

'''
    
}

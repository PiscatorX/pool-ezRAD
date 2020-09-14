#!/usr/bin/env nextflow

NXF_ANSI_LOG=false

//script parameters
//params.reads	 = "/opt/DB_REF/Clinid.ezRAD/"
//params.reads   = "/home/drewx/Documents/pool-ezRAD/Test/"
params.reads     = "/home/drewx/Documents/pool-ezRAD/test1"
//params.reads_ext = "fastq.gz"
//params.pattern 	 =  "*R{1,2}_001.fastq.gz"
params.pattern 	 =  "*R{1,2}.fastq.gz"
params.output   =  "$PWD/SimRAD"
params.hcpu	=  4
params.mcpu	=  4
params.barcodes =  "/home/drewx/Documents/pool-ezRAD/test1/SimRAD.barcodes"
params.demultiplex = true
params.min_cov  = 4
params.min_sample  = 4
params.clu_perc = 0.8


reads_pattern 	=  params.reads + "/*" + params.pattern

rad_tags = Channel.fromFilePairs(reads_pattern)
                  .ifEmpty{ exit 1, "params.reads empty no reads found" }
		  

if (params.demultiplex){

barcodes = Channel.fromPath(params.barcodes)
                  .ifEmpty{ exit 1, "params.barcodes empty no reads found" }


process process_radtags{

    //echo true
    publishDir  params.output + "/RadTags", mode: 'copy'
    
    input:
        set val(sample),  file(reads) from rad_tags
	file barcodes

    output:
	file "process_radtags.log"
	file "*.{F,R}.fq.gz"  into demultiplexed_reads 


script:
  (read1, read2) = reads


"""
    cut -f 2 \
        ${barcodes} > \
        barcodes
    
    process_radtags \
        -1 ${read1} \
        -2 ${read2} \
        -b barcodes \
        -e ecoRI \
        --renz_2 mspI \
        -r \
        -i gzfastq

    rm *rem*

    Rename_for_dDocent.sh \
        SimRAD.barcodes

"""

}

    demultiplexed_reads = demultiplexed_reads.flatten().collate(2, false)

}
else{

   demultiplexed_reads = rad_tags.map{ it[1] }


}


process uniq_RAD{

    
    input:
        set foward, reverse from demultiplexed_reads

    output:
	file "${sample}.uniq_seq"  into (uniq_seqs1, uniq_seqs2, uniq_seqs3, uniq_seqs4) 


script:
    if(params.demultiplex){
       sample=foward.getName().replace(".F.fq.gz",'')
    }
    else{
       sample=foward.getName().replace("R1_001.fastq.gz",'')

    }
    
    PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'

"""

   zcat ${foward} | mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | mawk '!/>/'   > ${sample}.forward
   zcat ${reverse} |mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' | mawk '!/>/'  > ${sample}.reverse

   paste -d '-' ${sample}.forward  ${sample}.reverse | mawk '!/NNN/' |  sed 's/-/NNNNNNNNNN/' | perl -e  '$PERLT' > ${sample}.uniq_seq
  
 
"""


}



process sample_coverage{
    
    publishDir params.output + "/coverage", mode: 'copy'
    
    input:
        file uniqseqs from  uniq_seqs1
	
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
   -j !{params.mcpu} \
   "echo -n {}xxx && mawk -v x={} '\\$1 >= x' !{uniqseqs} | wc -l " \
   | mawk  '{gsub("xxx","\t",\$0); print;}' \
   | sort -g > !{uniqseqs.baseName}.coverage
   
'''

}



process combined_coverage{
    
    publishDir  params.output + "/coverage", mode: 'copy'
    
    input:
        file uniqseqs from  uniq_seqs2.collect()
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
    -j !{params.mcpu} \
    "echo -n {}xxx && mawk -v x={} '\\$1 >= x' !{uniqseqs} | wc -l " \
    | mawk  '{gsub("xxx","\t",\$0); print;}' \
    | sort -g > combined.coverage

'''
//TO DO
//Simplify this shell pipe


}



process coverage_filtering{
    
   
    publishDir params.output + "/uniqRAD/"
	
    input:
        file uniqseqs from uniq_seqs3.collect()
	
    output:
        file "combined.uniqCperindv"
	
script:
     PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'

	
"""

    parallel --no-notice  \
       -j ${params.mcpu} \
       mawk -v x=${params.min_cov} \\''\$1 >= x'\\' ::: ${uniqseqs} \
       | cut -f2 \
       | perl -e '${PERLT}' > combined.uniqCperindv

    wc -l combined.uniqCperindv
   
"""
	
}



process seq_filter{

    echo true
    publishDir params.output + "/coverage", pattern: "uniqseq.peri.data", mode: 'copy'
    
    
    input:
       file uniqseqs from  uniq_seqs4.collect()

   output:
      //file  "uniqCperindv"
      file  "uniqseq.peri.data"
      //file  "count_k${params.min_sample}_c${params.min_cov}.seqs"
      file  "uniq_k${params.min_sample}_c${params.min_cov}.seqs" into uniq_filtrd


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

   mawk -v x=!{params.min_sample} '\$1 >= x' uniqCperindv > uniq_k!{params.min_sample}_c!{params.min_cov}.seqs

   wc -l uniq_k!{params.min_sample}_c!{params.min_cov}.seqs | tee  - count_k!{params.min_sample}_c!{params.min_cov}.seqs

'''


}



process get_contigs{

    
    publishDir params.output + "/contigs"
    
    input:
        file uniq_seqs from uniq_filtrd

    output:
	 file "totaluniqseq" into totaluniqseq
	 file "uniq.fasta" into uniq_fasta
	 file "uniq_Fwd.fasta" into uniq_FWD_reads

shell:
'''
	
    cut -f 2 !{uniq_seqs}  > totaluniqseq
    
    txt2fasta.py totaluniqseq 

    sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f 1 > uniq_Fwd.fasta

'''

}



process cd_hit_FWD{


     //publishDir params.output + "/CD_Hit", mode: 'copy'
     
     input:
         file uniq_FWD from uniq_FWD_reads

    output:
        file "${uniq_FWD.baseName}_cdhit"
	file "${uniq_FWD.baseName}_cdhit.clstr" into cd_hit_clusters
 	

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

   
    //publishDir params.output + "/clusters", mode: 'copy'
    
    input:
         file Fwd_clusters from cd_hit_clusters
	 file totaluniqseq
	 file uniq_fasta
	 
    output:
        file "sort_contig.cluster_ids" 
        file "Rclusters" into rainbow_clusters

shell:
"""

   mawk '{if (\$1 ~ /Cl/) clus = clus + 1; else  print \$3 "\t" clus}'  ${Fwd_clusters}  \
      | sed 's/[>Contig_,...]//g' \
      | sort -g -k 1  > sort_contig.cluster_ids
   
  contig2clstr.py \
      -s  ${uniq_fasta} \
      -c sort_contig.cluster_ids | sed -e 's/NNNNNNNNNN/\\t/g' > Rclusters 

"""

//sort_contig.cluster_ids is tsv file
//col1=contig_ID
//col2=Cluster_ID
//NB: Clustered become 1-indexed
//Rclusters 
//Read_ID	Cluster_ID	Forward_Read	Reverse_Read
}



process rainbow_div{

    
    publishDir params.output + "/clusters", mode: 'copy'
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

    
    publishDir params.output + "/contigs", mode: 'copy'
    
    input:
        file   rainbow_div
	

    output:
         file "rainbow_assembly.rbasm" into rainbow_assembly
	 file "rainbow_contigs.fasta"  into rainbow_contigs	 
shell:
'''
   
    rainbow \
        merge \
        -i !{rainbow_div} \
        -r 1 \
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



process cd_hit_contigs{


    publishDir params.output + "/ezRAD_contigs", mode: 'copy'
    
    input:
        file merged_contigs  from rainbow_contigs
	
    output:
        file "${merged_contigs.baseName}*" into filterd_contigs 
	    
script:
"""

     cd-hit-est \
	 -i ${merged_contigs} \
	 -o ${merged_contigs.baseName}_cdhit \
	 -c ${params.clu_perc} \
	 -T ${params.hcpu} \
	 -M 0 

"""

}

//TO DO
//write a script to describe the sequence distribution
//plot sequence distribution

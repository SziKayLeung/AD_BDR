################################################################################################
## Define functions for All_Demultiplex_Final_Functions.sh
# 1) run_CCS <input_ccs_bam> <prefix_output_name> 
# 2) run_LIMA <sample> <"no_multiplex"/"multiplex">
# 3) run_REFINE <sample> 
# 4) run_targeted_REFINE <sample> <input_config_file>
# 5) merging_at_refine <output_name> <samples.....>
# 6) run_CLUSTER $Sample 
# 7) run_map_cupcakecollapse <sample_prefix_input/output_name> 
# 8) demux_targeted <sample>
# 9) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
# 10) rnaseq_merge_fastq <sample>
# 11) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# 12) run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna>
# 13) TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna>
# 14) TAMA_sqanti_filter <sample> <mode=basic/full/nokallisto/lncrna>
# 15) parse_stats_per_sample <sample>
# 16) on_target_rate <input_probe_bed.file> <input_polished.hq.fasta> <input_mappped.fastq.sam> <output_file_name>
# 17) run_target_rate 
# 18) counts_subset_4tappas <input_class> <output_class> <type_genes>
# 19) TAMA_tappas_input <sample> <mode=basic/full/nokallisto/lncrna>
################################################################################################

#************************************* DEFINE VARIABLES
module load Miniconda2/4.3.21

# Listing versions
source activate isoseq3
ccs --version #ccs 5.0.0 (commit v5.0.0)
lima --version #lima 2.0.0 (commit v2.0.0)
isoseq3 --version #isoseq3 3.4.0 (commit v3.4.0)

source activate sqanti2_py3
echo "Minimap2 version:" $(minimap2 --version) # echo version


################################################################################################
#*************************************  Isoseq3 [Function 1,2,3,4,5,6]
# run_CCS <input_ccs_bam> <prefix_output_name> 
run_CCS(){
  source activate isoseq3
  
  mkdir -p $WKD_ROOT/1_ccs; cd $WKD_ROOT/1_ccs
  echo "Processing Sample $2 from Functions script"
  echo "Processing $1"
  echo "Checking if $2.ccs.bam exists"

    if [ -f $2.ccs.bam ]; then
      echo "$2.ccs.bam file already exists; CCS no need to be processed on Sample $2"
    else
      echo "$2.ccs.bam file does not exist"
      echo "Processing CCS for sample $2"
      # ccs <input.subreads.bam> <output.ccs.bam>
      time ccs $1 $2.ccs.bam --minPasses=1 --min-rq 0.9 --reportFile $2_ccs_report.txt
      echo "CCS for Sample $2 successful"
      ls *$2*
    fi
  source deactivate
}

# run_LIMA <sample> <"no_multiplex"/"multiplex">
# removed --no pbi as this is needed for downstream polishing
run_LIMA(){
  source activate isoseq3

  echo "Processing $1 file for demultiplexing"
  
  if [ -f $WKD_ROOT/2_lima/$1.fl.json ]; then
    echo "$1.fl.bam file already exists; LIMA no need to be processed on Sample $1"
  elif [ $2 = "multiplex" ]; then
    echo "Multiplex: use Targeted_FASTA"
    #lima <input.ccs.merged.consensusreadset.xml> <input.primerfasta> <output.fl.bam>
    mkdir -p $WKD_ROOT/2_lima; cd $WKD_ROOT/2_lima
    time lima $WKD_ROOT/1_ccs/$1.ccs.bam $TARGETED_FASTA $1.fl.bam --isoseq --dump-clips --dump-removed --peek-guess
    echo "lima $1 successful"
    ls $1.fl*
  else
    echo "No-Multiplex: use FASTA"
    mkdir -p $WKD_ROOT/2b_lima_batches; cd $WKD_ROOT/2b_lima_batches
    time lima $WKD_ROOT/1_ccs/$1.ccs.bam $FASTA $1.fl.bam --isoseq --dump-clips --dump-removed
    echo "lima $1 successful"
    ls $1.fl*
  fi
  source deactivate
}

# run_REFINE <sample>
run_REFINE(){
  source activate isoseq3
  
  mkdir -p $WKD_ROOT/3b_refine_batches; cd $WKD_ROOT/3b_refine_batches

  echo "Processing $1 file for refine"
  if [ -f $1.flnc.bam ]; then
    echo "$1.flnc bam file already exists; Refine no need to be processed on Sample $1"
  else
    #refine --require-polya <input.lima.consensusreadset.xml> <input.primer.fasta> <output.flnc.bam>
    time isoseq3 refine $WKD_ROOT/2b_lima_batches/$1.fl.primer_5p--primer_3p.bam $FASTA $1.flnc.bam --require-polya
    echo "refine $1 successful"
    ls $1.flnc*
  fi

  source deactivate
}

# run_targeted_REFINE
run_targeted_REFINE(){
  source activate isoseq3
  
  mkdir -p $WKD_ROOT/3_refine; cd $WKD_ROOT/3_refine
  
  # barocode.configs file contain 3 columns, tab separated
  # column 1 = sample name output
  # column 2 = barcode number
  # column 3 = batch name 
  grep "^[^#;]" $BARCODE_CONFIG | while read p; do
    # tab separate cut and paste to variable 
    sample=$(echo "$p" | cut -d $'\t' -f 1)
    barcode=$(echo "$p" | cut -d $'\t' -f 2)
    batch=$(echo "$p" | cut -d $'\t' -f 3)
    
    # append the batch and barcode name to extract corresponding file from lima 
    lima_file=$batch."fl.primer_5p--BC100"$barcode"_3p.bam"
    echo "Working with sample $sample = $lima_file"
    
    # isoseq refine for individual files
    time isoseq3 refine $WKD_ROOT/2_lima/$lima_file $TARGETED_FASTA $sample.flnc.bam --require-polya
    echo "refine $sample successful"
  done

  source deactivate
}

# merging_at_refine <output_name> <samples.....>
# aim: merging bam files from refine onwards (similar to run_isoseq3_2_1_merge, but no need to rerun from ccs)
# <output_name> = output name for merged files
# samples.... = list of the sample names
merging_at_refine(){

  source activate isoseq3

  ###********************* Merging at REFINE
  # Define variable "Merge_Samples" as a list of all samples, in order to find the specified flnc.bam (for dataset create ConsensusReadSet)
  # Define variable "all_flnc_bams" for merged path-directory of all flnc samples (for dataset create ConsensusReadSet)
  Merge_Samples=$(echo "${@:2}")

  echo "Merging flnc of samples $Merge_Samples"
  all_flnc_bams=$(
      for i in ${Merge_Samples[@]}; do
          flnc_bam_name=$(find $WKD_ROOT/3_refine -name "*.flnc.bam" -exec basename \{} \; | grep ^$i )
          flnc_bam=$(find $WKD_ROOT/3_refine -name "*.flnc.bam" | grep "$flnc_bam_name" )
          echo $flnc_bam
      done
  )

  cd $WKD_ROOT/5_merged_cluster
  printf '%s\n' "${all_flnc_bams[@]}" > $1.flnc.fofn
  cat $1.flnc.fofn

  ###*********************

  isoseq3 cluster $1.flnc.fofn $1.clustered.bam --verbose --use-qvs
  gunzip *.gz*

  source deactivate
}

# run_CLUSTER <sample>
run_CLUSTER(){
  source activate isoseq3

  mkdir -p $WKD_ROOT/4b_cluster_batches; cd $WKD_ROOT/4b_cluster_batches
  echo "Processing $1 file for cluster"
  if [ -f $1.unpolished.bam ]; then
    echo "$1.unpolished.bam file already exists; Cluster no need to be processed on Sample $1"
  else
    # cluster <input.flnc.bam> <output.unpolished.bam>
    time isoseq3 cluster $WKD_ROOT/3b_refine_batches/$1.flnc.bam $1.clustered.bam --verbose --use-qvs 2> $1.cluster.log
    echo "cluster $1 successful"
    ls $1.clustered*
  fi

  source deactivate
}


################################################################################################
#************************************* Post_Isoseq3 (Minimap2, Cupcake, Demultiplex) [Function 7,8]
# run_map_cupcakecollapse <output_name> 
# Prerequisite: mm10 cage peak
run_map_cupcakecollapse(){


    source activate sqanti2_py3

    # convert between fasta and fastq for downstream process
    echo "fasta2fastq conversion"
    python $SEQUENCE/fa2fq.py $WKD_ROOT/5_merged_cluster/$1.clustered.hq.fasta

    samtools --version # echo version
    echo "Minimap2 version:" $(minimap2 --version) # echo version

    echo "Processing Sample $1 for Minimap2 and sort"
    mkdir -p $WKD_ROOT/6_minimap; cd $WKD_ROOT/6_minimap 
    minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $GENOME_FASTA $WKD_ROOT/5_merged_cluster/$1.clustered.hq.fastq > $1.sam 2> $1.map.log
    samtools sort -O SAM $1.sam > $1.sorted.sam

    # Alignment stats
    mkdir PAF; cd PAF
    source activate nanopore
    htsbox samview -pS $WKD_ROOT/6_minimap/$1.sorted.sam > $1.paf
    awk -F'\t' '{if ($6="*") {print $0}}' $1.paf > $1.allread.paf # all reads
    awk -F'\t' '{if ($6=="*") {print $0}}' $1.paf > $1.notread.paf
    awk -F'\t' '{if ($6!="*") {print $0}}' $1.paf > $1.filtered.paf
    awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $1.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $1"_reads_with_alignment_statistics.txt"
    echo "Number of mapped transcripts to genome:"
    wc -l $1.filtered.paf
    echo "Number of ummapped transcripts to genome:"
    wc -l $1.notread.paf

    echo "Processing Sample $1 for TOFU, with coverage 85% and identity 95%"
    source activate cupcake
    cd $WKD_ROOT/7_tofu
    collapse_isoforms_by_sam.py -c 0.85 -i 0.95 --input $WKD_ROOT/5_merged_cluster/$1.clustered.hq.fastq --fq -s $WKD_ROOT/6_minimap/$1.sorted.sam --dun-merge-5-shorter -o $1 &>> $1.collapse.log
    get_abundance_post_collapse.py $1.collapsed $WKD_ROOT/5_merged_cluster/$1.clustered.cluster_report.csv &>> $1.abundance.log

    source activate sqanti2_py3
    # convert rep.fq to rep.fa for SQANTI2 input
    seqtk seq -a $1.collapsed.rep.fq > $1.collapsed.rep.fa
    echo "Processing Sample $1 for TOFU successful"
    source deactivate
}

# run_minimap2 <output_name> <input_fasta> 
run_minimap2(){
    source activate sqanti2_py3
    cd $4
    echo "Processing Sample $1 for Minimap2 Mapping and sort"
  
    minimap2 -t 32 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $GENOME_FASTA $2 > $1.sam 2> $1.sam.log
    sort -k 3,3 -k 4,4n $1.sam > $1.sorted.sam
    
    source deactivate
}


# demux_targeted <output_name>
# read in read.stat file from cupcake collapse output and count abundance of each sample based on CCS ID obtained from each refine report of the sample
demux_targeted(){
    source activate sqanti2_py3

    # script.R <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
    Rscript $DEMUXFUNCTIONS/Demultiplex_Cupcake.R $WKD_ROOT/3_refine \
      $WKD_ROOT/5_merged_cluster/$1.clustered.cluster_report.csv \
      $WKD_ROOT/7_tofu/$1.collapsed.read_stat.txt \
      $WKD_ROOT/7_tofu/$1.Demultiplexed_Abundance.txt
}

################################################################################################
#************************************* RNAseq [Function 9, 10]
# run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
# Aim: Individually align RNASeq filtered fastq files to STAR-indexed genome using STAR, and subsequently input SJ.out.bed file into SQANTI2
run_star(){
    source activate sqanti2
    echo "STAR Version"
    STAR --version
    # samtools version: Version: 1.9 (using htslib 1.9)

    # extract official sample name from all_filtered directory
    # extract only files with "fastq.filtered", beginning with sample name, and R1/R2

    # cd to $J20/Tg4510_input_directory
    cd $2
    echo "Finding reads for Sample $1 for STAR"
    # find reverse and forward file, trimming "./" therefore only printing file name
    F=$(find . | grep "fastq" | grep $1 | grep "R1" | sed 's|^./||' )
    R=$(find . | grep "fastq" | grep $1 | grep "R2" | sed 's|^./||' )
    # save path directory of files as variable for later mapping
	  F_File=$(realpath $F)
    R_File=$(realpath $R)
    echo "Processing Forward Reads: $F_File"
    echo "Processing Reverse Reads: $R_File"
	
	  STAR_reference_input_directory=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/hg38
	
	  cd $3
    mkdir $1 
    cd $1
    
    # Parameters according to https://github.com/jennylsmith/Isoseq3_workflow/blob/master/shell_scripts/8__STAR_Junctions_ShortReads.sh
    # 2-pass mode alignment with chimeric read detection
    # at least 25 bp of one read of a pair maps to a different loci than the rest of the read pair
    # require 20 pb on either side of chimeric junction
    # include chimeric reads in the output BAM
    # don't include chimeras that map to reference genome contig with Ns
    # --outSAMtype BAM SortedByCoordinate, output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.
    # note STAR indexed with genocde vM22 primary assembly annotation gtf therefore no need to include into command line, otherwise create additional folders
    STAR --runMode alignReads --runThreadN 32 --genomeDir $STAR_reference_input_directory \
    --readFilesIn $F_File $R_File \
	  --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --chimSegmentMin 25 \
  	--chimJunctionOverhangMin 20 \
  	--chimOutType WithinBAM \
  	--chimFilter banGenomicN \
  	--chimOutJunctionFormat 1 \
  	--twopassMode Basic \
  	--twopass1readsN -1 &>> $1_star.log #use all reads
      #--sjdbGTFfile $4/gencode.vM22.primary_assembly.annotation.gtf \
  
    # to remove duplicates between samples
    picard MarkDuplicates INPUT=$3/$1".sorted.bam" OUTPUT=$1".sorted.nodup.bam" METRICS_FILE=$1".dup.txt" VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> $1.PicardDup.log
  
    # rename output files
    mv Aligned.sortedByCoord.out.bam $1.sorted.bam
    mv Log.final.out $1.Log.final.out
    mv Log.out $1.Log.out
    mv Log.progress.out $1.Log.progress.out
    mv SJ.out.tab $1.SJ.out.bed
    mv $1.SJ.out.bed ../
}

# rnaseq_merge_fastq <output_name>
# Aim: Merge forward and reverse of all RNASeq filtered fastq, for further downstream alignment to IsoSeq Tofu output using Kallisto
rnaseq_merge_fastq(){
    
    cd $RNASEQ_FILTERED_DIR
    R1_READS=($(
        for i in ${RNASEQ_SAMPLES_NAMES[@]}; do
            F_name=$(find . | grep "fastq" | grep $i | grep "R1" | sed 's|^./||' )
            F_File=$(realpath $F_name)
            echo $F_File
        done
    ))

    R2_READS=($(
        for i in ${RNASEQ_SAMPLES_NAMES[@]}; do
            R_name=$(find . | grep "fastq" | grep $i | grep "R2" | sed 's|^./||' )
            R_File=$(realpath $R_name)
            echo $R_File
        done
    ))

    R1_READS_MERGE=$(echo "${R1_READS[@]}")
    R2_READS_MERGE=$(echo "${R2_READS[@]}")

    echo "Processing R1 READS"
    echo $R1_READS_MERGE
    echo "Processing R2 READS"
    echo $R2_READS_MERGE

    cd $RNASEQ_MAPPED_DIR
    zcat $R1_READS_MERGE > $1"_R1.fq"
    zcat $R2_READS_MERGE > $1"_R2.fq"
}


################################################################################################
#************************************* RNASeq & IsoSeq [Function 11]
# run_kallisto <sample_prefix_output_name> 
# Aim: Align merged RNASeq fastq files to IsoSeq Tofu output fasta using Kallisto
# Prerequisite:
    # run mouse_merge_fastq/any equivalent merging of all RNASeq fastq files (note sample_prefix_output name is the same)
# Output: <output_prefix>.mod.abundance.tsv for input into Sqanti_qc.py (not modified script to take in space delimited rather than tab file)
run_kallisto(){
    source activate sqanti2

    echo "Processing Kallisto for $1"
    mkdir -p $WKD_ROOT/8_kallisto; cd $WKD_ROOT/8_kallisto
    kallisto version
    time kallisto index -i $1_Kallisto.idx $WKD_ROOT/7_tofu/$1.collapsed.rep.fa 2> $1_Kallisto.index.log
    time kallisto quant -i $WKD_ROOT/8_kallisto/$1_Kallisto.idx --rf-stranded $RNASEQ_MAPPED_DIR/$1"_R1.fq" $RNASEQ_MAPPED_DIR/$1"_R2.fq" -o $WKD_ROOT/8_kallisto 2> $1_Kallisto.quant.log
    mv abundance.tsv $1.abundance.tsv

    # problem: retained co-ordinates, which does not input well into SQANTI2
    echo "Kallisto original $1.abundance.tsv"
    head $1.abundance.tsv
    # solution: retain the PB.ID
    while read line ; do
	    first=$( echo "$line" |cut -d\| -f1 ) # each read line, and remove the first part i.e. PacBio ID
	    rest=$( echo "$line" | cut -d$'\t' -f2-5 ) #save the remaining columns
	    echo $first $rest # concatenate
    done < $WKD_ROOT/8_kallisto/$1.abundance.tsv > $WKD_ROOT/8_kallisto/$1.temp.abundance.tsv

    header=$(head -n 1 $WKD_ROOT/8_kallisto/$1.abundance.tsv)
	  sed -i '1d' $WKD_ROOT/8_kallisto/$1.temp.abundance.tsv # remove header of temp.file to be replaced
    echo $header > foo
	  cat foo $WKD_ROOT/8_kallisto/$1.temp.abundance.tsv > $WKD_ROOT/8_kallisto/$1.mod.abundance.tsv

    echo "Kallisto $1.mod.abundance.tsv"
    head $WKD_ROOT/8_kallisto/$1.mod.abundance.tsv
    rm $1.temp.abundance.tsv
	  rm foo

    source deactivate
}



################################################################################################
#************************************* SQANTI2 [Function 12]

# run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna>
run_sqanti3(){

    source activate sqanti2_py3
    
    # variable 
    sample=$1.collapsed
    gtf=$WKD_ROOT/7_tofu/$1.collapsed.gff
    abundance=$WKD_ROOT/7_tofu/$1.Demultiplexed_Abundance.txt
    
    
    # create directory
    mkdir -p $WKD_ROOT/9_sqanti3 
    mkdir -p $WKD_ROOT/9_sqanti3/$2; cd $WKD_ROOT/9_sqanti3/$2
    
    # copy STAR output SJ.out.bed files
    SJ_OUT_BED=($(
        for rnaseq in ${RNASEQ_SAMPLES_NAMES[@]}; do
            name=$(find $RNASEQ_MAPPED_DIR -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
            File=$(find $RNASEQ_MAPPED_DIR -name "$name")
            echo $File
        done
        ))
    for file in ${SJ_OUT_BED[@]}; do cp $file $WKD_ROOT/9_sqanti3/with_junc ;done
    
    # prepare kallisto file 
    if [ $2 == "full" ] || [ $2 == "lncrna" ]; then 
      # create tab separated file from kallisto
      # Rscript script.R <input.file_fromkallisto> <output.file_intosqanti>
      Rscript $KALLSTOINPUT $WKD_ROOT/8_kallisto $WKD_ROOT/9_sqanti3/$2/$sample".mod2.abundance.tsv"
      kallisto_expfile=$WKD_ROOT/9_sqanti3/$2/$sample".mod2.abundance.tsv"
      echo "Using $kallisto_expfile"
    fi
    
    # sqanti qc
    echo "Processing Sample $sample for SQANTI3 QC"
    python $SQANTI3_DIR/sqanti3_qc.py -v
    echo $GENOME_GTF
    echo $GENOME_FASTA
    
    if [ $2 == "basic" ]; then
      echo "Processing basic commands"
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
      --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report pdf &> $sample.sqanti.qc.log

    elif [ $2 == "full" ]; then
      echo "Processing full commands"
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
      -c ./"*SJ.out.bed" --expression $kallisto_expfile --fl_count $abundance \
      --genename --isoAnnotLite --gff3 $GFF3 --report pdf &> $sample.sqanti.qc.log

    elif [ $2 == "nokallisto" ]; then 
      echo "Processing basic commands with RNA-Seq bed files but not kallisto"
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
      -c "./*SJ.out.bed" \
      --genename --isoAnnotLite --gff3 $GFF3 --report pdf &> $sample.sqanti.qc.log

    elif [ $2 == "lncrna" ]; then
      echo "Processing with $LNCRNA_GTF for genome annotation "
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $LNCRNA_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
      -c "./*SJ.out.bed" --skipORF --fl_count $abundance \
      --genename --isoAnnotLite --gff3 $GFF3 --report pdf &> $sample.sqanti.qc.log

    else
      echo "2nd argument required"
    fi

    echo "Processing Sample $sample for SQANTI2 filter"
    python $SQANTI3_DIR/sqanti3_RulesFilter.py $sample"_classification.txt" $sample"_corrected.fasta" $sample"_corrected.gtf" -a 0.6 -c 3 &> $1.sqanti.filter.log
    
    if [ $2 != "basic" ]; then 
      # remove temp SJ.out bed files
      rm *SJ.out.bed
    fi
    
    source deactivate
}



################################################################################################
#************************************* TAMA [Function 13,14]
# TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna>
# remove short fragments from post tofu
# Prerequisite: Require TAMA_prepare.R to change column 4 of bed12 to gene_name: transcript_name for correct file TAMA format
#TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_remove_fragments(){
  
    mkdir -p $WKD_ROOT/9b_filter_cont $WKD_ROOT/9b_filter_cont/tama
    mkdir -p $WKD_ROOT/9b_filter_cont/tama/$2; cd $WKD_ROOT/9b_filter_cont/tama/$2
    
    # variables 
    sample=$1.collapsed
    io_dir=$WKD_ROOT/9b_filter_cont/tama/$2

    source activate sqanti2

    # convert gtf to bed12
	  gtfToGenePred $WKD_ROOT/9_sqanti3/$2/$sample"_classification.filtered_lite.gtf" $1.genepred
	  genePredToBed $1.genepred $1.bed12
    awk -F'\t' '{print $1,$2,$3,$4,"40",$6,$7,$8,"255,0,0",$10,$11,$12}' $1.bed12| sed s/" "/"\t"/g|sed s/",\t"/"\t"/g|sed s/",$"/""/g > Tama_$1.bed12

    # Rscript script.R <name>.bed12 <input_dir>
	  Rscript $TAMAMERGE Tama_$1 $io_dir
	  python $TAMA_DIR/tama_remove_fragment_models.py -f Tama_$1_mod.bed12 -o $io_dir

	  rm *Tama*

    echo "Number of isoforms filtered by TAMA:"
    wc -l $1"_discarded.txt"

	  source deactivate
}

# TAMA_sqanti_filter <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_sqanti_filter(){
  
  # variables
  sqname=$1.collapsed_classification.filtered_lite
  sq_dir=$WKD_ROOT/9_sqanti3/$2
  io_dir=$WKD_ROOT/9b_filter_cont/tama/$2
  
  cd $WKD_ROOT/9b_filter_cont/tama/$2
  
  source activate sqanti2_py3
  # Rscript .R <path/tama_filtered_output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_junc_txt> <output_prefix_name> <output_dir>
  Rscript $TAMASUBSET $io_dir/$1.bed $sq_dir $sqname"_classification.txt" $sqname".gtf" $sqname"_junctions.txt" $1 $io_dir
  
  # extract fasta sequence based on the pbid (https://www.biostars.org/p/319099/)
  # script.py <path/sqanti_filtered.fasta> <path/retained_pbid_tama.txt> <path/output.fasta>
  awk '{ print $4 }' $1| cut -d ";" -f 2  > $io_dir/tama_retained_pbid.txt
  python $TAMASUBSETFASTA $sq_dir/$sqname".fasta" $io_dir/tama_retained_pbid.txt $io_dir/$1"_sqantifiltered_tamafiltered_classification.fasta"

  Rscript $SQ_Report $1"_sqantitamafiltered.classification.txt" $1"_sqantitamafiltered.junction.txt"
}

# remove_3ISM <sample> <mode=basic/full/nokallisto/lncrna>
remove_3ISM(){
  
  mkdir -p $WKD_ROOT/9b_filter_cont $WKD_ROOT/9b_filter_cont/no3ISM
  mkdir -p $WKD_ROOT/9b_filter_cont/no3ISM/$2; cd $WKD_ROOT/9b_filter_cont/no3ISM/$2
  
  # variables 
  sq_dir=$WKD_ROOT/9_sqanti3/$2
  
  # Rscript script.R <input.classfile> <input.gtf> <input.junc> <output.dir> <prefix>
  Rscript $ISMREMOVE $WKD_ROOT/9_sqanti3/$2/$1.collapsed_classification.filtered_lite_classification.txt \
    $sq_dir/$1.collapsed_classification.filtered_lite.gtf \
    $sq_dir/$1.collapsed_classification.filtered_lite_junctions.txt $WKD_ROOT/9b_filter_cont/no3ISM/$2 $1

}


################################################################################################
#************************************* QC [Function 15, 16, 17]
# parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
parse_stats_per_sample(){

    source activate sqanti2_py3
    
    python $CCSQC $WKD_ROOT/1_ccs $WKD_ROOT/1_ccs $1
    python $LIMAQC $WKD_ROOT/2_lima $WKD_ROOT/2_lima $1
    
    source deactivate
}

# on_target_rate <input_probe_bed.file> <input_polished.hq.fasta> <input_mappped.fastq.sam> <output_file_name>
on_target_rate(){
    source activate sqanti2_py3
 
    echo "Processing input fasta: $2"
    echo "Processing input mapped fasta: $3"
    echo "Writing output: $4"
    python $CUPCAKE/targeted/calc_probe_hit_from_sam.py $1 $2 $3 --start_base 0 --end_base 0 --o $4
    
    source deactivate
}

# run_target_rate 
run_target_rate(){
  
  mkdir -p $WKD_ROOT/6_minimap/Individual 
  mkdir -p $WKD_ROOT/6b_target_rate
  
  for i in ${ALL_SAMPLES_NAMES[@]}; do 
    
    # run_minimap2 <output_name> <input_fasta> 
    gunzip $WKD_ROOT/4_cluster/$i.clustered.hq.fasta.gz
    cd $WKD_ROOT/6_minimap/Individual; run_minimap2 $i $WKD_ROOT/4_cluster/$i.clustered.hq.fasta
    
    # on_target_rate
    cd $WKD_ROOT/6b_target_rate; on_target_rate $PROBES $WKD_ROOT/4_cluster/$i.clustered.hq.fasta \
    $WKD_ROOT/6_minimap/Individual/$i.clustered.hq.fastq.sam \
    $WKD_ROOT/6b_target_rate/$i".fasta.sam.probe_hit.txt"

  done
}


################################################################################################
# run_cpat <input_fasta> <output_name>
run_cpat(){
  mkdir -p $WKD_ROOT/10_characterise
  mkdir -p $WKD_ROOT/10_characterise/CPAT; cd $WKD_ROOT/10_characterise/CPAT
  
  source activate sqanti2_py3
  cpat.py -x $HEXAMER -d $LOGITMODEL -g $1 --min-orf=50 --top-orf=50 -o $2"_cpat" 2> $2"_cpat.e"
}

# extract_best_orf <sample>
# extract the best ORF from CPAT for further analysis of ORF predictions for predicted NMD
extract_best_orf(){
  
  # variables
  io_dir=$WKD_ROOT/10_characterise/CPAT
  
  cd $io_dir
  python $BESTORF --fa $io_dir/$1"_cpat.ORF_seqs.fa" --orf $io_dir/$1"_cpat.ORF_prob.best.tsv" --o_name $1"_bestORF" --o_dir $io_dir &> orfextract.log
}


# convert_gtf_bed12 <input_gtf> 
convert_gtf_bed12(){
    
    # variables 
    output_dir="$(dirname $1)" 
    sample=${1%.gtf} # removes .gtf

  	source activate sqanti2_py3
    cd $output_dir

    gtfToGenePred $1 $sample.genePred
    genePredToBed $sample.genePred > $sample.bed12
    sort -k1,1 -k2,2n $sample.bed12 > $sample"_sorted.bed12"
    rm $sample.genePred $sample.bed12
    # Rscript script.R <input.classfile> <input_output_dir of bed file> <prefix>
    #Rscript $GENERAL/annotate_uscs_tracks.R $SQANTI $WKD $sample
    #bedToBigBed -extraIndex=name $sample $sample"_Modified.bed12" $REFERENCE/mm10.chrom.filtered.sizes $sample.bb

    #source deactivate
    #./bedToBigBed -as=bedExample2.as -type=bed9+3 -extraIndex=name $sample"_Modified.bed12" $REFERENCE/mm10.chrom.filtered.sizes $sample.bb
}

# colour_by_abundance <sample> <input_gtf>
colour_by_abundance(){
  
  convert_gtf_bed12 $2
  bed12=${2%.gtf}_sorted.bed12
  sample="$(basename $2)" 
  output_bed12=${sample%.gtf}_sorted_coloured.bed12

  # Merge FL across samples for abundance file
  # Rscript script.R <input.abundance> <output_dir> <prefix>
  Rscript $MERGEFL $WKD_ROOT/7_tofu/$1.Demultiplexed_Abundance.txt $WKD_ROOT/7_tofu $1
  
  python $COLOURTRANS \
  --bed $bed12 \
  --cpat $WKD_ROOT/10_characterise/CPAT/$1"_cpat.ORF_prob.best.tsv" \
  --noORF $WKD_ROOT/10_characterise/CPAT/$1"_cpat.no_ORF.txt" \
  --a $WKD_ROOT/7_tofu/$1"_Merged_Abundance.csv" \
  --o $WKD_ROOT/10_characterise/$output_bed12
  
}

# subset_gene_reference 
subset_gene_reference(){
  
  mkdir -p $WKD_ROOT/10_characterise/TargetGenesRef
  
  python $TGENEPREP --r=$GENOME_GTF --glist ${TGENES[@]} --o $WKD_ROOT/10_characterise/TargetGenesRef
  #TARGETS_GENCODE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Utilities/Targets_Gencode
  #for g in ${TargetGene[@]}; do grep_reference_gene $g $TARGETS_GENCODE; done
}


# full_characterisation <gene>
full_characterisation(){
  
  mkdir -p $WKD_ROOT/10_characterise/TargetGenes
  
  # variables 
  ref_dir=$WKD_ROOT/10_characterise/TargetGenesRef/
  input_bed=$WKD_ROOT/9_sqanti3/basic/$1.collapsed_corrected_sorted.bed12
  input_gtf=$WKD_ROOT/9_sqanti3/basic/$1.collapsed_classification.filtered_lite.gtf
  noISM_path=$WKD_ROOT/9b_filter_cont/no3ISM/basic/$1"_ISMrem.classification.txt"
  ORF_dir=$WKD_ROOT/10_characterise/CPAT/$1"_cpat.ORF_prob.best.tsv"  
  
  source activate sqanti2_py3 
  python $FICLE --gene=$2 --ref=$ref_dir \
  --i_bed=$input_bed \
  --i_gtf=$input_gtf \
  --noISM=$noISM_path \
  --orf=$ORF_dir \
  --o_dir=$WKD_ROOT/10_characterise/TargetGenes/ #&> $WKD_ROOT/10_characterise/TargetGenes/$2"_characterise.log"
  
}
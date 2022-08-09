################################################################################################
#************************************* Prepare Input files for tappAS [Function 18, 19]
subset_targets(){
  # Rscript script.R <output_dir/output_file>
  Rscript $SUBSETTARGET $ISO_WKD_ROOT/$SQNAME"_classification.txt" $SPECIES $1/$NAME"_TargetTrans.txt"
}

# counts_subset_4tappas <input_class> <output_class> <type_genes>
# regenerate expression matrix for tappas from long reads FL read counts
counts_subset_4tappas(){
  
  # variables
  input_class=$1
  output_class=$2
  type_genes=$3
  
  # Rscript script.R <input.classfile> <output.classfile> <type>
  # type = AD or Non-AD
  source activate sqanti2_py3
  Rscript $SUBSETEXP $input_class $output_class $type_genes
}

# TAMA_tappas_input <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_tappas_input(){
  
  # variable
  sample=$1"_sqantitamafiltered
  
  # isoannolite_generate <input_dir> <input_gtf> <input_class> <input_junc> <species> <output_dir> <output_name>
  # generate IsoAnnotLite output required for TAPPAS after sqanti (final files )
  isoannolite_generate cd $WKD_ROOT/9b_filter_cont/$2 $sample.classification.final.gtf" $sample.classification.txt $sample.junction.txt Human $WKD_ROOT/9b_filter_cont/$2 $1"_tappasannot.gff3"
  
  # counts_subset_4tappas <input_class> <output_class>
  counts_subset_4tappas $WKD_ROOT/9b_filter_cont/$2/$sample.classification.txt $WKD_ROOT/9b_filter_cont/$sample.expression.txt
}

# generate_rnaseq_counts <input_dir>
# input_dir = directory containing kallisto output files from alignment of RNA-Seq to Iso-Seq annotation file
generate_rnaseq_counts(){
  
  source activate nanopore
  # Rscript script.R <input.dir> <output.file> <type=Whole/Targeted/WholeTargeted> <targeted.class.files>
  Rscript $RNASEQCOUNT $1 $NAME"_RNASeq.expression.txt" Whole NA
  
}

# run_kallisto_1sample <input_RNASEQ_rawdir> <sample> <input_ref_name_idx> <output_dir>
run_kallisto_1sample(){
    source activate sqanti2_py3

    # variables
    input_dir=$1
    sample=$2
    input_ref_name=$3
    output_dir=$4
    
    # cd to $J20/Tg4510_input_directory
    cd $input_dir
    echo "Finding reads for Sample $sample for STAR"
    # find reverse and forward file, trimming "./" therefore only printing file name
    F=$(find . | grep "fastq" | grep $sample | grep "R1" | sed 's|^./||' )
    R=$(find . | grep "fastq" | grep $sample | grep "R2" | sed 's|^./||' )
    # save path directory of files as variable for later mapping
	  R1_READS=$(realpath $F)
    R2_READS=$(realpath $R)

    echo "Processing Kallisto for $sample"
    echo $R1_READS
    echo $R2_READS

    cd $output_dir
    #kallisto version
    time kallisto quant -i $input_ref_name --rf-stranded $R1_READS $R2_READS -o $sample 2> $2"_Kallisto.quant.log"
 
    source deactivate
}


# trim_and_run_kallisto <input_dir> <sample> <input_ref_name> <output_dir_trim> <output_dir_kallisto>
trim_and_run_kallisto(){
  source activate sqanti2
  
  # variables
  input_dir=$1
  sample=$2
  input_ref_name=$3
  output_dir_trim=$4
  output_dir_kallisto=$5
  
  echo "Trimming $sample"
  # find reverse and forward file, trimming "./" therefore only printing file name
  R1_READS=$(find $input_dir | grep "fastq" | grep $sample | grep "r1" | sed 's|^./||' )
  R2_READS=$(find $input_dir | grep "fastq" | grep $sample | grep "r2" | sed 's|^./||' )
  echo $R1_READS
  echo $R2_READS
  trim_galore --gzip --paired $R1_READS $R2_READS -o $output_dir_trim
  
  echo "Processing Kallisto for $sample"
  # "val" files are the final "validation" files generated from trim_galore
  R1_TRIMMED_READS=$(find $output_dir_trim | grep "val" | grep $sample | grep "r1" | sed 's|^./||' )
  R2_TRIMMED_READS=$(find $output_dir_trim | grep "val" | grep $sample | grep "r2" | sed 's|^./||' )
  echo $R1_TRIMMED_READS
  echo $R2_TRIMMED_READS
  cd $output_dir_kallisto
  time kallisto quant $R1_TRIMMED_READS $R2_TRIMMED_READS -i $input_ref_name --rf-stranded -o $sample 2> $sample"_Kallisto.quant.log"
  echo "All Done!"
  source deactivate
}

#!/bin/bash               

arg_s=1
arg_e=3
arg_ee=3
arg_C=1
RScrt="/home/16S_pipeline/code/Longitudinal.R"
arg_r=5
arg_M="LMM_NBMM_ZINBMM_mTMAT"
arg_n="Control_Case"
Tree_EZ="/usr/local/16S_pipeline/file/eztaxon_sina_fastmp.tre"
Tree_Silva="/usr/local/16S_pipeline/file/silva_132_97_otus.tre"
Tree_gg="/usr/local/16S_pipeline/file/97_otus_gg.tree"
arg_q=0.05
arg_l=2000

setDirectories(){
    INPUT_DIR=$(readlink -f $INPUT_DIR)
    OUTPUT_DIR=$(readlink -f $OUTPUT_DIR)
}

help() {
    echo "process [OPTIONS] FILE"
    echo "    -h         Help"
    echo "    -d ARG     Input dirname (S1)"
    echo "    -D ARG     Output dirname(S1-3)"
    echo "    -t ARG     Threads (S2)"
    echo "    -z ARG     Database choice: ez, silva, gg (S2-3, for S3 you can put it like this: ez_silva)" 
    echo "    -s ARG     Start step"
    echo "    -e ARG     End step"
    echo "    -M ARG     Methods choice. Default: LMM_NBMM_ZINBMM_mTMAT"
    echo "    -p ARG     Number of Permutation (S2)"
    echo "    -c ARG     Number of Columns of Plots (S2-3)"
    echo "    -P ARG     Phenotype file name. Just a name. Not a path. (S1)"
    echo "    -C ARG     Is Phenotype Continuous or Ordinal? If yes, 1. If not (Binary or categorical), 0. (S2)"
    echo "    -r ARG     Toxonomy rank. 1: phylum 2: Class 3: Order 4: Family 5: Genus 6: Species"
    echo "    -n ARG     Names of phenotype from 0 to k k=Number of levels of Phenotype
                         ex) -n Control_Case => 0:Control, 1:Case (S3)"
    echo "    -q ARG     FDR-q value for ANCOM"
    echo "    -l ARG     Library size cutoff"
    echo "    -V ARG     Phenotype variable"
    echo "    -v ARG     Covariates separated by comma. ex) -v race,sex "
    echo "    -Z ARG     Covariates for zero modeling for zero-inflated model. ex) -Z race,sex "
    echo "    -T ARG     Time variable by which you want to see the trend in the plot."
    echo ""
        exit 0
}

while getopts "d:D:t:z:s:e:M:p:c:P:C:r:n:q:l:V:v:Z:T:" opt
do
    case $opt in
        d) arg_d=$OPTARG
            INPUT_DIR=$arg_d
            ;;
	r) arg_r=$OPTARG
	    ;;
        D) arg_D=$OPTARG
	    OUTPUT_DIR=$arg_D
	    ;;
        t) arg_t=$OPTARG
            THREADS=$arg_t
            ;;
        z) arg_z=$OPTARG
            ;;
        s) arg_s=$OPTARG
           ;;
        e) arg_e=$OPTARG
           ;;
	M) arg_M=$OPTARG
           ;;
	p) arg_p=$OPTARG
           ;;
        c) arg_c=$OPTARG
           ;;
	P) arg_P=$OPTARG
	   ;;
	C) arg_C=$OPTARG
	   ;;
	n) arg_n=$OPTARG
	   ;;
	q) arg_q=$OPTARG
	   ;;
        l) arg_l=$OPTARG
	   ;;
        V) arg_V=$OPTARG
	   ;;
        v) arg_v=$OPTARG
	   ;;
        Z) arg_Z=$OPTARG
	   ;;
        T) arg_T=$OPTARG
	   ;;
        h) help ;;
        ?) help ;;
    esac
done


run(){
	
	if [ "$arg_z" == "gg" ];  then
		Tree=$Tree_gg
		input_DIR=$INPUT_DIR/otu_table_mc2_new_gg.txt
        elif [ "$arg_z" == "silva" ];  then
		Tree=$Tree_Silva
		input_DIR=$INPUT_DIR/otu_table_mc2_new_silva.txt
	elif [ "$arg_z" == "ez" ];  then
		Tree=$Tree_EZ
		input_DIR=$INPUT_DIR/otu_table_mc2_new_ez.txt
        fi

	pheno_DIR=$INPUT_DIR/$arg_P

        Rscript $RScrt $OUTPUT_DIR $input_DIR $arg_M $Tree $pheno_DIR $arg_z $arg_r $arg_C $arg_l $arg_V $arg_v $arg_Z $arg_T
}


mkdir -p $OUTPUT_DIR
setDirectories
run



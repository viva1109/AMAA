#!/bin/bash                                                                       

INPUT_DIR=""
TRIM_DIR=""
MERGE_DIR=""
FILTER_DIR=""
STAT_DIR=""
CHIM_DIR=""
TAX_DIR=""

MINLEN=350
MAXLEN=550 

THREADS=1
PERL=$(which perl)

PRIMER_FOR=""
PRIMER_REV=""

FILE_IDENTIFER_F="_1"
FILE_IDENTIFER_R="_2"

METHOD="silva"

S_STEP=1
E_STEP=12

CURRENT_NSEQ=0
CURRENT_NERR=0

CURRENT_STATUS=0

CONFIG_FILE=""
OUTFILE_ID="otu_table"

ALL_OUTPUT_DIR=""
ALL_LOG_DIR=""
ALL_LOG_FILE="all.log"
ALL_ERRLOG_FILE="all.error.log"

RAREFY_STRING="10,100"
OTU_UID="F"

OTU_PICK_METHOD="denovo"
TMP_OTU_BIOM_DIR=$OTU_PICK_METHOD/$METHOD
SKIP_STEP_TEN="F"

readConfig(){
    CONFIG_FILE=$(readlink -f $CONFIG_FILE)
    source $CONFIG_FILE

    MAKE_OTU=$(readlink -f $MAKE_OTU)
    MATCH_OTU=$(readlink -f $MATCH_OTU)
    PERL_SCRIPT=$(readlink -f $PERL_SCRIPT)
    BOKULICH=$(readlink -f $BOKULICH)
    ID_TO_FULL=$(readlink -f $ID_TO_FULL) 
    RE_REPL=$(readlink -f $RE_REPL)
    DADA2_R=$(readlink -f $DADA2_R)

    DB_SILVA_128=$(readlink -f $DB_SILVA_128)
    TAX_ID_SILVA_128=$(readlink -f $TAX_ID_SILVA_128)
    TAX_REF_SILVA_128=$(readlink -f $TAX_REF_SILVA_128)

    DB_SILVA_132=$(readlink -f $DB_SILVA_132)
    TAX_ID_SILVA_132=$(readlink -f $TAX_ID_SILVA_132)
    TAX_REF_SILVA_132=$(readlink -f $TAX_REF_SILVA_132)

    DB_EZ=$(readlink -f $DB_EZ)
    TAX_ID_EZ=$(readlink -f $TAX_ID_EZ)
    TAX_REF_EZ=$(readlink -f $TAX_REF_EZ)

    DB_GG=$(readlink -f $DB_GG)
    TAX_ID_GG=$(readlink -f $TAX_ID_GG)
    TAX_REF_GG=$(readlink -f $TAX_REF_GG)

    CHIMERA_REF=$(readlink -f $CHIMERA_REF)

}


setDirectories(){
    INPUT_DIR=$(readlink -f $INPUT_DIR)
    TRIM_DIR=$(dirname $INPUT_DIR)/"1.trimmed" 
    MERGE_DIR=$(dirname $INPUT_DIR)/"2.merged"    
    FILTER_DIR=$(dirname $INPUT_DIR)/"3.filter"   
    #STAT_DIR=$(dirname $INPUT_DIR)/"stat"   
    DEREP_DIR=$(dirname $INPUT_DIR)/"4.derep"  
    CHIM_DIR=$(dirname $INPUT_DIR)/"5.chimera"    
    OTU_DIR=$(dirname $INPUT_DIR)/"6.otu"
    if [ -z "$ALL_OUTPUT_DIR" ];then
        ALL_OUTPUT_DIR=$(dirname $INPUT_DIR)/"output"
        mkdir -p $ALL_OUTPUT_DIR
    else 
     #   ALL_OUTPUT_DIR=$(readlink -f $ALL_OUTPUT_DIR)
     #   echo "path $ALL_OUTPUT_DIR"
        mkdir -p $ALL_OUTPUT_DIR
        if [ ! -d $ALL_OUTPUT_DIR ];then
             echo "No directory : "$ALL_OUTPUT_DIR
             exit 1
	fi
        ALL_OUTPUT_DIR=$(readlink -f $ALL_OUTPUT_DIR)
    fi
    ALL_LOG_DIR=$ALL_OUTPUT_DIR/"log"
    mkdir -p $ALL_LOG_DIR
    ALL_LOG_FILE=$ALL_LOG_DIR/$ALL_LOG_FILE
    ALL_ERRLOG_FILE=$ALL_LOG_DIR/$ALL_ERRLOG_FILE
    echo -n "" > $ALL_LOG_FILE
    #touch $ALL_LOG_FILE
    #touch $ALL_ERRLOG_FILE
    echo -n "" > $ALL_ERRLOG_FILE
}

checkError()
{
if [ $CURRENT_STATUS -gt 0 ]; then
    
    echo '
   ***************
   *** ABORTED ***
   ***************
'
    echo "   An error occurred in Step $1"
    exit 1
fi
}



unzip(){ 
    cd $INPUT_DIR
    count=$(find ./ -maxdepth 1 -name "*.gz" | wc -l)
    
    printf '  number of gz files: %d\n' "$count"
    currentCount=1
    for file in $(find ./ -maxdepth 1 -type f);do
        if [[ "$file" == *.gz ]]; then
            gzip -d $file;
            abpath=$(readlink -e $file)
            printf '  %d/%d :%s\n\r' "$currentCount" "$count" "$abpath" 
            echo -e "$currentCount/$count $abpath" >> $ALL_LOG_FILE
            currentCount=$((currentCount+1))
        fi
    done 
}

obtainChimeraRef(){
    echo Obtaining Gold reference database for chimera detection

    if [ ! -e $CHIMERA_REF ]; then
        cd $(dirname $CHIMERA_REF)
        if [ ! -e Silva.gold.bacteria.zip ]; then
            wget https://www.mothur.org/w/images/f/f1/Silva.gold.bacteria.zip
    fi

    echo Decompressing and reformatting...

    echo "unzip -p Silva.gold.bacteria.zip silva.gold.align | sed -e "s/[.-]//g" > $(basename $CHIMERA_REF)"
    unzip -p Silva.gold.bacteria.zip silva.gold.align | sed -e "s/[.-]//g" > $(basename $CHIMERA_REF)
    fi
}

checkFQformat(){
    cd $INPUT_DIR
    count=$(find ./ -maxdepth 1  \( -name "*.fq" -or -name "*.fastq" \) | wc -l)
    printf '  number of FASTQ files: %d\n' "$count" 
    #cat test.fq | paste - - - - | awk -F "\t" '{ if (length($2) != length($4)) print $0 }' | wc -l
    currentCount=1
    for file in $(find ./ -maxdepth 1 -type f  \( -name "*.fq" -or -name "*.fastq" \) );do
	errorLine=$(cat $file | paste - - - - | awk -F "\t" '{ if (length($2) != length($4)) print $0 }' | wc -l)
        abpath=$(readlink -e $file)

        printf '\r                                                                                                  '
        printf '\r  %d/%d : %s - %d' "$currentCount" "$count" "$abpath" "$errorLine"
        echo -e "$currentCount/$count $abpath $errorLine"  3>&1 1>>${ALL_LOG_FILE} 2>&1
        #exec 3>&1 1>>${LOG_FILE} 2>&1

        if [ $errorLine -gt 0 ];then
            printf '\n    invalid FASTQ file : %s\n' "$abpath" 
            echo -e "   invalid FASTQ file $abpath" >> $ALL_ERRLOG_FILE 
            CURRENT_STATUS=$(($CURRENT_STATUS+1))
        fi
        currentCount=$((currentCount+1))
    done 

    printf '\r                                                                                                  '
}

checkPair(){
    cd $INPUT_DIR
    for file in $(find ./ -maxdepth 1 -type f  \( -name "*.fq" -or -name "*.fastq" \) );do
    #    echo $file
        if [[ "$file" == *$FILE_IDENTIFER_F* ]]; then
            local str=$(basename $file)
            local TMP_ID=${str//$FILE_IDENTIFER_F*/""}
            local FILE_EXT=${str//*$FILE_IDENTIFER_F/""}
            local REV_FILE=./$TMP_ID$FILE_IDENTIFER_R$FILE_EXT

            if [ ! -e $REV_FILE ]; then
                echo "PAIR ERROR : $REV_FILE not found."
                echo -e "PAIR ERROR : $REV_FILE not found." >> $ALL_ERRLOG_FILE
                CURRENT_STATUS=$(($CURRENT_STATUS+1))
            fi
        elif [[ "$file" == *$FILE_IDENTIFER_R* ]]; then
            local str=$(basename $file)
            local TMP_ID=${str//$FILE_IDENTIFER_R*/""}
            local FILE_EXT=${str//*$FILE_IDENTIFER_R/""}
            local FOR_FILE=./$TMP_ID$FILE_IDENTIFER_F$FILE_EXT

            if [ ! -e $FOR_FILE ]; then
                echo "PAIR ERROR : $FOR_FILE not found."         
                echo -e "PAIR ERROR : $FOR_FILE not found." >>$ALL_ERRLOG_FILE
                CURRENT_STATUS=$(($CURRENT_STATUS+1))
            fi
        else 
            echo "PAIR ERROR : $file has no pair" 
            echo -e "PAIR ERROR : $file has no pair" >> $ALL_ERRLOG_FILE 
        fi
    done
}

trimAdapter(){

    cd $INPUT_DIR
    local OUT_DIR=$TRIM_DIR
    mkdir -p $OUT_DIR
    mkdir -p $OUT_DIR/log    
    local LOG_DIR=$OUT_DIR/log 
    local totalcount=$(find ./ -maxdepth 1  \( -name "*.fq" -or -name "*.fastq" \) | wc -l)
    totalcount=$((totalcount/2)) 
    local currentCount=1

    for file in $(find ./ -maxdepth 1 -type f  \( -name "*.fq" -or -name "*.fastq" \) );do
        if [[ $file == *$FILE_IDENTIFER_F* ]]; then
     
            local str=$(basename $file)
            local TMP_ID=${str//$FILE_IDENTIFER_F*/""}
            local FILE_EXT=${str//*$FILE_IDENTIFER_F/""}
            local REV_FILE=./$TMP_ID$FILE_IDENTIFER_R$FILE_EXT
            
            local OUTFILE_FWD=$OUT_DIR/$TMP_ID$FILE_IDENTIFER_F$FILE_EXT
            local OUTFILE_REV=$OUT_DIR/$TMP_ID$FILE_IDENTIFER_R$FILE_EXT
            
	    if [ -f $OUTFILE_REV ]; then
          	  echo "file exist "$OUTFILE_REV
		  currentCount=$((currentCount+1))

		 continue
	    fi  
            printf '\r                                                                          '
            printf '\r  %d/%d :%s \n' "$currentCount" "$totalcount" "$TMP_ID"
            echo -e "$currentCount/$totalcount :$TMP_ID" >>$ALL_LOG_FILE
            #echo "$CUTADAPT -g $PRIMER_FOR -G $PRIMER_REV -o $OUTFILE_FWD -p $OUTFILE_REV $file $REV_FILE -O 11 -e 0.15 -m 10" 
            $CUTADAPT -g $PRIMER_FOR -G $PRIMER_REV -o $OUTFILE_FWD -p $OUTFILE_REV $file $REV_FILE -O 11 -e 0.15 -m 10 > $LOG_DIR/$TMP_ID".log" 2> $LOG_DIR/$TMP_ID".err"            
            #echo $CUTADAPT -g $PRIMER_FOR -G $PRIMER_REV -o $OUTFILE_FWD -p $OUTFILE_REV $file $REV_FILE -O 11 -e 0.15 -m 10
            cat $LOG_DIR/$TMP_ID".log" >> "$ALL_LOG_FILE"
            cat $LOG_DIR/$TMP_ID".err" >> "$ALL_ERRLOG_FILE"
 
            
            if [ -s $LOG_DIR/$TMP_ID".err" ]; then 
                printf ""
                echo "  FILE ERROR : "$TMP_ID  
                if grep -q improperly $LOG_DIR/$TMP_ID".err"; then
  		    echo  
		else 
                    CURRENT_STATUS=$(($CURRENT_STATUS+1))
                fi
                grep "cutadapt: error:" $LOG_DIR/$TMP_ID".err"
            else
                rm  $LOG_DIR/$TMP_ID".err"
            fi
            currentCount=$((currentCount+1))
        fi
    done
}

merge(){
    cd $TRIM_DIR
    local OUT_DIR=$MERGE_DIR
    mkdir -p $OUT_DIR
    mkdir -p $OUT_DIR/log
    local LOG_DIR=$OUT_DIR/log
    local totalcount=$(find ./ -maxdepth 1  \( -name "*.fq" -or -name "*.fastq" \) | wc -l)
    totalcount=$((totalcount/2))
    local currentCount=1

    for file in $(find ./ -maxdepth 1 -type f  )
    do
        if [[ $file == *$FILE_IDENTIFER_F* ]]; then
            local str=$(basename $file)
            local TMP_ID=${str//$FILE_IDENTIFER_F*/""}
            local FILE_EXT=${str//*$FILE_IDENTIFER_F/""}
            local REV_FILE=./$TMP_ID$FILE_IDENTIFER_R$FILE_EXT
            local OUT_FILE=$OUT_DIR/$TMP_ID
            printf '  %d/%d : %s \n' "$currentCount" "$totalcount" "$TMP_ID"
            echo -e "$currentCount/$totalcount : $TMP_ID" >>$ALL_LOG_FILE
            $CASPER $file $REV_FILE -t $THREADS -g 0.27 -o $OUT_FILE  > $LOG_DIR/$TMP_ID".log" 2> $LOG_DIR/$TMP_ID".err"
            #to default jellyfish error, when the error happens there is no jellyfish error 
            
            cat $LOG_DIR/$TMP_ID".log" >> "$ALL_LOG_FILE"
            cat $LOG_DIR/$TMP_ID".err" >> "$ALL_ERRLOG_FILE"

            if [ -s $LOG_DIR/$TMP_ID".err" ]; then
                rm  $LOG_DIR/$TMP_ID".err"
            else
                echo ">>ERROR : "$TMP_ID  
                cat  $LOG_DIR/$TMP_ID".log"
                CURRENT_STATUS=$(($CURRENT_STATUS+1))
            fi
            currentCount=$((currentCount+1)) 
        fi
    done
}

filter(){
    cd $MERGE_DIR
    local OUT_DIR=$FILTER_DIR
    mkdir -p $OUT_DIR
    
    local totalcount=$(find ./ -maxdepth 1  \( -name "*.fq" -or -name "*.fastq" \) | wc -l)
    local currentCount=1

    for file in $(find ./ -maxdepth 1 -type f  \( -name "*.fq" -or -name "*.fastq" \) );do
        printf '  %d/%d : %s \n' "$currentCount" "$totalcount" "$file"
        echo -e "$currentCount/$totalcount :$file" >>$ALL_LOG_FILE
        if [ -s $file ]; then
            local str=$(basename $file)
            local OUT_FILE=$OUT_DIR/$str
            python $BOKULICH -i $file -o $OUT_FILE -m $MINLEN -M $MAXLEN 
	else
            echo $file " size 0"
            echo -e "$file : file size 0" >>$ALL_ERRLOG_FILE
        fi
        currentCount=$((currentCount+1))
    done
}

dereplication(){
    cd $FILTER_DIR
    local OUT_DIR=$DEREP_DIR
    mkdir -p $OUT_DIR
    mkdir -p $OUT_DIR/log
    local LOG_DIR=$OUT_DIR/log
    local totalcount=$(find ./ -maxdepth 1  \( -name "*.fq" -or -name "*.fastq" \) | wc -l)
    local currentCount=1

    for file in $(find ./ -maxdepth 1 -type f  \( -name "*.fq" -or -name "*.fastq" \) );do
        printf '  %d/%d : %s \n' "$currentCount" "$totalcount" "$file"
        echo -e "$currentCount/$totalcount :$file" >>$ALL_LOG_FILE

        local TMP_ID=$(basename "$file" | cut -d. -f1)
        $VSEARCH --threads $THREADS \
        --derep_fulllength $file \
        --strand plus \
        --output $OUT_DIR/$TMP_ID.derep.fasta  \
        --sizeout \
        --uc $OUT_DIR/$TMP_ID.derep.uc \
        --relabel $TMP_ID. \
        --fasta_width 0 &> $LOG_DIR/$TMP_ID".log"
        cat $LOG_DIR/$TMP_ID".log" >> "$ALL_LOG_FILE"
        currentCount=$((currentCount+1))
    done
}


removeChiRef(){
    cd $DEREP_DIR
    local OUT_DIR=$CHIM_DIR
    mkdir -p $OUT_DIR
    mkdir -p $OUT_DIR/log
    local LOG_DIR=$OUT_DIR/log
    local totalcount=$(find ./ -maxdepth 1  -name "*.fasta" | wc -l)
    local currentCount=1

    for file in $(find ./ -maxdepth 1 -type f  -name "*.fasta" )
    do
        if [[ $file == *"fasta"* ]]; then
            printf '  %d/%d : %s \n' "$currentCount" "$totalcount" "$file"
            echo -e "$currentCount/$totalcount :$file" >>$ALL_LOG_FILE
            local TMP_ID=$(basename "$file" | cut -d. -f1)
            $VSEARCH --threads $THREADS \
                --uchime_ref $file \
                --db $CHIMERA_REF \
                --sizein \
                --sizeout \
                --fasta_width 0 \
                --nonchimeras $OUT_DIR/$TMP_ID.derep.chim.fasta  &> $LOG_DIR/$TMP_ID".log"
            cat $LOG_DIR/$TMP_ID".log" >> "$ALL_LOG_FILE"
        fi
        currentCount=$((currentCount+1))
    done

    rm -f all.nonchimeras.derep.fasta
    cat $OUT_DIR/*.derep.chim.fasta > $OUT_DIR/all.nonchimeras.derep.fasta
}


redoreplication(){    
    cd $CHIM_DIR
    local OUT_DIR=$CHIM_DIR
    cat *.derep.chim.fasta > $OUT_DIR/all.nonchimeras.derep.fasta

    local INPUT_FILE=$OUT_DIR/all.nonchimeras.derep.fasta
    local FASTA_FILE=$(readlink -f $INPUT_FILE)

    $VSEARCH --threads $THREADS \
    --derep_fulllength $FASTA_FILE \
    --strand plus \
    --output all.tbc.nonchi.derep.fasta  \
    --sizein \
    --sizeout \
    --uc all.tbc.nonchi.derep.uc \
    --fasta_width 0   

    echo Processed file all.tbc.nonchi.derep.fasta
}


cluster(){
    cd $CHIM_DIR
    local OUT_DIR=$OTU_DIR
    mkdir -p $OUT_DIR
    mkdir -p $OUT_DIR/log 
    local LOG_DIR=$OUT_DIR/log
    # mkdir -p $TMP_OTU_BIOM_DIR
    echo "   clustering : $OTU_PICK_METHOD "

    echo -e "clustering : $OTU_PICK_METHOD " >>$ALL_LOG_FILE  
    if [ "$OTU_PICK_METHOD" == "denovo" ];  then
    	$VSEARCH --threads $THREADS \
        	--cluster_size all.nonchimeras.derep.fasta \
	        --id 0.97 \
	        --strand plus \
	        --sizein \
	        --fasta_width 0 \
	        --uc $OUT_DIR/all.clustered.uc \
	        --relabel OTU_ \
	        --centroids $OUT_DIR/all.otu.centroids.fasta \
	        --otutabout $OUT_DIR/all.otutab.txt 
        echo         $VSEARCH --threads $THREADS \
                --cluster_size all.nonchimeras.derep.fasta \
                --id 0.97 \
                --strand plus \
                --sizein \
                --fasta_width 0 \
                --uc $OUT_DIR/all.clustered.uc \
                --relabel OTU_ \
                --centroids $OUT_DIR/all.otu.centroids.fasta \
                --otutabout $OUT_DIR/all.otutab.txt 
        continue
    elif [ "$OTU_PICK_METHOD" == "closed" ] ||  [ "$OTU_PICK_METHOD" == "open" ];  then
        local DB=""
        local TAX=""
        local TAX_REF=""
#        local BIOM=$ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom
#        local OUT_DIR=$TMP_OTU_BIOM_DIR
        python $RE_REPL -i all.nonchimeras.derep.fasta -o all.nonchi.fasta
        source activate qiime1        
        if [ "$METHOD" == "silva" ];  then
            DB=$DB_SILVA_132
            TAX_REF=$TAX_REF_SILVA_132
	    TAX=$TAX_ID_SILVA_132
        elif [ "$METHOD" == "silva128" ];  then
           DB=$DB_SILVA_128
           TAX_REF=$TAX_REF_SILVA_128
           TAX=$TAX_ID_SILVA_128
	    elif [ "$METHOD" == "ez" ];  then
           DB=$DB_EZ
           TAX_REF=$TAX_REF_EZ
           TAX=$TAX_ID_EZ
        elif [ "$METHOD" == "gg" ];  then
            DB=$DB_GG
            TAX_REF=$TAX_REF_GG
            TAX=$TAX_ID_GG
         fi

         if [ "$OTU_PICK_METHOD" == "open" ];  then
              echo "pick_otus:threads	$THREADS" > $OUT_DIR/open_setting.txt
              pick_open_reference_otus.py -m usearch61 -i all.nonchi.fasta -o $OUT_DIR -a --jobs_to_start  $THREADS  -r $DB --suppress_align_and_tree --suppress_taxonomy_assignment -f -p $OUT_DIR/open_setting.txt --min_otu_size 1 > $LOG_DIR/tmp.log 2> $LOG_DIR/tmp.err
	      mv $OUT_DIR/rep_set.fna $OUT_DIR/all.otu.centroids.fasta
	      biom convert -i $OUT_DIR/otu_table_mc1.biom  -o $OUT_DIR/all.otutab.txt  --to-tsv
              rm $OUT_DIR/otu_table_mc1.biom 
       
         elif [ "$OTU_PICK_METHOD" == "closed" ]; then
              parallel_pick_otus_uclust_ref.py -i all.nonchi.fasta -o $OUT_DIR/closed_otus -r $DB -Z 2  --jobs_to_start $THREADS  > $LOG_DIR/tmp.log 2> $LOG_DIR/tmp.err
              make_otu_table.py -i $OUT_DIR/closed_otus/all.nonchi_otus.txt -t $TAX_REF -o  $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom > $LOG_DIR/tmp2.log 2> $LOG_DIR/tmp2.err
              biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.txt --to-tsv --header-key taxonomy
              cat $LOG_DIR/tmp2.log >> $LOG_DIR/tmp.log
              cat $LOG_DIR/tmp2.err >> $LOG_DIR/tmp.err
              
              if [ "$OTU_UID" == "T" ]; then
                   make_otu_table.py -i $OUT_DIR/closed_otus/all.nonchi_otus.txt -t $TAX -o  $ALL_OUTPUT_DIR/"$OUTFILE_ID".biom > $LOG_DIR/tmp3.log 2> $LOG_DIR/tmp3.err
                   cat $LOG_DIR/tmp3.log >> $LOG_DIR/tmp.log
                   cat $LOG_DIR/tmp3.err >> $LOG_DIR/tmp.err
              fi              
         fi
         cat $LOG_DIR/tmp.log >> "$ALL_LOG_FILE"
         cat $LOG_DIR/tmp.err >> "$ALL_ERRLOG_FILE"
    elif [ "$OTU_PICK_METHOD" == "asv" ];  then
        echo "redo replication for asv"
        redoreplication
    fi
}

assign_tax(){
    local DB=""
    local TAX=""
    source activate qiime1
    local INPUT_FILE=$OTU_DIR/all.otu.centroids.fasta

    echo -e "database : $METHOD" >>$ALL_LOG_FILE
    if [ "$OTU_PICK_METHOD" == "closed" ]; then
       echo "	no need for closed picking"
       continue
    fi

    if [ "$METHOD" == "silva" ];  then
        echo database is silva
        DB=$DB_SILVA_132
        TAX=$TAX_ID_SILVA_132
        TAX_DIR=$OTU_DIR/"$METHOD"
    elif [ "$METHOD" == "silva128" ];  then
        echo database is silva
        DB=$DB_SILVA_128
        TAX=$TAX_ID_SILVA_128
        TAX_DIR=$OTU_DIR/"$METHOD"
    elif [ "$METHOD" == "ez" ];  then
        echo database is ez
        DB=$DB_EZ
        TAX=$TAX_ID_EZ
        TAX_DIR=$OTU_DIR/"$METHOD"
    elif [ "$METHOD" == "gg" ];  then
        echo database is green gene
        DB=$DB_GG
        TAX=$TAX_ID_GG
        TAX_DIR=$OTU_DIR/"$METHOD"
    else 
	echo "Invalied database : $METHOD"
        exit 
    fi
    
    if [ "$OTU_PICK_METHOD" == "asv" ]; then
        INPUT_FILE=$CHIM_DIR/"all.tbc.nonchi.derep.fasta"
    fi

    parallel_assign_taxonomy_uclust.py \
                                -i $INPUT_FILE \
                                -o $TAX_DIR \
                                --reference_seqs_fp $DB \
                                --id_to_taxonomy_fp $TAX \
                                -O $THREADS \
                                --uclust_max_accepts 1 \
                                -Z 20 -T \
                                --min_consensus_fraction 1
}


post_otu(){
    cd $OTU_DIR
    local TAX_DIR=$OTU_DIR/"$METHOD"
    local FILE_ID="all.otu.centroids"
    source activate qiime1
    echo "make OTU table"
    echo -e "make OTU table" >>$ALL_LOG_FILE
    OTU_UID=$(echo $OTU_UID | tr 'a-z' 'A-Z')   
    
    if [ "$OTU_PICK_METHOD" == "closed" ]; then
        filter_otus_from_otu_table.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -n 2
        biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.txt --to-tsv --header-key taxonomy
        echo "summarize taxa"
        echo -e "summarize taxa" >>$ALL_LOG_FILE
        summarize_taxa.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -o  $ALL_OUTPUT_DIR"/summary" -L 2,3,4,5,6,7  
	 
        if [ "$OTU_UID" == "T" ]; then
             filter_otus_from_otu_table.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID".biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2.biom -n 2
             biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2.txt --to-tsv --header-key taxonomy
	    fi 

    elif [ "$OTU_PICK_METHOD" == "asv" ];  then
            local TAX_REF=$TAX_REF_SILVA_132
	    if [ "$METHOD" == "gg" ];  then
                TAX_REF=$TAX_REF_GG
            elif [ "$METHOD" == "silva128" ];  then
                TAX_REF=$TAX_REF_SILVA_128
	    elif [ "$METHOD" == "ez" ];  then
                TAX_REF=$TAX_REF_EZ
            fi
            FILE_ID="all.tbc.nonchi.derep"
            python $ID_TO_FULL -t $TAX_DIR/"$FILE_ID"_tax_assignments.txt -r $TAX_REF -o  $TAX_DIR/"$FILE_ID"_tax_assignments_w_tax.txt
            TAX_FILE=$TAX_DIR/"$FILE_ID"_tax_assignments_w_tax.txt
            UC_FILE=$CHIM_DIR/"all.tbc.nonchi.derep.uc"
            FASTA_FILE=$CHIM_DIR/"all.tbc.nonchi.derep.fasta"
            echo  python $MAKE_OTU -u $UC_FILE -f $FASTA_FILE -t $TAX_FILE -c $THREADS -o $ALL_OUTPUT_DIR
            python $MAKE_OTU -u $UC_FILE -f $FASTA_FILE -t $TAX_FILE -c $THREADS -o $ALL_OUTPUT_DIR
            biom convert -i $ALL_OUTPUT_DIR/"otu_table_w.txt" -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy
            biom convert -i $ALL_OUTPUT_DIR/"otu_table_mc2_w.txt" -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy
            return 
    else
        local TAX_REF=$TAX_REF_SILVA_132
        if [ "$METHOD" == "gg" ];  then 
           TAX_REF=$TAX_REF_GG
	    elif [ "$METHOD" == "silva128" ];  then
           TAX_REF=$TAX_REF_SILVA_128
        elif [ "$METHOD" == "ez" ];  then
           TAX_REF=$TAX_REF_EZ
        fi
    
        python $ID_TO_FULL -t $TAX_DIR/"$FILE_ID"_tax_assignments.txt -r $TAX_REF -o  $TAX_DIR/"$FILE_ID"_tax_assignments_w_tax.txt
    
        if [ "$OTU_UID" == "T" ]; then
             echo "otu table will be generated with taxonomy uid"
            python $MATCH_OTU -t $TAX_DIR/"$FILE_ID"_tax_assignments.txt -m all.otutab.txt  -o $ALL_OUTPUT_DIR/"$OUTFILE_ID".txt
            biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID".txt -o $ALL_OUTPUT_DIR/"$OUTFILE_ID".biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy
            filter_otus_from_otu_table.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID".biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2.biom -n 2
            biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2.txt --to-tsv --header-key taxonomy

            python $MATCH_OTU -t $TAX_DIR/"$FILE_ID"_tax_assignments_w_tax.txt -m all.otutab.txt  -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.txt
            biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.txt -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy
            filter_otus_from_otu_table.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -n 2           
            biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.txt --to-tsv --header-key taxonomy

            echo "summarize taxa"
            echo -e "summarize taxa" >>$ALL_LOG_FILE
            summarize_taxa.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -o  $ALL_OUTPUT_DIR"/summary" -L 2,3,4,5,6,7   
            rm $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom 
        else    
            python $MATCH_OTU -t $TAX_DIR/"$FILE_ID"_tax_assignments_w_tax.txt -m all.otutab.txt  -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.txt
            biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.txt -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy
            filter_otus_from_otu_table.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_w.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -n 2
            biom convert -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -o $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.txt --to-tsv --header-key taxonomy  

            echo "summarize taxa"
            echo -e "summarize taxa" >>$ALL_LOG_FILE
            summarize_taxa.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -o  $ALL_OUTPUT_DIR"/summary" -L 2,3,4,5,6,7
        fi
    fi
    
    echo "make DATA summary"
    echo -e "make DATA summary" >>$ALL_LOG_FILE
    biom table-ids -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom  > $ALL_OUTPUT_DIR/list    
    biom summarize-table -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -o $ALL_OUTPUT_DIR/"table_summary.txt"

    # local MAX_LINE=$(grep -i  "Max:" $ALL_OUTPUT_DIR/"table_summary.txt") 
    # MAX_LINE=${MAX_LINE%.*}
    # local MAX_COUNT=$(echo $MAX_LINE | sed 's/[^0-9]*//g')
    # local INTERVAL=$(( (MAX_COUNT - 10)/10 ))
    
    # declare -a avalue=(10)
    # for ((i=1;i<=10;i++)); do
    #     avalue[$i]=$((10+(i * INTERVAL)))
    # done

    # avalue=(${avalue[@]} 500 1500 3000 5000)
    # sorted=($( printf "%s\n" "${avalue[@]}" | sort -un ))
    # RAREFY_STRING=$(printf ",%s" "${sorted[@]}")
    # RAREFY_STRING=${RAREFY_STRING:1}
    
    # echo "calculate chao1"
    # echo -e "calculate chao1" >>$ALL_LOG_FILE
    # Alpha rarefaction command 
    # multiple_rarefactions.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID".biom -m 10 -x $MAX_COUNT -s $INTERVAL -o $ALL_OUTPUT_DIR"/arare_max100/rarefaction/"  
#     multiple_rarefactions2.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -r $RAREFY_STRING -o $ALL_OUTPUT_DIR"/arare_max100/rarefaction/"   
#     #echo multiple_rarefactions2.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID".biom -r $RAREFY_STRING -o $ALL_OUTPUT_DIR"/arare_max100/rarefaction/"

# # Alpha diversity on rarefied OTU tables command 
#     alpha_diversity.py -i $ALL_OUTPUT_DIR"/arare_max100/rarefaction/" -o $ALL_OUTPUT_DIR"/arare_max100/alpha_div/"  -m chao1

#     # Collate alpha command 
#     collate_alpha.py -i $ALL_OUTPUT_DIR"/arare_max100/alpha_div/" -o $ALL_OUTPUT_DIR"/arare_max100/alpha_div_collated/" 

#     # Removing intermediate files command 
#     rm -r $ALL_OUTPUT_DIR"/arare_max100/rarefaction/" $ALL_OUTPUT_DIR"/arare_max100/alpha_div/"

#    echo "summarize taxa"
#    echo -e "summarize taxa" >>$ALL_LOG_FILE
#    summarize_taxa.py -i $ALL_OUTPUT_DIR/"$OUTFILE_ID".biom -o  $ALL_OUTPUT_DIR"/summary" -L 2,3,4,5,6,7

#    echo "parallel align seqs_pynast"
#    parallel_align_seqs_pynast.py -i all.otu.centroids.fasta -o pynast_aligned_seqs -T --jobs_to_start $THREADS

#    echo "make phylogeny"
#    qiime2
#    export OMP_NUM_THREADS=$THREADS
#    $FASTREE -fastest -nt pynast_aligned_seqs/"$FILE_ID"_aligned.fasta > rep_set.tre    
}

draw_tree(){
    cd $OTU_DIR
    local FILE_ID="all.otu.centroids"
    
    source activate qiime1
    filter_fasta.py -f "$FILE_ID".fasta -b $ALL_OUTPUT_DIR/"$OUTFILE_ID"_mc2_w.biom -o "$FILE_ID"_mc2.fasta 

    echo "parallel align seqs_pynast"
    parallel_align_seqs_pynast.py -i "$FILE_ID"_mc2.fasta -o pynast_aligned_seqs -T --jobs_to_start $THREADS

    echo "make phylogeny"
#    qiime2
    export OMP_NUM_THREADS=$THREADS
    $FASTREE -fastest -nt pynast_aligned_seqs/"$FILE_ID"_mc2_aligned.fasta > $ALL_OUTPUT_DIR"/rep_set_mc2.tre"
}
           
run1(){
    echo "Step 1 : decompress by gzip"
    #echo -e "\nStep 1 : decompress by gzip" >> $ALL_LOG_FILE
    unzip
}
run2(){
    echo "Step 2 : check fastq format"
    #echo -e "\nStep 2 : check fastq format" >> $ALL_LOG_FILE
    checkFQformat
}    
run3(){
    echo "Step 3 : check pair"
    #echo -e "\nStep 3: check pair" >> $ALL_LOG_FILE
    checkPair
}
run4(){
    echo "Step 4 : trim adapter"
    #echo -e "\nStep 4: trim adapter" >> $ALL_LOG_FILE
    trimAdapter
}
run5(){
    echo "Step 5 : merge paired reads"
    #echo -e "\nStep 4: trim adapter" >> $ALL_LOG_FILE
    merge
}
run6(){
    echo "Step 6 : filter by length"
    filter
}
run7(){
    echo "Step 7 : dereplication"
    dereplication
}
run8(){
    echo "Step 8 : remove chimera"
    removeChiRef
}
run9(){
    echo "Step 9 : clustering"
    cluster
}
run10(){ 
    echo "Step 10 : assign taxonomy"
    assign_tax
}
run11(){
    echo "Step 11 : post analysis"
    post_otu
}
run12(){
    if [ "$OTU_PICK_METHOD" == "open" ] ||  [ "$OTU_PICK_METHOD" == "denovo"  ]; then
	echo "Step 12 : draw OTU phylogeny "
        draw_tree
    elif [ "$OTU_PICK_METHOD" == "closed" ];  then
	echo "Step 12 : No need to draw OTU phylogeny. Use the reference Tree."
    fi
}

dadarun5(){
    echo "Step 5 : execute dada2 "
    Rscript $DADA2_R $TRIM_DIR $FILE_IDENTIFER_F $FILE_IDENTIFER_R $METHOD $THREADS 1 $ALL_OUTPUT_DIR $DADA_fMaxEE $DADA_rMaxEE $DADA_fTL $DADA_rTL
}
dadarun6(){
    echo "Step 6 : taxonomy assignment dada2 "
    Rscript $DADA2_R $TRIM_DIR $FILE_IDENTIFER_F $FILE_IDENTIFER_R $METHOD $THREADS 2 $ALL_OUTPUT_DIR $DADA_fMaxEE $DADA_rMaxEE $DADA_fTL $DADA_rTL
}

run(){
    mkdir -p /tmp/jobs
    chmod 777 -R /tmp/jobs
    
    if [ "$OTU_PICK_METHOD" == "dada" ]; then
        E_STEP=6
        dadasteps=0
        while [ $dadasteps -le $E_STEP ]
            do
            dadasteps=$(($dadasteps+1))
            if [ $dadasteps -ge $S_STEP ] && [ $dadasteps -le $E_STEP ]; then
                if [ $dadasteps -le 4 ]; then 
                    run$dadasteps
                else
                    dadarun$dadasteps
	        fi 
            echo ""
            fi
        done 
    else
        steps=0
        while [ $steps -le $E_STEP ]
            do
            steps=$(($steps+1))
            if [ $steps -ge $S_STEP ] && [ $steps -le $E_STEP ]; then
                echo -e "\nStep $steps" >> $ALL_LOG_FILE
                echo -e "\nStep $steps" >> $ALL_ERRLOG_FILE
                run$steps
                checkError $steps
            echo ""
            fi
        done
    fi 
}

printOptions(){
    echo "Input parameters:"
    echo "   inputdir     i) $INPUT_DIR"
    echo "   outputdir    o) $ALL_OUTPUT_DIR"
    echo "   f primer     f) $PRIMER_FOR"
    echo "   r primer     r) $PRIMER_REV"
    echo "   f identifer  F) $FILE_IDENTIFER_F"
    echo "   r identifer  R) $FILE_IDENTIFER_R"
    echo "   threads      t) $THREADS"
    echo "   database     z) $METHOD"
    echo "   min length   m) $MINLEN"
    echo "   max length   M) $MAXLEN"
    echo "   output id    u) $OUTFILE_ID"
    echo "   start step   s) $S_STEP"
    echo "   end step     e) $E_STEP"
    echo "   otu picking  p) $OTU_PICK_METHOD"
#echo "   rarefy range a) $RAREFY_STRING"
    echo ""
}

help(){
    echo " Usage: ./mgp [options]  (-d & -c are mandatory)"
    echo "    -h         help"
    echo "    -d STR     input directory where FASTQ files are located "
    echo "    -o STR     output directory [output]  "
    echo "    -f STR     forward primer seq"
    echo "    -r STR     reverse primer seq"
    echo "    -F STR     identifier for FASTQ file with forward reads"
    echo "    -R STR     identifier for FASTQ file with reverse reads"
    echo "    -t INT     threads [1]"
    echo "    -z         database : silva, silva128, gg, ez [silva]"
    echo "    -m INT     minimum read length after merge [350]"
    echo "    -M INT     maximum read length after merge [550]"
    echo "    -u STR     output table name [otu_table]"
    echo "    -s INT     start step [1]"
    echo "    -p STR     pick otu method :denovo,closed,open,dada [closed]"
    echo "    -e INT     end step [12]"
    echo "    	  step  1: unzip  "
    echo "          step  2: check file format"
    echo "          step  3: check file pair"
    echo "    	  step  4: trim adapter"
    echo "    	  step  5: merge"
    echo "    	  step  6: filtering by length and quality"
    echo "    	  step  7: dereplacation"
    echo "    	  step  8: remove chimera"
    echo "    	  step  9: clustering"
    echo "    	  step 10: assign taxonomy"
    echo "    	  step 11: post analysis"
    echo "          step 12: draw phylogeny"
    echo "    -c FILE    config file name *"
    echo "    -w STR     db uid output [F]"
    #echo "    -a list    rarefy range ex.(10,20,100)"
    exit 0
}

while getopts "d:f:r:F:R:t:z:m:M:s:e:c:u:o:a:w:p:" opt
do
    case $opt in
        d) arg_d=$OPTARG
            INPUT_DIR=$arg_d
            ;;
        f) arg_f=$OPTARG
            PRIMER_FOR=$arg_f
            ;;
        r) arg_r=$OPTARG
            PRIMER_REV=$arg_r
            ;;
        F) arg_F=$OPTARG
            FILE_IDENTIFER_F=$arg_F
            ;;
        R) arg_R=$OPTARG
            FILE_IDENTIFER_R=$arg_R
            ;;    
        t) arg_t=$OPTARG
            THREADS=$arg_t
            ;;
        z) arg_z=$OPTARG
            METHOD=$arg_z
            ;;
        m) arg_m=$OPTARG
            MINLEN=$arg_m
            ;;
        M) arg_M=$OPTARG
            MAXLEN=$arg_M
            ;;
        s) arg_s=$OPTARG
            S_STEP=$arg_s
           ;;
        e) arg_e=$OPTARG
            E_STEP=$arg_e 
           ;;
        c) arg_c=$OPTARG
            CONFIG_FILE=$arg_c
           ;;        
        u) arg_u=$OPTARG
            OUTFILE_ID=$arg_u
            ;;
        o) arg_o=$OPTARG
            ALL_OUTPUT_DIR=$arg_o
            ;;
        a) arg_a=$OPTARG
            RAREFY_STRING=$arg_a
            ;;
        w) arg_w=$OPTARG
            OTU_UID=$arg_w
            ;; 
        p) arg_p=$OPTARG
            OTU_PICK_METHOD=$arg_p
            ;;       
        h) help ;;
        ?) help ;;
    esac
done

echo $arg_z
if [ -z "$arg_d" ] || [ -z "$arg_c" ]; then
    "d & c should be specified"
    help
fi

#source activate qiime1
readConfig
setDirectories
printOptions
run

# echo "$file"



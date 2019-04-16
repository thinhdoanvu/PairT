#!/usr/bin/env bash
export LC_ALL=en_US.UTF-8
export SHELL=bash

VERSION='Version 0.1'
#This script serves as an interactive bash wrapper to QC, map, and estimate coverage of sequence from PE data.
#It requires that your raw data are split up by tagged individual and follow the naming convention of:

#Pop_Sample1.F.fastq and Pop_Sample1.R.fastq
#Pop_Sample2.F.fastq and Pop_Sample2.R.fastq
#...
#Pop_Samplek.F.fastq and Pop_Samplek.R.fastq

#Prints out title and contact info
echo -e "PairT" $VERSION
echo -e "Contact thinhdv@ntu.edu.vn with any problems \n\n "

###Code to check for the required software
echo "Checking for required software: bwa samtools bamToBed bedtools parallel fastp"
SOF=(bwa samtools bamToBed bedtools parallel fastp)
NUMSOF=0
for i in "${SOF[@]}"
do
	if which $i &> /dev/null; then
		foo=0
	else
    		echo "The" $i "software is not installed or is not in your" '$PATH'"."
    		NUMSOF=$((NUMDEP + 1))
	fi
done

#Checking samtools
SAMV1=$(samtools 2>&1 >/dev/null | grep Ver | sed -e 's/Version://' | cut -f2 -d " " | sed -e 's/-.*//' | cut -c1)
SAMV2=$(samtools 2>&1 >/dev/null | grep Ver | sed -e 's/Version://' | cut -f2 -d " " | sed -e 's/-.*//' | cut -c3)
	if [ "$SAMV1"  -ge "1" ]; then
		if [ "$SAMV2"  -lt "2" ]; then
        	echo "The version of Samtools installed in your" '$PATH' "is not optimized."
        	echo "Please install at least version 1.2.0"
			echo -en "\007"
			echo -en "\007"
			exit 1
		fi
	
	else
		    echo "The version of Samtools installed in your" '$PATH' "is not optimized."
        	echo "Please install at least version 1.3.0"
			echo -en "\007"
			echo -en "\007"
			exit 1
	fi
#Checking bwa
BWAV=$(bwa 2>&1 | mawk '/Versi/' | sed 's/Version: //g' | sed 's/0.7.//g' | sed 's/-.*//g' | cut -c 1-2)
	if [ "$BWAV" -lt "13" ]; then
        	echo "The version of bwa installed in your" '$PATH' "is not optimized."
        	echo "Please install at least version 0.7.13"
        	exit 1
	fi

BTC=$( bedtools --version | mawk '{print $2}' | sed 's/v//g' | cut -f1,2 -d"." | sed 's/2\.//g' )
	if [ "$BTC" -ge "26" ]; then
		BEDTOOLSFLAG="NEW"
		elif [ "$BTC" == "23" ]; then
		BEDTOOLSFLAG="OLD"
		elif [ "$BTC" != "23" ]; then
		echo "The version of bedtools installed in your" '$PATH' "is not optimized."
		echo "Please install version 2.23.0 or version 2.26.0 and above"
		exit 1	
	fi
#Checking fastp
FASTP=$(fastp -v 2>&1 | cut -f2 -d " ")
FASTP1=$(echo $FASTP | cut -f1 -d ".")
FASTP2=$(echo $FASTP | cut -f2 -d ".")
FASTP3=$(echo $FASTP | cut -f3 -d ".")
	if [ "$FASTP1" -lt "2" ]; then
		if [ "$FASTP2" -lt "20" ]; then
			if [ "$FASTP2" -lt "5" ]; then
				echo "The version of fastp installed in your" '$PATH' "is not optimized for dDocent."
				echo "Please install version 0.19.5 or above"
				exit 1
			fi
		fi
	fi
	
if ! sort --version | fgrep GNU &>/dev/null; then
	sort=gsort
else
	sort=sort
fi

if [ $NUMSOF -gt 0 ]; then
	echo -e "\nPlease install all required software before running dDocent again."
	exit 1
else
	echo -e "\nAll required software is installed!"
fi

#This code checks for individual fastq files follow the correct naming convention and are gziped
TEST=$(ls *.fastq 2> /dev/null | wc -l )
if [ "$TEST" -gt 0 ]; then
	echo -e "\nPairT is now configured to work on compressed sequence files.  Please run gzip to compress your files."
	echo "Do you want me to compress files for you? yes/no"

	read ans_TEST

	if [ "$ans_TEST" == "yes" ] 
	then
		gzip *fastq
	else
		echo "This is as simple as 'gzip *.fastq'"
		echo "Please rerun PairT after compressing files."
		exit 1
	fi
fi
#Count number of individuals in current directory
NumInd=$(ls *.F.fastq.gz 2> /dev/null | wc -l)
NumInd=$(($NumInd - 0))

#Create list of sample names
ls *.F.fastq.gz > namelist 2> /dev/null
sed -i'' -e 's/.F.fastq.gz//g' namelist

#Create an array of sample names
NUMNAMES=$(mawk '/_/' namelist | wc -l)

if [ "$NUMNAMES" -eq "$NumInd" ]; then
	NAMES=( `cat "namelist" `)
else
	echo "Individuals do not follow the PairT naming convention."
	echo "Please rename individuals to: Pop_Individual.F.fastq.gz"
	echo "For example: PopA_001.F.fastq.gz"
	exit 1
fi

if [[ "$1" == "help" || "$1" == "-help" || "$1" == "--help" || "$1" == "-h" || "$1" == "--h" ]]; then

	echo -e "\nTo run PairT, simply type '"PairT"' and press [ENTER]"
	echo -e "\nAlternatively, PaiRT can be run with a configuration file.  Usuage is:"
	echo -e "\nPairT PairT.config\n\n"
	exit 0
fi

#checking reference.fasta file
reffile=$(ls reference.fasta 2> /dev/null | wc -l)
reffile=$(($reffile - 0))
if [ $reffile -gt 0 ]; then
	echo -e "Sucessfuly found " $reffile " refence.fasta file for PairT"
else
	echo "Could not found any reference file for PairT"
fi

#Wrapper for main program functions.  This allows the entire file to be read first before execution
main(){
##########User Input Section##########
#This code gets input from the user and assigns variables
######################################

#Sets a start time variable
STARTTIME=$(date)

echo -e "\nPairT run started" $STARTTIME "\n"

#Checks if a configuration file is being used, if not asks for user input
if [ -n "$1" ]; then
	CONFIG=$1
	if [ ! -f $CONFIG ]; then
		echo -e "\nThe configuration file $CONFIG does not exist."
		exit 1
	fi

	#Checking CPU
	NUMProc=$(grep -A1 Processor $CONFIG 2> /dev/null | tail -1 ) 
	if [[ $NUMProc -lt 20 && $NUMProc -gt 1 ]]; then 
		MAXMemory1=$(grep -A1 Memory $CONFIG | sed 's/[g,G]//g' | tail -1)
	else
		echo -e "\nConfiguration file is not properly configured. This is a link: "
		exit 1
	fi

	#Checking MEMORY
	MAXMemory=$(( $MAXMemory1 / $NUMProc ))G
	if [[ "$OSTYPE" == "linux"* ]]; then
		MAXMemory=0
		MAXMemory1=0
	fi


	#Checking trimming or not
	TRIM=$(grep -A1 Trim $CONFIG | tail -1)

	#Checking for mapping or not	
	MAP=$(grep -A1 Mapping_R $CONFIG | tail -1)
	optA=$(grep -A1 _Match $CONFIG | tail -1)
	optB=$(grep -A1 MisMatch $CONFIG | tail -1)
	optO=$(grep -A1 Gap $CONFIG | tail -1)
	SNP=$(grep -A1 SNP $CONFIG | tail -1)

else
	GetInfo
fi

#Creates (or appends to) a run file recording variables
echo "Variables used in PairT Run at" $STARTTIME >> PairT.runs
echo "Number of Processors" >> PairT.runs
echo $NUMProc >> PairT.runs
echo "Maximum Memory" >> PairT.runs
echo $MAXMemory1 | sed 's/[g,G]//g' >> PairT.runs
echo "Trimming" >> PairT.runs
echo $TRIM >> PairT.runs
echo "Mapping_Reads?" >> PairT.runs
echo $MAP >> PairT.runs
echo "Mapping_Match_Value" >> PairT.runs
echo $optA >> PairT.runs
echo "Mapping_MisMatch_Value" >> PairT.runs
echo $optB >> PairT.runs
echo "Mapping_GapOpen_Penalty" >> PairT.runs
echo $optO >> PairT.runs
echo "Calling_SNPs?" >> PairT.runs
echo $SNP >> PairT.runs

##DoIT
if [[ "$TRIM" == "yes" ]]; then
        echo -e "\nTrimming reads"
        TrimReads 2> trim.log
fi


#Checks to see if reads will be mapped.
if [ "$MAP" != "no" ]; then
	echo -e "\nUsing BWA to map reads"
	if [ reference.fasta -nt reference.fasta.fai ]; then
	        samtools faidx reference.fasta &> index.log
	        bwa index reference.fasta >> index.log 2>&1
	fi
	
	#PairT now checks for trimmed read files before attempting mapping
        if [[ "$MAP" != "no" && ! -f "${NAMES[@]:(-1)}".R1.fastq.gz ]]; then
        	echo "PairT cannot locate trimmed reads files"
        	echo "Please rerun PairT with quality trimming"
        	exit 1
        fi
	
#Next
	echo -e "Now PairT will map all pairs\n"
	for i in "${NAMES[@]}"; do
		if [ -f "$i.R2.fastq.gz" ]; then
			echo "mapping for " $i "individulas: "$i".R1.fastq.gz" "and "$i".R2.fastq.gz"
        		bwa mem -L 20,5 -t $NUMProc -a -M -T 10 -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" reference.fasta $i.R1.fastq.gz $i.R2.fastq.gz  2> bwa.$i.log | mawk '$6 !~/[2-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -@$NUMProc -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
			echo -e "done \n"
        	else
			echo "mapping for " $i "individulas: "$i".R1.fastq.gz only"
        		bwa mem -L 20,5 -t $NUMProc -a -M -T 10 -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" reference.fasta $i.R1.fastq.gz  2> bwa.$i.log | mawk '$6 !~/[2-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -@$NUMProc -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
			echo -e "done \n"
        	fi
        	
		samtools sort -@$NUMProc $i.bam -o $i.bam 2>>$i.bam.log
		mv $i.bam $i-RG.bam
		samtools index $i-RG.bam
        done
  	
fi

##Creating mapping intervals if needed, CreateIntervals function is defined later in script
#If mapping is being performed, intervals are created automatically
if [ "$MAP" != "no" ]; then
	echo -e "\nCreating alignment intervals"
	ls *-RG.bam >bamlist.list
	CreateIntervals 
fi


##SNP Calling Section of code

if [ "$SNP" != "no" ]; then
	#Create list of BAM files
	echo -e "\nCreating alignment intervals"
	ls *-RG.bam >bamlist.list
	#If mapping is not being performed, but intervals do not exist they are created
	if [[ "$MAP" == "no" && ! -f "cat-RRG.bam" ]]; then
		CreateIntervals 
	fi
	#Check for runs from older versions to ensure the recreation of cat-RRG.bam
	if [[ "$MAP" == "no" && -f "map.bed" ]]; then
		CreateIntervals 
	fi
	#Check to make sure interval files have been created
	if [[ "$MAP" == "no" && ! -f "mapped.bed" ]]; then
		bedtools merge -i cat-RRG.bam -bed >  mapped.bed
	fi

	cat namelist | parallel 'bedtools coverage -a mapped.bed -b {}-RG.bam -bed > {}.stast'
	ls *stast | awk '{print substr($0,1,index($0,".")-1)}' >name
	t=$(cat name)
	echo -e "Chrom\tStart\tEnd\t$t" | tr "\n" "\t" |awk '{print $0}' > Bedcov_merge.txt
	temp1=$(head -1 name | awk '{print $0".stast"}')
	for i in $temp1; do awk '{print $1"\t"$2"\t"$3}' $i >col1_3; done
	for i in $t; do awk '{print $4}' $i.stast > $i.col4; done
	paste *col4 >col_4
	paste col1_3 col_4 >> Bedcov_merge.txt
	rm col* | rm col1_3 | rm col_4
else
	exit 1
fi
}

#Function for trimming reads using fastp
trim_reads(){
	
	if [ -f $1.R.fastq.gz ]; then	
		# paired
		fastp --in1 $1.F.fastq.gz --in2 $1.R.fastq.gz --out1 $1.R1.fastq.gz --out2 $1.R2.fastq.gz --cut_by_quality5 20 --cut_by_quality3 20 --cut_window_size 5 --cut_mean_quality 15 --correction $TW -q 15 -u 50 -j $1.json -h $1.html --detect_adapter_for_pe &> $1.trim.log
	fi
}

export -f trim_reads

TrimReads () { 
	#STACKS adds a strange _1 or _2 character to the end of processed reads, this looks for checks for errant characters and replaces them.
	#This functionality is now parallelized and will run if only SE sequences are used.
	NAMES=( `cat "namelist" `)
	STACKS=$(cat namelist| parallel -j $NUMProc --no-notice "gunzip -c {}.F.fastq.gz | head -1" | mawk '$0 !~ /\/1$/ && $0 !~ /\/1[ ,	]/ && $0 !~ / 1:.*[A-Z]*/' | wc -l )
	FB1=$(( $NUMProc / 2 ))
	if [ $STACKS -gt 0 ]; then
		
		echo "Removing the _1 character and replacing with /1 in the name of every sequence"
		cat namelist | parallel -j $FB1 --no-notice "gunzip -c {}.F.fastq.gz | sed -e 's:_1$:/1:g' > {}.F.fastq"
		rm -f *.F.fastq.gz
		cat namelist | parallel -j $FB1 --no-notice "gzip {}.F.fastq"
	fi

	if [ -f "${NAMES[@]:(-1)}".R.fastq.gz ]; then
	
		STACKS=$(cat namelist| parallel -j $NUMProc --no-notice "gunzip -c {}.R.fastq.gz | head -1" | mawk '$0 !~ /\/2$/ && $0 !~ /\/2[ ,	]/ && $0 !~ / 2:.*[A-Z]*/'| wc -l )

		if [ $STACKS -gt 0 ]; then
			echo "Removing the _2 character and replacing with /2 in the name of every sequence"
			cat namelist | parallel -j $FB1 --no-notice "gunzip -c {}.R.fastq.gz | sed -e 's:_2$:/2:g' > {}.R.fastq"
			rm -f *.R.fastq.gz
			cat namelist | parallel -j $FB1 --no-notice "gzip {}.R.fastq"
		fi
	fi

	cat namelist | parallel -j $NUMProc "gunzip -c {}.F.fastq.gz | head -2 | tail -1 >> lengths.txt"
	MLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
    	MLen=$(($MLen / 2))
	TW="--length_required $MLen"	
	mkdir trim_reports &>/dev/null
	cat namelist | parallel --env trim_reads -j $FB1 trim_reads {}
}


##Creating mapping intervals if needed, CreateIntervals function is defined later in script
CreateIntervals()
{
	samtools merge -@$NUMProc -b bamlist.list -f cat-RRG.bam &>/dev/null
	samtools index cat-RRG.bam 
	wait
	bedtools merge -i cat-RRG.bam -bed >  mapped.bed
}


##############  GetInfo function
#This task will get all information for configuration file
GetInfo(){
	#checking number of individuals
	echo "$NumInd individuals are detected. Is this correct? Enter yes or no and press [ENTER]"
	read Indcorrect

	if [ "$Indcorrect" == "no" ]; then
        	echo "Please double check that all fastq files are named PopA_001.F.fastq.gz and PopA_001.R.fastq.gz"
        	exit 1
		elif [ "$Indcorrect" == "yes" ]; then
            		echo "Proceeding with $NumInd individuals"
		else
        		echo "Incorrect Input"
        		exit 1
	fi

	#Tries to get number of processors, if not asks user

	if [[ "$OSTYPE" == "linux"* ]]; then
		NUMProc=( `grep -c ^processor /proc/cpuinfo 2> /dev/null` ) 
	else
		NUMProc=( `sysctl hw.ncpu | cut -f2 -d " " `)
	fi

	NUMProc=$(($NUMProc + 0)) 

	echo "This machine has $NUMProc processors available on this system."
	echo "Please enter the maximum number of processors to use for this analysis."
        read NUMProc
        
	if [ $NUMProc -lt 1 ]; then
        	echo "Incorrect. Please enter the number of processing cores on this computer"
        	read NUMProc
	fi                
	if [ $NUMProc -lt 1 ]; then
        	echo "Incorrect input, exiting"
        	exit 1
	fi

	#Tries to get maximum system memory, if not asks user
	if [[ "$OSTYPE" == "linux"* ]]; then
		MAXMemory=$(($(grep -Po '(?<=^MemTotal:)\s*[0-9]+' /proc/meminfo | tr -d " ") / 1000000))
		echo "Your machine has $MAXMemory gigabytes of maximum memory available on this system."
		echo "Please enter the maximum memory to use for this analysis in gigabytes"
		echo "For example, to limit to ten gigabytes, enter 10"
	
		read MAXMemory1
		MAXMemory1=$( echo $MAXMemory1 | sed 's/[g,G]//g' )
		MAXMemory=$(( $MAXMemory1 / $NUMProc ))G

		while [[ -z $MAXMemory ]];
		do
			echo "Incorrect input"
			echo -e "Please enter the maximum memory to use for this analysis in gigabytes." 
			read MAXMemory1
			MAXMemory=$(( $MAXMemory1 / $NUMProc ))G
		done
	else
			MAXMemory=0
	fi

	#Asks if user wants to trim reads.  This allows this part of the pipeline to be skipped during subsequent analyses
	echo -e "\nDo you want to quality trim your reads?" 
	echo "Type yes or no and press [ENTER]?"

	read TRIM

	#Asks if user wants to map reads and change default mapping variables for BWA
	echo "Do you want to map reads?  Type yes or no and press [ENTER]"
	read MAP
	
	if [ "$MAP" == "no" ]; then
        	echo "Mapping will not be performed"
        	optA=1
    		optB=4
    		optO=6
        else
                echo "BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa."
                echo "Would you like to enter a new parameters now? Type yes or no and press [ENTER]"
                read optq

        	if [ "$optq" == "yes" ]; then
        		echo "Please enter new value for A (match score).  It should be an integer.  Default is 1."
        		read newA
        		optA=$newA
                	echo "Please enter new value for B (mismatch score).  It should be an integer.  Default is 4."
        		read newB
        		optB=$newB
                	echo "Please enter new value for O (gap penalty).  It should be an integer.  Default is 6."
        		read newO
        		optO=$newO
        	else
                	echo "Proceeding with default values for BWA read mapping."
                	optA=1
               		optB=4
                	optO=6
        	fi
	fi

	#Does user wish to call bedtools coverage?
	echo "Do you want to use bedtools to estimate coverage?  Please type yes or no and press [ENTER]"
	read SNP

	while [[ $SNP != "yes" && $SNP != "no" ]];
	do
		echo "Incorrect input"
		echo -e "Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]"
		read SNP
	done

}


#Execute directly main here
if [ -n "$1" ]; then
	main $1 2>&1 | tee temp.LOG #Log all output
else
	main 2>&1 | tee temp.LOG #Log all output
fi

mawk '!/#.*%/' temp.LOG >> PairT.LOG
rm temp.LOG


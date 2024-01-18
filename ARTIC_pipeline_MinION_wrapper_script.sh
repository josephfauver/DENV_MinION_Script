################################################################################################
###### NOTE: This script is intended to be run within an activated ARTIC conda environment #####
################################################################################################


################################################################################################
######################### SPECIFY VARIABLES NEEDED TO RUN THE PIPELINE #########################
################################################################################################

## provide path to output directory:
OUTDIR=/home/fauverlab/artic/fieldbioinformatics/analysis/DENV_ARTIC_test

## provide path to scheme directory:
SCHEME_DIR=/home/fauverlab/artic/fieldbioinformatics/scheme

## provide path to directory in which basecalled fastqs are located:
READS_DIR=/home/fauverlab/artic/fieldbioinformatics/pass_all_DENV_test_reads
	## this should be a directory containing one subdirectory each per sample, named by barcode (e.g. "barcode01")
	## each subdirectory should contain basecalled fastq.gz files for that sample, as output in a guppy "pass" directory

## provide comma separated list linking barcode names and sample names
LIST=/home/fauverlab/artic/fieldbioinformatics/analysis/DENV_ARTIC_test/samples.txt
	## for example:
		## barcode01,samplename01
		## barcode02,samplename02
		## barcode03,samplename03

## provide mininum and maximum read lengths for guppyplex to keep
MIN_READLEN=400
MAX_READLEN=700

## assign the appropriate medaka model based on sequencing and basecalling parameters
MEDAKA_MODEL=r1041_e82_400bps_sup_v4.2.0

################################################################################################
######################### SHOULD NOT NEED TO EDIT CODE BELOW THIS LINE #########################
################################################################################################



######################################################
#### reports how many samples are being processed ####
######################################################

echo ""
NUM_SAMPLES=$(cat ${LIST} | wc -l)
echo "***********************************************"
echo "** TOTAL NUMBER OF SAMPLES TO PROCESS: " ${NUM_SAMPLES}
echo "***********************************************"
echo ""

#################################################################
#### runs the pipeline, looping over each line in input list ####
#################################################################

while IFS="," read -r BARCODE SAMPLE
do
	echo "***********************************************"
	echo "** PROCESSING" ${SAMPLE} "with" ${BARCODE}"..."
	echo "***********************************************"
	echo ""

	cd ${READS_DIR}/${BARCODE}
	LOGIC_TEST=$(ls | grep -c "fastq.gz")
	if [[ ${LOGIC_TEST} -gt 0 ]]
	then
		echo "***********************************************"
		echo "** UNZIPPING FASTQS FOR" ${SAMPLE}"..."
		echo "***********************************************"
		echo ""
		gunzip *.fastq.gz
	else
		echo "***********************************************"
		echo "** ALL FASTQS ALREADY UNZIPPED"
		echo "***********************************************"
		echo ""
	fi

	echo "***********************************************"
	echo "** DETERMINING SEROTYPE FOR" ${SAMPLE}"..."
	echo "***********************************************"
	echo ""
	FASTQS=(*.fastq)
	FASTQ1=(${FASTQS[0]})
	FASTQ2=(${FASTQS[1]})
	FASTQ3=(${FASTQS[2]})
	FASTQ4=(${FASTQS[3]})
	FASTQ5=(${FASTQS[4]})
	FASTQ6=(${FASTQS[5]})
	FASTQ7=(${FASTQS[6]})
	FASTQ8=(${FASTQS[7]})
	FASTQ9=(${FASTQS[8]})
	FASTQ10=(${FASTQS[9]})
	FASTQ11=(${FASTQS[10]})
	FASTQ12=(${FASTQS[11]})
	FASTQ13=(${FASTQS[12]})
	FASTQ14=(${FASTQS[13]})
	FASTQ15=(${FASTQS[15]})
	DENV1_REF=${SCHEME_DIR}/DENV1/V1/DENV1.reference.fasta
	DENV2_REF=${SCHEME_DIR}/DENV2/V1/DENV2.reference.fasta
	DENV3_REF=${SCHEME_DIR}/DENV3/V1/DENV3.reference.fasta
	DENV4_REF=${SCHEME_DIR}/DENV4/V1/DENV4.reference.fasta
	minimap2 -ax asm20 ${DENV1_REF} ${FASTQ1} ${FASTQ2} ${FASTQ3} ${FASTQ4} ${FASTQ5} ${FASTQ6} ${FASTQ7} ${FASTQ8} ${FASTQ9} ${FASTQ10} ${FASTQ11} ${FASTQ12} ${FASTQ13} ${FASTQ14} ${FASTQ15} | samtools sort -o ${SAMPLE}_DENV1_sorted.bam -T ${SAMPLE}_DENV1.tmp
	minimap2 -ax asm20 ${DENV2_REF} ${FASTQ1} ${FASTQ2} ${FASTQ3} ${FASTQ4} ${FASTQ5} ${FASTQ6} ${FASTQ7} ${FASTQ8} ${FASTQ9} ${FASTQ10} ${FASTQ11} ${FASTQ12} ${FASTQ13} ${FASTQ14} ${FASTQ15} | samtools sort -o ${SAMPLE}_DENV2_sorted.bam -T ${SAMPLE}_DENV2.tmp
	minimap2 -ax asm20 ${DENV3_REF} ${FASTQ1} ${FASTQ2} ${FASTQ3} ${FASTQ4} ${FASTQ5} ${FASTQ6} ${FASTQ7} ${FASTQ8} ${FASTQ9} ${FASTQ10} ${FASTQ11} ${FASTQ12} ${FASTQ13} ${FASTQ14} ${FASTQ15} | samtools sort -o ${SAMPLE}_DENV3_sorted.bam -T ${SAMPLE}_DENV3.tmp
	minimap2 -ax asm20 ${DENV4_REF} ${FASTQ1} ${FASTQ2} ${FASTQ3} ${FASTQ4} ${FASTQ5} ${FASTQ6} ${FASTQ7} ${FASTQ8} ${FASTQ9} ${FASTQ10} ${FASTQ11} ${FASTQ12} ${FASTQ13} ${FASTQ14} ${FASTQ15} | samtools sort -o ${SAMPLE}_DENV4_sorted.bam -T ${SAMPLE}_DENV4.tmp
	DENV1=$(samtools view -c -F 260 ${SAMPLE}_DENV1_sorted.bam)
	echo ""
	echo "***********************************************"
	echo "** NO. READS FOR" ${SAMPLE}" MAPPED TO DENV1:" ${DENV1}
	echo "***********************************************"
	echo ""
	DENV2=$(samtools view -c -F 260 ${SAMPLE}_DENV2_sorted.bam)
	echo "***********************************************"
	echo "** NO. READS FOR" ${SAMPLE}" MAPPED TO DENV2:" ${DENV2}
	echo "***********************************************"
	echo ""
	DENV3=$(samtools view -c -F 260 ${SAMPLE}_DENV3_sorted.bam)
	echo "***********************************************"
	echo "** NO. READS FOR" ${SAMPLE}" MAPPED TO DENV3:" ${DENV3}
	echo "***********************************************"
	echo ""
	DENV4=$(samtools view -c -F 260 ${SAMPLE}_DENV4_sorted.bam)
	echo "***********************************************"
	echo "** NO. READS FOR" ${SAMPLE}" MAPPED TO DENV4:" ${DENV4}
	echo "***********************************************"
	echo ""

	maxvarname() {
				for i; do
					echo "${!i} $i"
			done | sort -nr | sed -n '1s/.* \(.*\)/\1/p'
			}
			SEROTYPE=$(maxvarname DENV1 DENV2 DENV3 DENV4)

	if (( ${SEROTYPE} == 0 ))
	then
		echo "********************************************************"
		echo "** NO READS MAPPED; ASSINGING SEROTYPE DENV1 BY DEFAULT."
		echo "********************************************************"
		echo ""
		SEROTYPE=DENV1
	else
		
		echo "      ***********************************"
		echo "      **" ${SAMPLE} "IS SEROTYPE" ${SEROTYPE}
		echo "      ***********************************"
		echo ""
	fi

	rm *.bam

	echo "***********************************************"
	echo "** RUNNING GUPPYPLEX FOR" ${SAMPLE}"..."
	echo "***********************************************"
	echo ""
	artic guppyplex \
		--skip-quality-check \
		--min-length ${MIN_READLEN} \
		--max-length ${MAX_READLEN} \
		--directory ${READS_DIR}/${BARCODE} \
		--prefix all_${SAMPLE}

	echo ""
	echo "***********************************************"
	echo "** STARTING ARTIC PIPELINE FOR" ${SAMPLE}"..."
	echo "***********************************************"
	echo ""
	cd ${OUTDIR}
	mkdir ${BARCODE}_${SAMPLE}
	cd ${BARCODE}_${SAMPLE}
	artic minion \
		--scheme-directory ${SCHEME_DIR} \
		--read-file ${READS_DIR}/${BARCODE}/all_${SAMPLE}_${BARCODE}.fastq \
		--normalise 400 \
		--medaka \
		--medaka-model ${MEDAKA_MODEL} \
		${SEROTYPE}/V1 ${SAMPLE}

	echo ""
	echo "***********************************************************************"
	echo "** ANALYSIS FOR SAMPLE" ${SAMPLE} "FINISHED. CHECK ARTIC OUTPUT FOR ERRORS."
	echo "***********************************************************************"
	echo ""
	echo ""
done < ${LIST}

echo "ALL ANALYSES COMPLETE. GOODBYE."
echo""


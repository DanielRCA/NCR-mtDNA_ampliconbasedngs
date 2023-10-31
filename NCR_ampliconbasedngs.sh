#### Script generated to analyse data from NCR-mtDNA for ancient samples obtained by amplicon-based NGS methods ####

## time bash NCR_ampliconbasedngs.sh 2>&1 | tee NCR_ampliconbasedngs.log


source activate env_ncrngsamplicons


main_dir=$(pwd)

#A directory to save the files generated is made
vis=visualization
if [[ -d ${vis} ]]; then
    echo "Directory ${vis} already exists"
else
    mkdir ${vis}
fi

#A directory to save the reference is made
ref=reference
seq_ref=rCRS_NCR_lineal.fasta


if [[ -d ${ref} ]]; then
    echo "Directory ${ref} already exist, it is not going to be created again"
    NCR=${main_dir}/${ref}/${seq_ref}
else
    echo "A directory is created to keep the reference sequence and its indexes"
    mkdir ${ref}
    mv ${seq_ref} ${ref}
fi
cd ${ref}
#Indexes would be generated only the first time.
if [[ -e ${seq_ref}.amb ]] && [[ -e ${seq_ref}.fai ]] ; then
    echo "
    Indexes were already created
    "
else
    echo "
    Indexes are being created
    "
    time bwa index ${seq_ref}
    time samtools faidx ${seq_ref}
    picard CreateSequenceDictionary R=${seq_ref} O=rCRS_NCR_lineal.dict
fi
cd ${main_dir}

NCR=${ref}/${seq_ref}



cd ${main_dir}

NCR=${main_dir}/${ref}/${seq_ref}

#A list with all the directories that have samples is made
for samples in *; do if  [[ -d ${samples} ]] && [[ ${samples} != ${vis} ]] && [[ ${samples} != ${ref} ]]; then  echo ${samples}; fi done > samples.txt

echo -e "Sample\tPosition\tAlternative_base\tN_total\tN_alternative\tPercentage" > base_mix_cov_5.txt
echo -e "Sample\tPosition\tAlternative_base\tN_total\tN_alternative\tPercentage" > base_mix_cov_10.txt

while read sample; do
      
    echo "
    Starting with sample ${sample}
    "
    cd "${sample}"
        
    sample_dir=$(pwd)
   
    #Quality is checked with a FastQC
    fastq_file=1_fastqc_file_quality_pre-trimming
    if [[ -d ${fastq_file} ]]; then
        echo "Directory ${fastq_file} already exists"
    else
        mkdir ${fastq_file}
        fastqc -o ${fastq_file} -f fastq *fastq.gz
    fi
   
   #First fastp. To trim adapters and remove duplicates
   fastp --dedup \
       --detect_adapter_for_pe \
       --trim_poly_x \
       --length_required 30 \
       --thread 8 \
       --json ${sample}_temp.json \
       --out1 ${sample}_fastp_filtered_R1_temp.fq.gz \
       --out2 ${sample}_fastp_filtered_R2_temp.fq.gz \
       --in1 *R1*.fastq.gz \
       --in2 *R2*.fastq.gz
        
    #Number reads from fastqs
    num_1=$(grep before_filtering -m 1 -A1 ${sample}_temp.json | grep total | sed "s/.*://g" | sed "s/,//g")
    echo "${sample} ${num_1}" >> ${main_dir}/tem_all.txt
    
    #Percentage of duplication
    dup_rate=$(echo "$(grep duplication -m 1 -A1 ${sample}_temp.json | grep rate | sed "s/.*://g" | sed "s/,//g")*100" | bc)    
    echo "${sample} ${dup_rate}" >> ${main_dir}/tem_dup_rate.txt

    
    #Second fastp. To trim primers (25pb from 5') and 5pb in 3'
	fastp --trim_front1 25 --trim_tail1 5 \
		--trim_front2 25 --trim_tail2 5 \
		--out1 ${sample}_nonprimer_R1_temp.fq.gz \
		--out2 ${sample}_nonprimer_R2_temp.fq.gz \
		--in1 ${sample}_fastp_filtered_R1_temp.fq.gz \
		--in2 ${sample}_fastp_filtered_R2_temp.fq.gz
		
	#Third fastp. To trim bad quality bases
	fastp --cut_tail --cut_tail_mean_quality 20 \
		--trim_poly_x \
		--merge --overlap_len_require 10 \
		--correction \
		--html ${sample}_temp.html \
		--json ${sample}_temp.json \
		--merged_out ${sample}_M.fq.gz \
		--out1 ${sample}_R1.fq.gz \
		--out2 ${sample}_R2.fq.gz \
		--in1 ${sample}_nonprimer_R1_temp.fq.gz \
		--in2 ${sample}_nonprimer_R2_temp.fq.gz
    rm *_temp*

    #Quality is checked with FastQC
    fastq_file=2_fastqc_file_quality_post-trimming
    if [[ -d ${fastq_file} ]]; then
        echo "Directory ${fastq_file} already exists"
    else
        mkdir ${fastq_file}
        fastqc -o ${fastq_file} -f fastq *.fq.gz
    fi
    
    
    #Alignment is done with bwa aln and bwa samse/sampe   
    bwa aln -t 9 -l 1024 -n 0.01 -o 2 ${NCR} ${sample}_R1.fq.gz > tem_aln_R1.sai
    bwa aln -t 9 -l 1024 -n 0.01 -o 2 ${NCR} ${sample}_R2.fq.gz > tem_aln_R2.sai
    bwa aln -t 9 -l 1024 -n 0.01 -o 2 ${NCR} ${sample}_M.fq.gz > tem_aln_M.sai
    
    
    bwa sampe ${NCR} tem_aln_R1.sai tem_aln_R2.sai ${sample}_R1.fq.gz ${sample}_R2.fq.gz | samtools view -Sbh -F4 -q30 | samtools sort -O bam -o tem_${sample}_NM_map.bam
    bwa samse ${NCR} tem_aln_M.sai ${sample}_M.fq.gz | samtools view -Sbh -F4 -q30 | samtools sort -O bam -o tem_${sample}_M_map.bam

    samtools merge -o tem_${sample}_map_unsort.bam tem_${sample}_NM_map.bam tem_${sample}_M_map.bam; samtools sort -O bam -o final_${sample}.bam tem_${sample}_map_unsort.bam
    bam=final_${sample}.bam
    samtools index ${bam}
        
    #The number of reads with good quality is determined
    num_2=$(samtools view -c ${bam})
    echo "${sample} ${num_2}" >> ${main_dir}/tem_useful.txt

    #Reads with post-mortem molecular (PMD) damage are extracted
    samtools view -h ${bam} | pmdtools --threshold 1 --header | samtools view -Sb - > pmd_final_${sample}.bam    
    pam=pmd_final_${sample}.bam
    samtools index ${pam}
        
    #The number of reads with PMD is determined
    num_3=$(samtools view -c ${pam})
    echo "${sample} ${num_3}" >> ${main_dir}/tem_useful_pmd.txt

    #Quality and mean depth coverage of bams is checked for the file with all reads
	qualimap bamqc -bam ${bam} -c -gd hg19 -outdir 3_qualimap_bamqc
	qual=$(grep "mean mapping quality" 3_qualimap_bamqc/genome_results.txt | sed 's/mean mapping quality = //')
	echo "${sample} ${qual}" >> ${main_dir}/tem_qual.txt
	cov=$(grep "mean coverageData" 3_qualimap_bamqc/genome_results.txt | sed 's/mean coverageData = //' | sed 's/X//' | sed 's/,//')
	echo "${sample} ${cov}" >> ${main_dir}/tem_cov.txt

    #Ranges for positions with a depth coverage over 5X and 10X are generated
    echo "
    Ranges of positions are being determined...
    "
	#Depth coverage over 5X in file with all reads
    samtools depth ${bam} | awk '$3>4 {print $2}' |  awk '{if ($1<670) print $1+15900; else print $1-669}' | awk  '$1 < 303 || $1 > 315 {print $1}' | sort -n > temp.txt
	range=$(python ${main_dir}/range.py temp.txt | sed "s/[',\[\)\( ]//g" | sed "s/]//g" | sed "s/;$//")
	#Depth coverage over 10X in file with all reads
	echo "${sample}_5 ${range}" >> ${main_dir}/tem_range_5.txt
	samtools depth ${bam} | awk '$3>9 {print $2}' |  awk '{if ($1<670) print $1+15900; else print $1-669}' | awk  '$1 < 303 || $1 > 315 {print $1}' | sort -n > temp.txt
	range=$(python ${main_dir}/range.py temp.txt | sed "s/[',\[\)\( ]//g" | sed "s/]//g" | sed "s/;$//")
 	#Depth coverage over 5X in file with damaged reads
	echo "${sample}_10 ${range}" >> ${main_dir}/tem_range_10.txt
    samtools depth ${pam} | awk -F "\t" '$3>4 {print($2)}' |  awk '{if ($1<670) print $1+15900; else print $1-669}' | awk  '$1 < 303 || $1 > 315 {print $1}' | sort -n > temp_pmd.txt
	range_pmd=$(python ${main_dir}/range.py temp_pmd.txt | sed "s/[',\[\)\( ]//g" | sed "s/]//g" | sed "s/;$//")
 	#Depth coverage over 10X in file with daamged reads
	echo "${sample}_5_pmd ${range_pmd}" >> ${main_dir}/tem_range_pmd_5.txt  
	samtools depth ${pam} | awk -F "\t" '$3>9 {print($2)}' |  awk '{if ($1<670) print $1+15900; else print $1-669}' | awk  '$1 < 303 || $1 > 315 {print $1}' | sort -n > temp_pmd.txt
	range_pmd=$(python ${main_dir}/range.py temp_pmd.txt | sed "s/[',\[\)\( ]//g" | sed "s/]//g" | sed "s/;$//")
	echo "${sample}_10_pmd ${range_pmd}" >> ${main_dir}/tem_range_pmd_10.txt
    rm temp_pmd.txt
    
    bases_with_good_cov=$(samtools depth ${bam} | awk '$3>4 {print $2"\t"1}' | awk '{sum += $2} END {print sum}')
    echo -e "${sample}\t${bases_with_good_cov}" >> ${main_dir}/tem_goodcovinsamples_5.txt
    bases_with_good_cov=$(samtools depth ${bam} | awk '$3>9 {print $2"\t"1}' | awk '{sum += $2} END {print sum}')
    echo -e "${sample}\t${bases_with_good_cov}" >> ${main_dir}/tem_goodcovinsamples_10.txt

    #Variant call is done with freebayes and vcfallelicprimitives
    echo "
    Starting the variant call
    "
	#All reads, depth coverage 5X
	freebayes -f ${NCR} -i -X -F 0.3 -C 2 --min-coverage 5 ${bam} | vcfallelicprimitives -kg > ${bam}_5.vcf 
	awk '$3 == "." && length($5) == 1 && length($4) == 1 {print $2"\t"$4"\t"$5"\t"$10}' ${bam}_5.vcf | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"100*$9/$5}' | awk '{if ($6 < 70) print $1"\t"$2$3; else print $1"\t"$3;}' |  sed 's/,.*//' | sed 's/[AGTC]\{4\}/N/g' | sed 's/[ACG]\{3\}/V/g' |sed 's/[ATC]\{3\}/H/g' | sed 's/[ATG]\{3\}/D/g' | sed 's/[GTC]\{3\}/B/g' | sed 's/[AC]\{2\}/M/g' | sed 's/[TG]\{2\}/K/g' | sed 's/[AT]\{2\}/W/g' | sed 's/[GC]\{2\}/S/g' | sed 's/[TC]\{2\}/Y/g' | sed 's/[AG]\{2\}/R/g' | awk '{if ($1<670) print $1+15900$2; else print $1-669$2}' |   sed -e 's/\t$//g' | awk -F"\t" '$1 < 303 || $1 > 315 {print $1$2}' > ${bam}_5.txt
	#Damaged reads, depth coverage 5X
	freebayes -f ${NCR} -i -X -F 0.3 -C 2 --min-coverage 5 ${pam} | vcfallelicprimitives -kg > ${pam}_5.vcf
	awk '$3 == "." && length($5) == 1 && length($4) == 1 {print $2"\t"$4"\t"$5"\t"$10}' ${pam}_5.vcf | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"100*$9/$5}' | awk '{if ($6 < 70) print $1"\t"$2$3; else print $1"\t"$3;}' | sed 's/[AGTC]\{4\}/N/g' | sed 's/[ACG]\{3\}/V/g' |sed 's/[ATC]\{3\}/H/g' | sed 's/[ATG]\{3\}/D/g' | sed 's/[GTC]\{3\}/B/g' | sed 's/[AC]\{2\}/M/g' | sed 's/[TG]\{2\}/K/g' | sed 's/[AT]\{2\}/W/g' | sed 's/[GC]\{2\}/S/g' | sed 's/[TC]\{2\}/Y/g' | sed 's/[AG]\{2\}/R/g' | awk '{if ($1<670) print $1+15900$2; else print $1-669$2}' | sed -e 's/\t$//g' | awk -F"\t" '$1 < 303 || $1 > 315 {print $1$2}' > ${pam}_5.txt
	#All reads, depth coverage 10X
	freebayes -f ${NCR} -i -X -F 0.3 -C 3 --min-coverage 10 ${bam} | vcfallelicprimitives -kg > ${bam}_10.vcf 
	awk '$3 == "." && length($5) == 1 && length($4) == 1 {print $2"\t"$4"\t"$5"\t"$10}' ${bam}_10.vcf | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"100*$9/$5}' | awk '{if ($6 < 70) print $1"\t"$2$3; else print $1"\t"$3;}' |  sed 's/,.*//' | sed 's/[AGTC]\{4\}/N/g' | sed 's/[ACG]\{3\}/V/g' |sed 's/[ATC]\{3\}/H/g' | sed 's/[ATG]\{3\}/D/g' | sed 's/[GTC]\{3\}/B/g' | sed 's/[AC]\{2\}/M/g' | sed 's/[TG]\{2\}/K/g' | sed 's/[AT]\{2\}/W/g' | sed 's/[GC]\{2\}/S/g' | sed 's/[TC]\{2\}/Y/g' | sed 's/[AG]\{2\}/R/g' | awk '{if ($1<670) print $1+15900$2; else print $1-669$2}' |   sed -e 's/\t$//g' | awk -F"\t" '$1 < 303 || $1 > 315 {print $1$2}' > ${bam}_10.txt
	#Damaged reads, depth coverage 10X
	freebayes -f ${NCR} -i -X -F 0.3 -C 3 --min-coverage 10 ${pam} | vcfallelicprimitives -kg > ${pam}_10.vcf
	awk '$3 == "." && length($5) == 1 && length($4) == 1 {print $2"\t"$4"\t"$5"\t"$10}' ${pam}_10.vcf | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"100*$9/$5}' | awk '{if ($6 < 70) print $1"\t"$2$3; else print $1"\t"$3;}' | sed 's/[AGTC]\{4\}/N/g' | sed 's/[ACG]\{3\}/V/g' |sed 's/[ATC]\{3\}/H/g' | sed 's/[ATG]\{3\}/D/g' | sed 's/[GTC]\{3\}/B/g' | sed 's/[AC]\{2\}/M/g' | sed 's/[TG]\{2\}/K/g' | sed 's/[AT]\{2\}/W/g' | sed 's/[GC]\{2\}/S/g' | sed 's/[TC]\{2\}/Y/g' | sed 's/[AG]\{2\}/R/g' | awk '{if ($1<670) print $1+15900$2; else print $1-669$2}' | sed -e 's/\t$//g' | awk -F"\t" '$1 < 303 || $1 > 315 {print $1$2}' > ${pam}_10.txt
		
	#Documents with mixtures are done for all reads
 	awk '$3 == "." && length($4) == 1 && length($5) == 1 {if ($2<670) print $2+15900"\t"$5"\t"$10; else print $2-669"\t"$5"\t"$10}' ${bam}_5.vcf | awk  '$1 < 303 || $1 > 315 {print $0}' | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$4"\t"$8"\t"100*$8/$4"\t"}' | awk -v vd=${sample} '$5 < 70 {print vd"\t"$0}' >> ${main_dir}/base_mix_cov_5.txt
	
	awk '$3 == "." && length($4) == 1 && length($5) == 1 {if ($2<670) print $2+15900"\t"$5"\t"$10; else print $2-669"\t"$5"\t"$10}' ${bam}_10.vcf | awk  '$1 < 303 || $1 > 315 {print $0}' | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$4"\t"$8"\t"100*$8/$4"\t"}' | awk -v vd=${sample} '$5 < 70 {print vd"\t"$0}' >> ${main_dir}/base_mix_cov_10.txt
 
    
    cp ${bam}_*.txt ${main_dir}/${vis}
    cp ${pam}_*.txt ${main_dir}/${vis}
    
    rm *tem*
    echo "
    Changing sample...
    " 
    cd ${main_dir}

done < samples.txt

#File in hsd format (the one requiered by HaploGrep2) is generated in "Visualization directory"
cd ${main_dir}/${vis}
echo "
Creating the hsd file to use in HaploGrep2
"
echo "Sample Range Haplogroup Haplotype" > tem_haplogrep_format.txt
#Variable with the same value for all samples is generated
haplogroup="?"
#The rest of the variables are generated line by line
for t in *final*bam_*.txt; do
    #Name is extracted
    sample="${t//*final_/}"
    if [[ ${t} == pmd* ]]
    then
        if [[ ${t} == *_5.txt ]]; then
			sample="${sample//.bam_5.*/_5_pmd}"
			ranges=$(awk -v SAMPLE=${sample} '$1==SAMPLE {print $2}' ${main_dir}/tem_range_pmd_5.txt)
		else
			sample="${sample//.bam_10.*/_10_pmd}"
			ranges=$(awk -v SAMPLE=${sample} '$1==SAMPLE {print $2}' ${main_dir}/tem_range_pmd_10.txt)
		fi
	else
		if [[ ${t} == *_5.txt ]]; then
			sample="${sample//.bam_5.*/_5}"
			ranges=$(awk -v SAMPLE=${sample} '$1==SAMPLE {print $2}' ${main_dir}/tem_range_5.txt)
		else
			sample="${sample//.bam_10.*/_10}"
			ranges=$(awk -v SAMPLE=${sample} '$1==SAMPLE {print $2}' ${main_dir}/tem_range_10.txt)
		fi
    fi
    #Fake tabulators are converted into real ones
    sed -i 's/    /\t/g' ${t}
    #Polimoprhisms are print in a temporal file (vertical)
    awk -F"\t" '{print $1 $2}' ${t} > temp_file.txt
    #Colum is transformed into a row and save in an object
    haplotype=$(for i in `< temp_file.txt`; do echo -en ${i}" "; done)
    #Row is added in the order determined by HaploGrep2 (with spaces instead of tabulators)
    echo "${sample} ${ranges} ${haplogroup} ${haplotype}" >> tem_haplogrep_format.txt
done

#Spaces are transformed into tabulators
sed 's/ /\t/g' tem_haplogrep_format.txt | sed 's/\t$//g' | awk -F "\t" '{ if ($4 == "") print $0"\t16519T"; else print $0}' > haplogrep_format.txt
cp haplogrep_format.txt ${main_dir}
rm *tem*

cd ${main_dir}

echo "
Starting HaploGrep2...
"
grep -v pmd haplogrep_format.txt > tem_haplogrep_format.txt
bash haplogrep classify --in tem_haplogrep_format.txt --format hsd --out tem_haplos.txt

echo "
Generating final table...
"
echo -e 'ID\tTotal_reads\tDuplication_rate\tUtil_reads\tPercentage_of_util_reads\tMean_depth_coverage\tMean_mapping_quality\tDamaged_util_reads\tPercentage_of_Damaged_reads\tMixed_bases_cov5\tHaplogrouop_cov5\tQuality_cov5\tPercentage_of_NCR_recovered_cov5\tRange_cov5\tHaplotype_cov5\tMixed_bases_cov10\tHaplogrouop_cov10\tQuality_cov10\tPercentage_of_NCR_recovered_cov10\tRange_cov10\tHaplotype_cov10' > final_table.txt

while read sample; do
	echo "${sample}" | awk -v v0=${sample} \
	-v v1=$(awk -v v1_1=${sample} '$1==v1_1 {print $2}' tem_all.txt) \
	-v v2=$(awk -v v2_1=${sample} '$1==v2_1 {print $2}' tem_dup_rate.txt) \
	-v v3=$(awk -v v3_1=${sample} '$1==v3_1 {print $2}' tem_useful.txt) \
	-v v4=$(awk -v v4_1=${sample} '$1==v4_1 {print $2}' tem_cov.txt) \
	-v v5=$(awk -v v5_1=${sample} '$1==v5_1 {print $2}' tem_qual.txt) \
	-v v6=$(awk -v v6_1=${sample} '$1==v6_1 {print $2}' tem_useful_pmd.txt) \
	-v v7=$(awk -v v7_1=${sample} '$1==v7_1 {print $0}' base_mix_cov_5.txt | wc -l) \
	-v v8=$(awk -v v8_1=\"${sample}_5\" '$1==v8_1 {print $2"\t"$4}' tem_haplos.txt | sed 's/"//g' | sed 's/\t/--/g') \
	-v v9=$(awk -v v9_1=${sample} '$1==v9_1 {print $2}' tem_goodcovinsamples_5.txt) \
	-v v10=$(awk -v v10_1=${sample}_5 '$1==v10_1 {print $2}' tem_range_5.txt) \
	-v v11=$(awk -v v11_1=${sample}_5 '$1==v11_1 {print $0}' tem_haplogrep_format.txt | cut -f 4- | sed 's/\t/__/g') \
	-v v12=$(awk -v v12_1=${sample} '$1==v12_1 {print $0}' base_mix_cov_10.txt | wc -l) \
	-v v13=$(awk -v v13_1=\"${sample}_10\" '$1==v13_1 {print $2"\t"$4}' tem_haplos.txt | sed 's/"//g' | sed 's/\t/--/g') \
	-v v14=$(awk -v v14_1=${sample} '$1==v14_1 {print $2}' tem_goodcovinsamples_10.txt) \
	-v v15=$(awk -v v15_1=${sample}_10 '$1==v15_1 {print $2}' tem_range_10.txt) \
	-v v16=$(awk -v v16_1=${sample}_10 '$1==v16_1 {print $0}' tem_haplogrep_format.txt | cut -f 4- | sed 's/\t/__/g') \
	'$1==v0 {print v0"\t"v1"\t"v2"\t"v3"\t"v3/v1*100"\t"v4"\t"v5"\t"v6"\t"v6/v3*100"\t"v7"\t"v8"\t"100*v9/1198"\t"v10"\t"v11"\t"v12"\t"v13"\t"100*v14/1198"\t"v15"\t"v16}' \
	| sed 's/--/\t/g' | sed 's/__/ /g' |  sed 's/16519T$//g' | sed 's/\t16519T\t/\t\t/g' |sed 's/,/\./g' | sort -h >> final_table.txt
done < samples.txt

rm *tem*.txt
rm -r ${vis}
            
echo "          
                THE END
     "

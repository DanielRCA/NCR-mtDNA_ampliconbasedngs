#### Script generated to analyse data from NCR-mtDNA for ancient samples obtained by amplicon-based NGS methods ####

## time bash NCR_ampliconbasedngs.sh 2>&1 | tee NCR_ampliconbasedngs.log


source activate ngs_analysis
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
seq_ref=rCRS_NCR_linear.fasta

if [[ -d ${ref} ]]; then
    echo "Directory ${ref} already exist, it is not going to be created again"
    NCR=$(pwd)/${ref}/${seq_ref}
else
    echo "A directory is created to keep the reference sequence and its indexes"
    mkdir ${ref}
    mv ${seq_ref} ${ref}
fi
cd ${ref}
#Si ya existen los índices, se enviará un mensaje, en caso contrartio, se crearán.
if [[ -e ${seq_ref}.amb ]] && [[ -e ${seq_ref}.fai ]] ; then
    echo "
    Indexes were already created
    "
    picard CreateSequenceDictionary R=${seq_ref} O=rCRS_NCR_linear.dict
    NCR=$(pwd)/${re2}
else
    echo "
    Indexes are being created
    "
    time bwa index ${seq_ref}
    time samtools faidx ${seq_ref}
    picard CreateSequenceDictionary R=${seq_ref} O=rCRS_NCR_linear.dict
fi
cd ${main_dir}

NCR=${ref}/${seq_ref}

#A list with all the directories that have samples is made
for sample in *; do if  [[ -d ${d} ]] && [[ ${d} != ${vis} ]] && [[ ${d} != ${ref} ]]; then  echo $sample; fi done > temp_samples.txt

echo -e "Sample\tPosition\tAlternative_base\tN_total\tN_alternative\tPercentage" > base_mix.txt

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
   
   #First fastp: to remove duplicates
   fastp --verbose \
         --dedup \
         --disable_length_filtering \
         --disable_quality_filtering \
         --json ${sample}_dedup_temp.json \
         --thread 8 \
         --out1 ${sample}_fastp_dedup_R1_temp.fq.gz \
         --out2 ${sample}_fastp_dedup_R2_temp.fq.gz \
         --in1 *R1*.fastq.gz \
         --in2 *R2*.fastq.gz
        
    #Number of reads with duplicates is determined
    num_1=$(grep before_filtering -m 1 -A1 ${sample}_dedup_temp.json | grep total | sed "s/.*://g" | sed "s/,//g")
    echo "Number of total reads ${sample}: ${num_1}" >> reads.txt
    echo "${sample} ${num_1}" >> ${main_dir}/tem_all.txt
    
    #Number of reads without duplicates is determined
    num_2=$(grep after_filtering -m 1 -A1 ${sample}_dedup_temp.json | grep total | sed "s/.*://g" | sed "s/,//g")
    echo "Number of total reads without duplicates ${sample}: ${num_2}" >> reads.txt
    echo "${sample} ${num_2}" >> ${main_dir}/tem_no_dups.txt
    
    #Second fastp: to trim adapters
    fastp --detect_adapter_for_pe \
          --trim_poly_x \
          --length_required 30 \
          --thread 8 \
          --out1 ${sample}_fastp_filtered_R1_temp.fq.gz \
          --out2 ${sample}_fastp_filtered_R2_temp.fq.gz \
          --in1 ${sample}_fastp_dedup_R1_temp.fq.gz \
          --in2 ${sample}_fastp_dedup_R2_temp.fq.gz
        
    #Third fastp: to trim primers (25pb from 5') and 5pb in 3'
    fastp --trim_front1 25 --trim_tail1 5 \
          --trim_front2 25 --trim_tail2 5 \
          --out1 ${sample}_nonprimer_R1_temp.fq.gz \
          --out2 ${sample}_nonprimer_R2_temp.fq.gz \
          --in1 ${sample}_fastp_filtered_R1_temp.fq.gz \
          --in2 ${sample}_fastp_filtered_R2_temp.fq.gz
        
    #Fourth fastp: to trim bad quality bases and to merge PE reads
    fastp --verbose \
          --cut_tail --cut_tail_mean_quality 20 \
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

    samtools merge -o tem_${sample}_map_unsort.bam tem_${sample}_NM_map.bam tem_${sample}_M_map.bam; samtools sort -O bam -o final_${Sample}.bam tem_${sample}_map_unsort.bam
    bam=final_${sample}.bam
    samtools index ${bam}
        
    #The number of reads with good quality is determined
    num_3=$(samtools view -c ${bam})
    echo "Number of useful reads ${sample}: ${num_3}" >> reads.txt
    echo "${sample} ${num_3}" >> ${main_dir}/tem_no_dups_utils.txt

    #Reads with post-mortem molecular (PMD) damage is extracted
    echo "
    Starting PMDtools...
    "
    samtools view -h ${b} | pmdtools --threshold 1 --header | samtools view -Sb - > pmd_final_${sample}.bam    
    pam=pmd_final_${sample}.bam
    samtools index ${pam}
        
    #The number of reads with PMD is determined
    num_4=$(samtools view -c ${pam})
    echo "Number of damaged reads ${sample}: ${num_4}" >> reads.txt
    echo "${sample} ${num_4}" >> ${main_dir}/tem_no_dups_utils_pmd.txt
    
    #Mean depth coverage is extracted
    echo "
    Starting qualimap...
    "
    qualimap bamqc -bam ${bam} -c -gd hg19 -outdir 3_${sample}_qualimap_bamqc_all_reads
    qualimap bamqc -bam ${pam} -c -gd hg19 -outdir 4_${sample}_qualimap_bamqc_damaged_reads
    cov_b=$(grep "mean coverageData" 3_${sample}_qualimap_bamqc_all_reads/genome_results.txt | sed 's/mean coverageData = //' | sed 's/X//' | sed 's/,//')
    cov_p=$(grep "mean coverageData" 4_${sample}_qualimap_bamqc_damaged_reads/genome_results.txt | sed 's/mean coverageData = //' | sed 's/X//' | sed 's/,//')
    echo "Mean depth coverage (all reads) for ${sample}: ${cov_b}" >> reads.txt
    echo "Mean depth coverage (all reads) for ${sample}: ${cov_p}" >> reads.txt  
    echo "${bam} ${cov_b}" >> ${main_dir}/tem_cov.txt
    echo "${pam} ${cov_p}" >> ${main_dir}/tem_cov.txt

    #Ranges for positions with a depth coverage over 10 are generated
    echo "
    Range of positions with a depth coverage higher than 10 reads is being determined...
    "
    samtools depth ${bam} | awk -F "\t" '$3>9 {print($2)}' |  awk '{if ($1<670) print $1+15900; else print $1-669}' | awk  '$1 < 303 || $1 > 315 {print $1}' | sort -n > temp.txt
    range=$(python ${main_dir}/range.py temp.txt | sed "s/[',\[\)\( ]//g" | sed "s/]//g" | sed "s/;$//")
    rm temp.txt
    echo "${bam} ${range}" >> ${main_dir}/tem_range.txt
    samtools depth ${pam} | awk -F "\t" '$3>9 {print($2)}' |  awk '{if ($1<670) print $1+15900; else print $1-669}' | awk  '$1 < 303 || $1 > 315 {print $1}' | sort -n > temp.txt
    range=$(python ${main_dir}/range.py temp.txt | sed "s/[',\[\)\( ]//g" | sed "s/]//g" | sed "s/;$//")
    rm temp.txt
    echo "${pam} ${range}" >> ${main_dir}/tem_pmd_range.txt
    
    bases_with_good_cov=$(samtools depth ${bam} | awk '$3>9 {print $2"\t"1}' | awk '{sum += $2} END {print sum}')
    echo -e "${sample} ${bases_with_good_cov}" >> ${main_dir}/tem_goodcovinsamples.txt

    #Variant call is done with freebayes and vcfallelicprimitives
    echo "
    Starting the variant call
    "
    freebayes -f ${NCR} -i -X -F 0.3 -C 3 --min-coverage 10 ${bam} | vcfallelicprimitives -kg > ${bam}.vcf 
    awk '$3 == "." && length($4) == 1 && length($5) == 1 {print $2"\t"$4"\t"$5"\t"$10}' ${bam}.vcf | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"100*$9/$5}' | awk '{if ($6 < 70) print $1"\t"$2$3; else print $1"\t"$3;}' |  sed 's/,.*//' | sed 's/[AGTC]\{4\}/N/g' | sed 's/[ACG]\{3\}/V/g' |sed 's/[ATC]\{3\}/H/g' | sed 's/[ATG]\{3\}/D/g' | sed 's/[GTC]\{3\}/B/g' | sed 's/[AC]\{2\}/M/g' | sed 's/[TG]\{2\}/K/g' | sed 's/[AT]\{2\}/W/g' | sed 's/[GC]\{2\}/S/g' | sed 's/[TC]\{2\}/Y/g' | sed 's/[AG]\{2\}/R/g' | awk '{if ($1<670) print $1+15900$2; else print $1-669$2}' |   sed -e 's/\t$//g' | awk -F"\t" '$1 < 303 || $1 > 315 {print $1$2}' > ${bam}.txt
    freebayes -f ${NCR} -i -X -F 0.3 -C 3 --min-coverage 10 ${pam} | vcfallelicprimitives -kg > ${pam}.vcf
    awk '$3 == "." && length($4) == 1 && length($5) == 1 {print $2"\t"$4"\t"$5"\t"$10}' ${pam}.vcf | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"100*$9/$5}' | awk '{if ($6 < 70) print $1"\t"$2$3; else print $1"\t"$3;}' | sed 's/,.*//' | sed 's/[AGTC]\{4\}/N/g' | sed 's/[ACG]\{3\}/V/g' |sed 's/[ATC]\{3\}/H/g' | sed 's/[ATG]\{3\}/D/g' | sed 's/[GTC]\{3\}/B/g' | sed 's/[AC]\{2\}/M/g' | sed 's/[TG]\{2\}/K/g' | sed 's/[AT]\{2\}/W/g' | sed 's/[GC]\{2\}/S/g' | sed 's/[TC]\{2\}/Y/g' | sed 's/[AG]\{2\}/R/g' | awk '{if ($1<670) print $1+15900$2; else print $1-669$2}' | sed -e 's/\t$//g' | awk -F"\t" '$1 < 303 || $1 > 315 {print $1$2}' > ${pam}.txt
 
    awk '$3 == "." && length($4) == 1 && length($5) == 1 {if ($2<670) print $2+15900"\t"$5"\t"$10; else print $2-669"\t"$5"\t"$10}' ${bam}.vcf | awk  '$1 < 303 || $1 > 315 {print $0}' | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$4"\t"$8"\t"100*$8/$4"\t"}' | awk -v vd=${sample} '$5 < 70 {print vd"\t"$0}' >> ${main_dir}/base_mix.txt
 
    cp ${bam}.txt ${main_dir}/${vis}
    cp ${pam}.txt ${main_dir}/${vis}
    
    rm *tem*
    echo "
    Changing sample...
    " 
    cd ${main_dir}

done < temp_samples.txt

#File in hsd format (the one requiered by HaploGrep2) is generated in "Visualization directory"
cd ${main_dir}/${vis}
echo "
Creating the hsd file to use in HaploGrep2
"
echo "Sample Range Haplogroup Haplotype" > tem_haplogrep_format.txt
#Variable with the same value for all samples is generated
haplogroup="?"
#The rest of the variables are generated line by line
for t in *final*bam.txt; do
    #Name is extracted
    sample="${t//*final_/}"
    name="${t//.txt/}"
    if [[ ${t} == pmd* ]]
    then
        sample="${sample//.*/_pmd}"
        ranges=$(awk -v sample=${name} '$1==sample {print $2}' ${main_dir}/tem_pmd_range.txt)
    else
        sample="${sample//.*/}"
        ranges=$(awk -v sample=${name} '$1==sample {print $2}' ${main_dir}/tem_range.txt)
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
echo -e 'ID\tTotal_reads\tReads_without_dups\tDuplication_rate\tUtil_reads\tPercentage_of_util_reads\tMean_depth_coverage\tDamaged_util_reads\tPercentage_of_Damaged_reads\tHeteroplamies\tHaplogrouop\tQuality\tPercentage_of_NCR_recovered\tRange_cov>10\tHaplotype' > final_table.txt

while read sample; do
    echo "${sample}" | awk -v v0=${sample} -v v1=$(awk -v v1_1=${sample} '$1==v1_1 {print $2}' tem_all.txt) -v v2=$(awk -v v2_1=${sample} '$1==v2_1 {print $2}' tem_no_dups.txt) -v v3=$(awk -v v3_1=${sample} '$1==v3_1 {print $2}' tem_no_dups_utils.txt) -v v4=$(grep -v pmd tem_cov.txt | sed 's/.bam//g' | sed 's/final_//g' | awk -v v4_1=${sample} '$1==v4_1 {print $2}') -v v5=$(awk -v v5_1=${sample} '$1==v5_1 {print $2}' tem_no_dups_utils_pmd.txt) -v v6=$(awk -F"\t" -v v6_1=${sample} '$1 == v6_1 {print $0}' base_mix.txt | wc -l) -v v7=$(awk -v v7_1=\"${sample}\" '$1==v7_1 {print $2"\t"$4}' tem_haplos.txt | sed 's/"//g' | sed 's/\t/--/g') -v v8=$(cat tem_range.txt | sed 's/.bam//g' | sed 's/final_//g' | awk -v v8_1=${sample} '$1==v8_1 {print $2}') -v v9=$(awk -v v9_1=${sample} '$1==v9_1 {print $0}' tem_haplogrep_format.txt | cut -f 4- | sed 's/\t/__/g') '$1==v0 {print v0"\t"v1"\t"v2"\t"(1-v2/v1)*100"\t"v3"\t"v3/v1*100"\t"v4"\t"v5"\t"v5/v3*100"\t"v6"\t"v7"\t"v8"\t"v9}' | sed 's/--/\t/g' | sed 's/__/ /g' |  awk -F"\t" -v vn1=$(awk -v vn1_1=${sample} '$1==vn1_1 {print $2}' tem_goodcovinsamples.txt)  '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"100*vn1/1198"\t"$13"\t"$14}' | sed 's/16519T$//g' | sort -h >> final_table.txt
done < temp_samples.txt

rm *tem*.txt
            
echo "          
                THE END
     "

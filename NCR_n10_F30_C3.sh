####SCRIPT PARA ANALIZAR  EL D-LOOP DE MUESTRAS ANTIGUAS (GUERRA CIVIL Y ARQUEOLÓGICAS) GENERADAS CON POWERSEQ (MEDIANTE AMPLICONES). EL SCRIPT TE GENERA EL ALINEAMIENTO, SE QUEDA SOLAMENTE CON LOS READS QUE MAPEAN, ELIMINA LOS DUPLICADOS Y REALINEA SEGÚN LOS INDELS. ADEMÁS, GENERA UN BAM CON LOS READS QUE PRESENTAN DAÑO MOLECULAR. A PARTIR DE ESTOS DOS ARCHIVOS, SI TIENEN UN BUEN COVERAGE (9 PARA LA MUESTRA ORIGINAL Y 4 PARA LA MUESTRA SOLAMENTE CON DAÑO) SACA LOS POLIMORFISMOS. EN EL CASO DE QUE LA MUESTRA CON DAÑO MOLECULAR TENGA UN BAJO COVERAGE SOLO SE QUEDA CON LA ORIGINAL, PERO SI LAS DOS TIENE BUEN COVERAGE EN LA MUESTRA CON DAÑO MOLECULAR SOLO SE QUEDARÁ CON LOS POLIMORFISMOS COMUNES (ELIMINANDO DE ESTA FORMA, LOS POSIBLE POLIMORFISMOS QUE APARECEN AL REDUCIR LOS UMBRALES GENERADOS POR EL DAÑO MOLECULAR). AL FINAL GENERARÁ UN ARCHIVO EN EL FORMATO NECESARIO PARA SUBIR DIRECTAMENTE AL HAPLOIGREP####

source activate ngs_analysis

#Guardaremos la ruta de la carpeta en la que se encuentran varios ejecutables
PICARD=/home/uab/anaconda3/envs/ngs_analysis/share/picard-2.26.2-0
PMDTOOLS=/home/uab/anaconda3/envs/ngs_analysis/share/pmdtools-0.60-3/pmdtools.0.60.py
echo "Start"

#Creamos una carpeta donde se guardarán los resumenes para usar multiqc (si no existe)
vis=visualization
if [[ -d ${vis} ]]; then
    echo "Ya exite el directorio ${vis}, por lo que no se volverá a crear"
else
    echo "Creamos un directorio para guardar los resúmenes"
    mkdir ${vis}
fi
#Se crea una subcarpeta para guardar los bams y bais finales
if [[ -d ${vis}/bambai ]]; then
    echo "Ya exite el directorio para guardar los bams y bais, por lo que no se volverá a crear"
else
    echo "Creamos un directorio para guardar los bams y los bais"
    mkdir ${vis}/bambai
fi
file=$(pwd)
bambai=${file}/${vis}/bambai

#Creamos una carpeta donde guardaremos la referencia, si no existe. Además, guardaremos la ruta en la que se encuentra la referencia
ref=referencia
re2=rCRS_dloop_aprox.fasta

#Si existe la carpeta envía un mensaje. Si no existe, creala y copia las referencias en ella
if [[ -d ${ref} ]]; then
    echo "Ya exite el directorio ${ref}, por lo que no se volverá a crear"
    dloop_aprox=$(pwd)/${ref}/${re2}
else
    echo "Creamos un directorio para guardar la secuencia de referencia"
    mkdir ${ref}
    cp ${re2} ${ref}
fi
#Entramos en la carpeta (para comprobar que estén las secuencias ya hechas)
cd ${ref}
#Si ya existen los índices, se enviará un mensaje, en caso contrartio, se crearán.
if [[ -e ${re2}.amb ]] && [[ -e ${re2}.fai ]] ; then
    echo "
    Ya esxiste el index de la secuencia de referencia
    "
    dloop_aprox=$(pwd)/${re2}
else
    echo "
    Se crea el index de la secuencia de referencia
    "
    time bwa index ${re2}
    time samtools faidx ${re2}
    time java -jar ${PICARD}/picard.jar CreateSequenceDictionary R=${re2} O=rCRS_dloop_aprox.dict
    
    dloop=$(pwd)/${re}
    dloop_aprox=$(pwd)/${re2}
    #conda deactivate
fi
cd ..

#Se actuará sobre todas las carpetas (salvo la carpeta visualización y la que contiene la referencia)
for d in *; do if  [[ -d ${d} ]] && [[ ${d} != ${vis} ]] && [[ ${d} != ${ref} ]]; then  echo $d; fi done > temp_samples.txt

echo -e "Sample\tPosition\tAlternative_base\tN_total\tN_alternative\tPercentage" > base_mix.txt

while read d; do
      
    echo "
    Empezando con la muestra ${d}
    "
    cd "$d"
        
    dfile=$(pwd)
   
    #Se hace un fastq
    echo "
    Hacemos un fastq
    "
    fastq_file=fastqc_file_quality
    if [[ -d ${fastq_file} ]]; then
        echo "Ya exite el directorio ${fastq_file}, por lo que no se volverá a crear"
    else
        echo "Creamos un directorio para guardar los resúmenes"
        mkdir ${fastq_file}
        time fastqc -o ${fastq_file} -f fastq *f*q.gz
    fi
   
   #Hacemos un primer fastp para eliminar los duplicados
   fastp --verbose \
        --dedup \
        --disable_length_filtering \
        --disable_quality_filtering \
        --json ${d}_dedup_temp.json \
        --thread 8 \
        --out1 ${d}_fastp_dedup_R1_temp.fq.gz \
        --out2 ${d}_fastp_dedup_R2_temp.fq.gz \
        --in1 *R1*.fastq.gz \
        --in2 *R2*.fastq.gz
        
    #Se cuenta el número de reads con duplicados
    num_1=$(grep before_filtering -m 1 -A1 ${d}_dedup_temp.json | grep total | sed "s/.*://g" | sed "s/,//g")
    echo "Number of total reads ${d}: ${num_1}" >> reads.txt
    echo "Number of total reads ${d}: ${num_1}" >> ${file}/reads.txt
    echo "${d} ${num_1}" >> ${file}/tem_all.txt
    
    #Se cuenta el número de reads sin duplicados
    num_2=$(grep after_filtering -m 1 -A1 ${d}_dedup_temp.json | grep total | sed "s/.*://g" | sed "s/,//g")
    echo "Number of total reads without duplicates ${d}: ${num_2}" >> reads.txt
    echo "Number of total reads without duplicates ${d}: ${num_2}" >> ${file}/reads.txt
    echo "${d} ${num_2}" >> ${file}/tem_no_dups.txt
    
    #Second fastp: to trim adapters
    fastp --verbose \
        --detect_adapter_for_pe \
        --trim_poly_x \
        --length_required 30 \
        --thread 8 \
        --out1 ${d}_fastp_filtered_R1_temp.fq.gz \
        --out2 ${d}_fastp_filtered_R2_temp.fq.gz \
        --in1 ${d}_fastp_dedup_R1_temp.fq.gz \
        --in2 ${d}_fastp_dedup_R2_temp.fq.gz
        
    #Third fastp: to trim primers (25pb from 5') and 5pb in 3'
    time fastp --verbose \
        --trim_front1 25 --trim_tail1 5 \
        --trim_front2 25 --trim_tail2 5 \
        --out1 ${d}_nonprimer_R1_temp.fq.gz \
        --out2 ${d}_nonprimer_R2_temp.fq.gz \
        --in1 ${d}_fastp_filtered_R1_temp.fq.gz \
        --in2 ${d}_fastp_filtered_R2_temp.fq.gz
        
    #Fourth fastp: To trim bad quality bases
    time fastp --verbose \
        --cut_tail --cut_tail_mean_quality 20 \
        --trim_poly_x \
        --merge --overlap_len_require 10 \
        --correction \
        --html ${d}_temp.html \
        --json ${d}_temp.json \
        --merged_out ${d}_M.fq.gz \
        --out1 ${d}_R1.fq.gz \
        --out2 ${d}_R2.fq.gz \
        --in1 ${d}_nonprimer_R1_temp.fq.gz \
        --in2 ${d}_nonprimer_R2_temp.fq.gz
    
    rm *_temp*

    #Se hace un fastq
    echo "
    Hacemos un fastq
    "
    fastq_file=fastqc_post_adapter_trimming
    if [[ -d ${fastq_file} ]]; then
        echo "Ya exite el directorio ${fastq_file}, por lo que no se volverá a crear"
    else
        echo "Creamos un directorio para guardar los resúmenes"
        mkdir ${fastq_file}
        time fastqc -o ${fastq_file} -f fastq *.fq.gz
    fi
        
    #Se hace el alineamiento    
    echo "
    Empezamos el alineamiento
    "
    echo "
    Alineamiento (1/2)
    "
    bwa aln -t 9 -l 1024 -n 0.01 -o 2 ${dloop_aprox} ${d}_R1.fq.gz > aln_R1.sai
    bwa aln -t 9 -l 1024 -n 0.01 -o 2 ${dloop_aprox} ${d}_R2.fq.gz > aln_R2.sai
    bwa aln -t 9 -l 1024 -n 0.01 -o 2 ${dloop_aprox} ${d}_M.fq.gz > aln_M.sai
    echo "
    Alineamiento (2/2)
    "
    bwa sampe ${dloop_aprox} aln_R1.sai aln_R2.sai ${d}_R1.fq.gz ${d}_R2.fq.gz | samtools view -Sbh -F4 -q30 | samtools sort -O bam -o tem_${d}_NM_map.bam
    bwa samse ${dloop_aprox} aln_M.sai ${d}_M.fq.gz | samtools view -Sbh -F4 -q30 | samtools sort -O bam -o tem_${d}_M_map.bam

    samtools merge -o tem_${d}_map_unsort.bam tem_${d}_NM_map.bam tem_${d}_M_map.bam; samtools sort -O bam -o final_${d}.bam tem_${d}_map_unsort.bam
    b=final_${d}.bam
    samtools index ${b}
        
    #Se calcula el número de reads
    num_3=$(samtools view -c ${b})
    echo "Number of util reads ${b}: ${num_3}" >> reads.txt
    echo "Number of util reads ${d}: ${num_3}" >> ${file}/reads.txt
    echo "${d} ${num_3}" >> ${file}/tem_no_dups_utils.txt

    #Obtenemos los reads que presentan daño molecular
    echo "
    Nos quedamos con los reads que tienen daño molecular
    "
    time samtools view -h ${b} | python ${PMDTOOLS} --threshold 1 --header | samtools view -Sb - > pmd_final_${d}.bam
            
    p=pmd_final_${d}.bam
    samtools index ${p}
        
    #Contamos los reads que tienen daño molecular
    num_4=$(samtools view -c ${p})
    echo "Number of damaged reads ${p}: ${num_4}" >> reads.txt
    echo "Number of damaged reads pmd_${d}: ${num_4}" >> ${file}/reads.txt
    echo "${d} ${num_4}" >> ${file}/tem_no_dups_utils_pmd.txt
       
    for f in *final_${d}.bam; do
        #Miramos la calidad de los bams
        echo "
        Fem qualimap
        "
        time qualimap bamqc -bam ${f} -c -gd hg19 -outdir ${f}_qualimap_bamqc  
        echo "Calidad de la muestra ${f}: $(grep "mean mapping quality" ${f}_qualimap_bamqc/genome_results.txt | sed 's/mean mapping quality = //' )" >> reads.txt
        echo "Calidad de la muestra ${f}: $(grep "mean mapping quality" ${f}_qualimap_bamqc/genome_results.txt | sed 's/mean mapping quality = //')" >> ${file}/reads.txt
        cov=$(grep "mean coverageData" ${f}_qualimap_bamqc/genome_results.txt | sed 's/mean coverageData = //' | sed 's/X//' | sed 's/,//')
        echo "Coverage medio de la muestra ${f}: ${cov}" >> reads.txt             
        echo "Coverage medio de la muestra ${f}: ${cov}" >> ${file}/reads.txt
        echo "${f} ${cov}" >> ${file}/tem_cov.txt
    done

    echo "
    Se calcula el rango de posiciones con más de 10 reads
    "
    for X in final*bam; do
        echo $X
        #Se crean los rangos en los que nos encontramos más de 10 reads
        samtools depth ${X} | awk -F "\t" '$3>9 {print($2)}' |  awk '{if ($1<670) print $1+15900; else print $1-669}' | awk  '$1 < 303 || $1 > 315 {print $1}' | sort -n > temp.txt
        range=$(python ${file}/range.py temp.txt | sed "s/[',\[\)\( ]//g" | sed "s/]//g" | sed "s/;$//")
        rm temp.txt
        echo "${X} ${range}" >> ${file}/tem_range.txt
    done
           
    #Se mirará si hay suficiente coverage de muestra, teniendo solo en cuenta los reads con daño molecular, que será cuando el "depth coverage" sea mayor que 0
    if [[ $(grep "mean coverageData" ${p}_qualimap_bamqc/genome_results.txt | sed 's/\.[0-9]*//' | sed 's/,//' | sed 's/[a-z]//ig' | sed "s/=//") -lt 0 ]]; then
    
        echo "
        Pasamos del bam que solo contiene daño molecular porque tiene mal coverage: $(grep "mean coverageData" ${p}_qualimap_bamqc/genome_results.txt | sed 's/mean coverageData = //' | sed 's/X//' | sed 's/,//')
        "
        echo "Bad coverage for sample ${b}_pmd: $(grep "mean coverageData" ${p}_qualimap_bamqc/genome_results.txt | sed 's/mean coverageData = //' | sed 's/X//' | sed 's/,//')" >> ${file}/error.txt
        cp ${p} ${bambai}
        cp ${p}.bai ${bambai} 
    else
        echo "
        Buscamos los polimporfismos porque tiene buen coverage: $(grep "mean coverageData" ${p}_qualimap_bamqc/genome_results.txt | sed 's/mean coverageData = //' | sed 's/X//' | sed 's/,//')
        "
        time freebayes -f ${dloop_aprox} -i -X -F 0.3 -C 3 --min-coverage 10 ${p} | vcfallelicprimitives -kg > ${p}.vcf
        awk '$3 == "." && length($5) == 1 && length($4) == 1 {print $2"\t"$4"\t"$5"\t"$10}' ${p}.vcf | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"100*$9/$5}' | awk '{if ($6 < 70) print $1"\t"$2$3; else print $1"\t"$3;}' | sed 's/[AGTC]\{4\}/N/g' | sed 's/[ACG]\{3\}/V/g' |sed 's/[ATC]\{3\}/H/g' | sed 's/[ATG]\{3\}/D/g' | sed 's/[GTC]\{3\}/B/g' | sed 's/[AC]\{2\}/M/g' | sed 's/[TG]\{2\}/K/g' | sed 's/[AT]\{2\}/W/g' | sed 's/[GC]\{2\}/S/g' | sed 's/[TC]\{2\}/Y/g' | sed 's/[AG]\{2\}/R/g' | awk '{if ($1<670) print $1+15900$2; else print $1-669$2}' | sed -e 's/\t$//g' | awk -F"\t" '$1 < 303 || $1 > 315 {print $1$2}' > ${p}.txt
        
        
        for X_PMD in pmd_final*bam; do
            echo $X_PMD
            #Se crean los rangos en los que nos encontramos más de 10 reads
            samtools depth ${X_PMD} | awk -F "\t" '$3>9 {print($2)}' |  awk '{if ($1<670) print $1+15900; else print $1-669}' | awk  '$1 < 303 || $1 > 315 {print $1}' | sort -n > temp_pmd.txt
            range=$(python ${file}/range.py temp_pmd.txt | sed "s/[',\[\)\( ]//g" | sed "s/]//g" | sed "s/;$//")
            rm temp_pmd.txt
            echo "${X_PMD} ${range}" >> ${file}/tem_pmd_range.txt
        done
    fi
        
    echo "
    Buscamos los polimporfismos porque tiene buen coverage: $(grep "mean coverageData" ${b}_qualimap_bamqc/genome_results.txt | sed 's/mean coverageData = //' | sed 's/X//' | sed 's/,//')
    "
    time freebayes -f ${dloop_aprox} -i -X -F 0.3 -C 3 --min-coverage 10 ${b} | vcfallelicprimitives -kg > ${b}.vcf 
    awk '$3 == "." && length($5) == 1 && length($4) == 1 {print $2"\t"$4"\t"$5"\t"$10}' ${b}.vcf | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"100*$9/$5}' | awk '{if ($6 < 70) print $1"\t"$2$3; else print $1"\t"$3;}' |  sed 's/,.*//' | sed 's/[AGTC]\{4\}/N/g' | sed 's/[ACG]\{3\}/V/g' |sed 's/[ATC]\{3\}/H/g' | sed 's/[ATG]\{3\}/D/g' | sed 's/[GTC]\{3\}/B/g' | sed 's/[AC]\{2\}/M/g' | sed 's/[TG]\{2\}/K/g' | sed 's/[AT]\{2\}/W/g' | sed 's/[GC]\{2\}/S/g' | sed 's/[TC]\{2\}/Y/g' | sed 's/[AG]\{2\}/R/g' | awk '{if ($1<670) print $1+15900$2; else print $1-669$2}' |   sed -e 's/\t$//g' | awk -F"\t" '$1 < 303 || $1 > 315 {print $1$2}' > ${b}.txt
    

    awk '$3 == "." && length($4) == 1 && length($5) == 1 {if ($2<670) print $2+15900"\t"$5"\t"$10; else print $2-669"\t"$5"\t"$10}' ${b}.vcf | awk  '$1 < 303 || $1 > 315 {print $0}' | sed 's/:/\t/g' | awk '{print $1"\t"$2"\t"$4"\t"$8"\t"100*$8/$4"\t"}' | awk -v vd=$d '$5 < 70 {print vd"\t"$0}' >> ${file}/base_mix.txt

    
    for b2 in *final_${d}.bam; do
        cp ${b2}.txt ${file}/${vis}
    done
    
    rm *tem*
    echo "
    Cambio de muestra
    " 
    cd ${file}

done < temp_samples.txt

#Se generas el archivo para Haplogrep en la carpeta de resultados llamada visualización
cd ${file}/${vis}
echo "
Generamos el archivo para el haplogrep
"
#Creamos el archivo para el haplogrep con los envabezados de las columnas
echo "Sample Range Haplogroup Haplotype" > tem_haplogrep_format.txt
#Generamos las variables que son comunes para todas las muestras
haplogroup="?"
#Iremos añadiendo línea a línea para cada muestra
for t in *final*bam.txt; do
    #Nos quedamos con el nombre
    sample="${t//*final_/}"
    name="${t//.txt/}"
    if [[ ${t} == pmd* ]]
    then
        sample="${sample//.*/_pmd}"
        ranges=$(awk -v sample=${name} '$1==sample {print $2}' ${file}/tem_pmd_range.txt)
    else
        sample="${sample//.*/}"
        ranges=$(awk -v sample=${name} '$1==sample {print $2}' ${file}/tem_range.txt)
    fi
    #Cambiamos los falsos tabuladores de plink por tabuladores reales
    sed -i 's/    /\t/g' ${t}
    #Imprimimos en un archivo temporal los polimorfismos (en vertical)
    awk -F"\t" '{print $1 $2}' ${t} > temp_file.txt
    #Trasnformamos la columna del archivo temporal en una fila y la guardamos en una variable
    haplotype=$(for i in `< temp_file.txt`; do echo -en ${i}" "; done)
    #Añadimos la fila en el formato hsd del haplogrep (separados por espacios)
    echo "${sample} ${ranges} ${haplogroup} ${haplotype}" >> tem_haplogrep_format.txt
done
#Trasnformamos los espacios en tabuladores, para que sea el formato correcto de haplogrep
sed 's/ /\t/g' tem_haplogrep_format.txt | sed 's/\t$//g' | awk -F "\t" '{ if ($4 == "") print $0"\t16519T"; else print $0}' > haplogrep_format.txt
cp haplogrep_format.txt ${file}
#Se borran los archivos temporales
#rm *tem*

cd ${file}

grep -v pmd haplogrep_format.txt > tem_haplogrep_format.txt
bash haplogrep classify --in tem_haplogrep_format.txt --format hsd --out tem_haplos.txt

echo "
Generamos una tabla final
"
echo -e 'ID\tTotal_reads\tReads_without_dups\tDuplication_rate\tUtil_reads\tPercentage_of_util_reads\tMean_depth_coverage\tDamaged_util_reads\tPercentage_of_Damaged_reads\tHeteroplamies\tHaplogrouop\tQuality\tPercentage_of_NCR_recovered\tRange_cov>10\tHaplotype' > final_table.txt

while read d; do
bases_buen_cov=$(samtools depth ${file}/${d}/final_${d}.bam | awk '$3>10 {print $2"\t"1}' | awk '{sum += $2} END {print sum}')
echo -e "$d\t$bases_buen_cov" >> tem_goodcovinsamples.txt
done < temp_samples.txt

for d in *; do
    if [[ -d ${d} ]] && [[ ${d} != "visualization" ]] && [[ ${d} != "referencia" ]] && [[ ${d} != "final_analysis" ]]; then
        echo "${d}" | awk -v v0=${d} -v v1=$(awk -v v1_1=${d} '$1==v1_1 {print $2}' tem_all.txt) -v v2=$(awk -v v2_1=${d} '$1==v2_1 {print $2}' tem_no_dups.txt) -v v3=$(awk -v v3_1=${d} '$1==v3_1 {print $2}' tem_no_dups_utils.txt) -v v4=$(grep -v pmd tem_cov.txt | sed 's/.bam//g' | sed 's/final_//g' | awk -v v4_1=${d} '$1==v4_1 {print $2}') -v v5=$(awk -v v5_1=${d} '$1==v5_1 {print $2}' tem_no_dups_utils_pmd.txt) -v v6=$(awk -F"\t" -v v6_1=${d} '$1 == v6_1 {print $0}' base_mix.txt | wc -l) -v v7=$(awk -v v7_1=\"${d}\" '$1==v7_1 {print $2"\t"$4}' tem_haplos.txt | sed 's/"//g' | sed 's/\t/--/g') -v v8=$(grep -v pmd tem_range.txt | sed 's/.bam//g' | sed 's/final_//g' | awk -v v8_1=${d} '$1==v8_1 {print $2}') -v v9=$(awk -v v9_1=${d} '$1==v9_1 {print $0}' tem_haplogrep_format.txt | cut -f 4- | sed 's/\t/__/g') '$1==v0 {print v0"\t"v1"\t"v2"\t"(1-v2/v1)*100"\t"v3"\t"v3/v1*100"\t"v4"\t"v5"\t"v5/v3*100"\t"v6"\t"v7"\t"v8"\t"v9}' | sed 's/--/\t/g' | sed 's/__/ /g' |  awk -F"\t" -v vn1=$(awk -v vn1_1=${d} '$1==vn1_1 {print $2}' tem_goodcovinsamples.txt)  '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"100*vn1/1198"\t"$13"\t"$14}' | sed 's/16519T$//g' | sort -h >> final_table.txt
    fi
done            

rm *tem*.txt
            
echo "          
                THE END
     "

time bash pmds.sh

## time bash NCR_n10_F30_C3.sh 2>&1 | tee NCR_n10_F30_C3.log

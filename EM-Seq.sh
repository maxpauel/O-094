#!/bin/bash

for file in ./Fastq/*R1_001.fastq.gz; do
    new_name="${file%%_S_*}_R1.fastq.gz"
    mv "$file" "$new_name"
done



# Пути (замените на свои)
INPUT_DIR="./Fastq"
OUTPUT_DIR="./Fastq_trim"
ADAPTERS="/media/maxpauel/6372EF8A3B0879C1/Methylom/adapters/adapters.fa"
THREADS=10

# Проверяем, существует ли выходная директория
mkdir -p "$OUTPUT_DIR"


for R1_FILE in "$INPUT_DIR"/*_R1.fastq.gz; do
    # Получаем соответствующий R2 файл
    R2_FILE="${R1_FILE/_R1.fastq.gz/_R2.fastq.gz}"


    # Имя образца (без _R1.fastq.gz)
    SAMPLE_NAME=$(basename "$R1_FILE" _R1.fastq.gz)

    # Запускаем Trimmomatic
    java -jar /usr/share/java/trimmomatic.jar PE \
        -threads "$THREADS" \
        "$R1_FILE" "$R2_FILE" \
        "$OUTPUT_DIR/${SAMPLE_NAME}_FP.fastq.gz" "$OUTPUT_DIR/${SAMPLE_NAME}_FU.fastq.gz" \
        "$OUTPUT_DIR/${SAMPLE_NAME}_RP.fastq.gz" "$OUTPUT_DIR/${SAMPLE_NAME}_RU.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 >> ./Fastq_trim/trimmomatic.log 2>&1
    echo "-----------------------------------" >> ./Fastq_trim/trimmomatic.log
done

# reformat trimmomatic log
grep -E '^Processing|^Input Reads' ./Fastq_trim/trimmomatic.log > ./Fastq_trim/trimmed_summary.tsv



# Пути (замените на свои)
REFERENCE="/media/maxpauel/6372EF8A3B0879C1/Methylom/bwa_meth_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
FASTQ_DIR="./Fastq_trim/"
BAM_DIR="./bam/"
THREADS=10

for R1_FASTQ in "$FASTQ_DIR"/*_FP.fastq.gz; do
    R2_FASTQ="${R1_FASTQ/_FP.fastq.gz/_RP.fastq.gz}"
    SAMPLE_NAME=$(basename "$R1_FASTQ" _FP.fastq.gz)
    bwameth.py -t "$THREADS" --reference "$REFERENCE" \
        "$R1_FASTQ" "$R2_FASTQ" > "$BAM_DIR/$SAMPLE_NAME.sam"
    # 2. samtools view (конвертация SAM → BAM)
    samtools view -@ "$THREADS" -bS "$BAM_DIR/$SAMPLE_NAME.sam" > "$BAM_DIR/$SAMPLE_NAME.bam"
    # 3. samtools sort (сортировка BAM)
    samtools sort -@ "$THREADS" "$BAM_DIR/$SAMPLE_NAME.bam" -o "$BAM_DIR/$SAMPLE_NAME.sorted.bam"
    # 4. Удаляем временные файлы
    rm "$BAM_DIR/$SAMPLE_NAME.sam" "$BAM_DIR/$SAMPLE_NAME.bam"
done


# Проходим по всем файлам *_FU.fastq.gz (непарные forward reads)
for FU_FASTQ in "$FASTQ_DIR"/*_FU.fastq.gz; do
    RU_FASTQ="${FU_FASTQ/_FU.fastq.gz/_RU.fastq.gz}"
    SAMPLE_NAME=$(basename "$FU_FASTQ" _FU.fastq.gz)
    bwameth.py -t "$THREADS" --reference "$REFERENCE" \
        "$FU_FASTQ,$RU_FASTQ" > "$BAM_DIR/$SAMPLE_NAME.u.sam"
    samtools view -@ "$THREADS" -bS "$BAM_DIR/$SAMPLE_NAME.u.sam" > "$BAM_DIR/$SAMPLE_NAME.u.bam"
    samtools sort -@ "$THREADS" "$BAM_DIR/$SAMPLE_NAME.u.bam" -o "$BAM_DIR/$SAMPLE_NAME.u.sorted.bam"
   rm "$BAM_DIR/$SAMPLE_NAME.u.sam" "$BAM_DIR/$SAMPLE_NAME.u.bam"
done

for SORTED_BAM in "$BAM_DIR"/*.sorted.bam; do
    # Пропускаем уже объединенные файлы
    if [[ "$SORTED_BAM" == *".u.sorted.bam" ]]; then
        continue
    fi

    SAMPLE_NAME=$(basename "$SORTED_BAM" ".sorted.bam")
    UNPAIRED_BAM="$BAM_DIR/$SAMPLE_NAME.u.sorted.bam"

    if [[ -f "$UNPAIRED_BAM" ]]; then
        samtools merge -@ "$THREADS" \
            -o "$BAM_DIR/$SAMPLE_NAME.merged.bam" \
            "$SORTED_BAM" "$UNPAIRED_BAM"
        
        # Опционально: удаляем исходные отсортированные файлы после объединения
        # rm "$SORTED_BAM" "$UNPAIRED_BAM"
    fi
done

METHYL_DACKEL="./MethylDackel"

for MERGED_BAM in "$BAM_DIR"/*.merged.bam; do
    # Получаем имя образца (без .merged.bam)
    SAMPLE_NAME=$(basename "$MERGED_BAM" .merged.bam)
    
    # 1. MethylDackel mbias

    "$METHYL_DACKEL" mbias "$REFERENCE" "$MERGED_BAM" "$BAM_DIR/$SAMPLE_NAME"
    
    # 2. MethylDackel extract

    "$METHYL_DACKEL" extract "$REFERENCE" "$MERGED_BAM" -o "$BAM_DIR/$SAMPLE_NAME"
    
done




# Configuration
THREADS=10
CONTROL_REF="/media/maxpauel/6372EF8A3B0879C1/Methylom3/control_library/controls.fa"
MAIN_BAM_DIR="./bam"
CONTROL_OUTPUT_DIR="/media/maxpauel/6372EF8A3B0879C1/Methylom3/control_library"
METHYL_DACKEL="./MethylDackel2"
CHROMOSOMES_FILE="/media/maxpauel/6372EF8A3B0879C1/Methylom3/control_library/bam/chromosomes.txt"  # Should contain Lambda_NEB and pUC19



# Process each sorted BAM file
for SORTED_BAM in "${MAIN_BAM_DIR}"/*.sorted.bam; do
    # Skip if no files found
    [ -e "$SORTED_BAM" ] || continue
    
    # Get sample name
    SAMPLE_NAME=$(basename "${SORTED_BAM}" .sorted.bam)
    
    
    # 1. Extract unmapped reads

    UNMAPPED_BAM="${CONTROL_OUTPUT_DIR}/${SAMPLE_NAME}.unmapped.bam"
    samtools view -@ ${THREADS} -f 4 -o "${UNMAPPED_BAM}" "${SORTED_BAM}"
    
    # 2. Convert unmapped BAM to FASTQ
    FQ1="${CONTROL_OUTPUT_DIR}/trash/${SAMPLE_NAME}.unmapped_f.fq"
    FQ2="${CONTROL_OUTPUT_DIR}/trash/${SAMPLE_NAME}.unmapped_r.fq"
    bedtools bamtofastq -i "${UNMAPPED_BAM}" -fq "${FQ1}" -fq2 "${FQ2}"
    
    # 3. Align to control sequences

    SAM_OUTPUT="${CONTROL_OUTPUT_DIR}/bam/${SAMPLE_NAME}_sh.sam"
    bwameth.py -t 10 --reference "${CONTROL_REF}" "${FQ1}" "${FQ2}" > "${SAM_OUTPUT}"
    
    # 4. Convert SAM to sorted BAM

    BAM_OUTPUT="${CONTROL_OUTPUT_DIR}/bam/${SAMPLE_NAME}_sh.bam"
    SORTED_OUTPUT="${CONTROL_OUTPUT_DIR}/bam/${SAMPLE_NAME}_sh.sorted.bam"
    
    samtools view -@ ${THREADS} -bS "${SAM_OUTPUT}" > "${BAM_OUTPUT}"
    samtools sort -@ ${THREADS} "${BAM_OUTPUT}" -o "${SORTED_OUTPUT}"
    samtools index "${SORTED_OUTPUT}"
    

    rm "${SAM_OUTPUT}" "${BAM_OUTPUT}"
    

    while read -r chr; do
        CHR_BAM="${CONTROL_OUTPUT_DIR}/bam/${chr}.${SAMPLE_NAME}.bam"
        samtools view -b -o "${CHR_BAM}" "${SORTED_OUTPUT}" "${chr}"
        
        # Run MethylDackel for each control chromosome
        "${METHYL_DACKEL}" extract "${CONTROL_REF}" "${CHR_BAM}" \
            -o "${CONTROL_OUTPUT_DIR}/bam/${chr}.${SAMPLE_NAME}"
    done < "${CHROMOSOMES_FILE}"
    

done



#!/bin/bash

INPUT_DIR="./Fastq"
OUTPUT_DIR="./Fastq_trim"
ADAPTERS="./adapters.fa"
THREADS=8

# Переходим в директорию с входными файлами
cd "$INPUT_DIR"

# Объединение файлов по префиксам (если нужно)
for prefix in $(ls *.fastq.gz | cut -d '_' -f 1 | sort -u); do
  cat "${prefix}"_*.fastq.gz > "${prefix}.fastq.gz"
done

# Запись количества ридов до обработки
for file in *.fastq.gz; do
    echo -n "$file: " >> read_counts.txt
    zcat "$file" | wc -l | awk '{print $1/4}' >> read_counts.txt
done

# Создаем выходную директорию, если её нет
mkdir -p "$OUTPUT_DIR"

# Обработка файлов с помощью Trimmomatic
for input_file in *.fastq.gz; do
    # Определяем имя выходного файла
    output_file="${OUTPUT_DIR}/${input_file%.fastq.gz}.trimmed.fastq.gz"
    
    # Запускаем Trimmomatic
    echo "Processing $input_file..." >> trimmomatic.log
    java -jar /usr/share/java/trimmomatic.jar SE \
        -threads "$THREADS" \
        "$input_file" \
        "$output_file" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:35 >> trimmomatic.log 2>&1
    echo "-----------------------------------" >> trimmomatic.log
done

# Создаем сводную таблицу результатов
grep -E '^Processing|^Input Reads' trimmomatic.log > "${OUTPUT_DIR}/trimmed_summary.tsv"

# Определение переменных
HT_INDEX="./Index/Index"
INPUT_DIR="./Fastq_trim"
BAM_OUTPUT_DIR="./bam"
COUNTS_OUTPUT_DIR="./counts"
SPLICE_SITES="./Homo_sapiens.GRCh38.101.splicesites.txt"
THREADS=8

# Создание выходных директорий
mkdir -p "$BAM_OUTPUT_DIR"
mkdir -p "$COUNTS_OUTPUT_DIR"

# Картирование с помощью HISAT2
for input_file in "$INPUT_DIR"/*.fastq.gz.trim.fastq.gz; do
    # Определение имен выходных файлов
    base_name=$(basename "$input_file" .fastq.gz.trim.fastq.gz)
    sam_file="$BAM_OUTPUT_DIR/$base_name.sam"
    bam_file="$BAM_OUTPUT_DIR/$base_name.bam"
    sorted_bam="$BAM_OUTPUT_DIR/$base_name.sorted.bam"
    novel_splice="${BAM_OUTPUT_DIR}/${base_name}.novel_splicesites.txt"
    
    echo "Processing $input_file..."
    
    # Запуск HISAT2
    hisat2 -p "$THREADS" -x "$HT_INDEX" -U "$input_file" -S "$sam_file" \
        --known-splicesite-infile "$SPLICE_SITES" \
        --novel-splicesite-outfile "$novel_splice" \
        --dta-cufflinks \
        --rna-strandness R \
        --summary-file "$BAM_OUTPUT_DIR/${base_name}.summary.txt"
    
    # Обработка BAM файлов
    samtools view -@ "$THREADS" -bSo "$bam_file" "$sam_file"
    samtools sort -@ "$THREADS" "$bam_file" -o "$sorted_bam"
    samtools index "$sorted_bam"
    
    # Удаление промежуточных файлов
    rm "$sam_file"
    rm "$bam_file"
    
    echo "Finished processing $input_file"
    echo "-----------------------------------"
done

echo "All files processed successfully"

# Подготовка к подсчету ридов
echo "Preparing for read counting..."
ls "$BAM_OUTPUT_DIR"/*.sorted.bam > "$BAM_OUTPUT_DIR/bam_list.txt"
GTF_FILE="./Homo_sapiens.GRCh38.101.gtf"

# Подсчет ридов с помощью Rsubread
echo "Running featureCounts in R..."
Rscript --vanilla - << 'EOF'
# Чтение параметров из переменных среды
bam_dir <- Sys.getenv("BAM_OUTPUT_DIR")
counts_dir <- Sys.getenv("COUNTS_OUTPUT_DIR")
gtf_file <- Sys.getenv("GTF_FILE")
threads <- as.numeric(Sys.getenv("THREADS"))

# Чтение списка BAM файлов
bam_files <- readLines(file.path(bam_dir, "bam_list.txt"))

# Подсчет ридов
library(Rsubread)
counts <- featureCounts(bam_files,
                       isGTFAnnotationFile = TRUE,
                       annot.ext = gtf_file,
                       GTF.featureType = "exon",
                       GTF.attrType = "gene_id",
                       strandSpecific = 2,
                       countMultiMappingReads = FALSE,
                       allowMultiOverlap = FALSE,
                       nthreads = threads)

# Чтение аннотации генов
library(rtracklayer)
gene_annot <- readGFF(gtf_file)
gene_annot <- gene_annot[gene_annot$type == "gene", c("gene_id", "gene_name", "gene_biotype")]

# Объединение результатов
result <- merge(gene_annot, counts$counts, by.x = "gene_id", by.y = "row.names", all = TRUE)

# Сохранение результатов
write.table(result, 
            file = file.path(counts_dir, "counts.txt"), 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE)

# Сохранение сырых счетчиков
write.table(counts$counts,
            file = file.path(counts_dir, "raw_counts.txt"),
            quote = FALSE,
            sep = "\t",
            col.names = NA)

# Сохранение статистики
write.table(counts$stat,
            file = file.path(counts_dir, "counts_stats.txt"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
EOF

echo "Analysis completed successfully. Results saved to $COUNTS_OUTPUT_DIR"



cd /media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/Fastq/

for prefix in $(ls *.fastq.gz | cut -d '_' -f 1 | sort -u); do
  cat "${prefix}"_*.fastq.gz > "${prefix}.fastq.gz"
done



for file in *.fastq.gz; do
    echo -n "$file: " >> read_counts.txt
    zcat "$file" | wc -l | awk '{print $1/4}' >> read_counts.txt
done


# Loop through each FASTQ file and run Trimmomatic
for input_file in *.fastq.gz; do
    # Define output filename (replace .fastq.gz with .trimmed.fastq.gz)
    output_file="/media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/Fastq_trimmed/${input_file%.fastq.gz}.trimmed.fastq.gz"
    
    # Run Trimmomatic and save output to log
    echo "Processing $input_file..." >> trimmomatic.log
    java -jar /usr/share/java/trimmomatic.jar SE \
        -threads 10 \
        "$input_file" \
        "$output_file" \
        ILLUMINACLIP:/media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/adapters.fa:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:15 >> trimmomatic.log 2>&1
    echo "-----------------------------------" >> trimmomatic.log
done

# reformat trimmomatic log
grep -E '^Processing|^Input Reads' /media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/Fastq/trimmomatic.log > /media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/Fastq/trimmed_summary.tsv


# prepare mature homo sapiens fasta
awk '/^>/ && /Homo sapiens/ {print; getline; gsub(/U/, "T"); print}' /media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/mature.fa > /media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/HS.mature.fa

hisat2-build /media/maxpauel/6372EF8A3B0879C1/miRNA/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.miRNAmasked_genome.fasta /media/maxpauel/6372EF8A3B0879C1/miRNA/hisat2_index_primiRNAinGenome/index
bowtie-build /media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/HS.mature.fa /media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/HS.mature_index/index


ht_index='/media/maxpauel/6372EF8A3B0879C1/miRNA/hisat2_index_primiRNAinGenome/index'
path_in='/media/maxpauel/01B5FA3D077CA3E9/miRNA_analysis/Fastq_trimmed/reads_le30/'
path_out='/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/'
for i in $path_in/*.fastq.gz; do {
res_file=${i%.fastq.gz*}.sam
res_file=${res_file/$path_in/}
res_file=$path_out$res_file

hisat2 -p 12 -x $ht_index -U $i  -S $res_file \
--rna-strandness F \
--known-splicesite-infile /media/maxpauel/01B5FA3D077CA3E9/Dry_immersion/GRCh38.104.splicesites.txt \
--summary-file ${res_file%.fastq.gz*}.txt;

samtools view -@ 12 -bSo ${res_file%.sam*}.bam $res_file;
samtools sort -@ 12 ${res_file%.sam*}.bam -o ${res_file%.sam*}.sorted.bam;
samtools index ${res_file%.sam*}.sorted.bam;
rm $res_file;
rm ${res_file%.sam*}.bam;
} ; done

samtools merge -@ 12 /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/merged.bam /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/*.bam
samtools sort -@ 12 -o /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/merged.sorted.bam  /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/merged.bam
samtools index /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/merged.sorted.bam
rm  /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/merged.bam



cd /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/
# extract unmapped reads
for bam in /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/*.bam; do
    # Extract base name (e.g., "sample" from "sample.bam")
    base="${bam%.bam}"
    
    # Extract unmapped reads and convert to FASTQ.GZ
    samtools view -@ 10 -b -f 4 "$bam" | \
    samtools fastq - | \
    gzip > "${base}_unmapped.fastq.gz"
done

mv /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_primary_miRNA/*.fastq.gz /media/maxpauel/6372EF8A3B0879C1/miRNA/fastq_unmapped_1
mv /media/maxpauel/6372EF8A3B0879C1/miRNA/fastq_unmapped_1/reads_gt30/*.fastq.gz /media/maxpauel/6372EF8A3B0879C1/miRNA/fastq_unmapped_1

cd /media/maxpauel/6372EF8A3B0879C1/miRNA/fastq_unmapped_1/

for prefix in $(ls *.fastq.gz | cut -d '_' -f 1 | sort -u); do
  cat "${prefix}"_*.fastq.gz > "${prefix}.fastq.gz"
done


cd /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1

ht_index='/media/maxpauel/01B5FA3D077CA3E9/Dry_immersion/hisat2/index'
path_in='/media/maxpauel/6372EF8A3B0879C1/miRNA/fastq_unmapped_1/'
path_out='/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1/'
for i in $path_in/*.fastq.gz; do {
res_file=${i%.fastq.gz*}.sam
res_file=${res_file/$path_in/}
res_file=$path_out$res_file

hisat2 -p 12 -x $ht_index -U $i  -S $res_file \
--known-splicesite-infile /media/maxpauel/01B5FA3D077CA3E9/Dry_immersion/GRCh38.104.splicesites.txt \
#--rna-strandness F \
--summary-file ${i%.fastq.gz*}.txt;

samtools view -@ 10 -bSo ${res_file%.sam*}.bam $res_file;
samtools sort -@ 10 ${res_file%.sam*}.bam -o ${res_file%.sam*}.sorted.bam;
rm $res_file;
rm ${res_file%.sam*}.bam;
} ; done

samtools merge -@ 12 /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1/merged.bam /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1/*.bam
samtools sort -@ 12 -o /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1/merged.sorted.bam  /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1/merged.bam
rm  /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1/merged.bam

for bam in /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1/*.bam; do
    samtools index "$bam"
done


# Input and output directories
bam_dir="/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1"
# Count primary alignments for each BAM
for bam in "$bam_dir"/*.bam; do
    count=$(samtools view -F 2304 -c "$bam")
    echo "$(basename "$bam"): $count primary alignments"
done > /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_1/alignment_counts.txt


# alignment to primary miRNA
for fq in /media/maxpauel/6372EF8A3B0879C1/miRNA/fastq_unmapped_2/*.fastq.gz; do
    # Extract base filename (e.g., "sample1" from "sample1.fastq.gz")
    base=$(basename "$fq" .fastq.gz)
    
    # Run Bowtie
    #bowtie -q -n 2 --best --strata -p 10 -m 1 -S  /media/maxpauel/6372EF8A3B0879C1/miRNA_old/Homo_sapiens_bowtie/Homo_sapiens.GRCh38.dna.primary_assembly \
    #bowtie -q -n 2 --best --strata -p 10 -a -S /media/maxpauel/6372EF8A3B0879C1/miRNA/bowtie_index_primiRNA/index \
    bowtie -q -n 2 --best -p 10 -S  /media/maxpauel/6372EF8A3B0879C1/miRNA_old/Homo_sapiens_bowtie/Homo_sapiens.GRCh38.dna.primary_assembly \
        <(zcat "$fq") \
        -S "/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_2/${base}.sam" \
        2> "/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_2/${base}.summary.txt"
done

# sam to sorted bam
for sam in /media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_2/*.sam; do
    base=$(basename "$sam" .sam)
    
    # 1. Convert SAM → BAM
    samtools view -@ 10 -b -o "/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_2/${base}.bam" "$sam" && \
    
    # 2. Sort BAM (required for indexing)
    samtools sort -@ 10 -o "/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_2/${base}.sorted.bam" "/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_2/${base}.bam" && \
    
    # 3. Index sorted BAM
    samtools index "/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_2/${base}.sorted.bam" && \
    
    # 4. Cleanup: Remove unsorted BAM and SAM
    rm -v "/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_2/${base}.bam" "$sam"  
done


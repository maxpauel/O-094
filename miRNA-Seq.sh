
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



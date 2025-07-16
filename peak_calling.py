import pysam

def find_local_max_peaks(bam_path, output_bed, min_max_cov=144, min_relative_cov=0.5, min_len=16, max_len=80):
    """
    Detect peaks by:
    1. Finding local maxima (coverage >= min_max_cov).
    2. Extending peaks to include bases with coverage >= (min_relative_cov * local max).
    3. Filtering by length (min_len <= peak length <= max_len).
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    peaks = []
    peak_counter = 1  # To number peaks
    
    for chrom in bam.references:
        # Initialize coverage arrays for both strands
        coverage_plus = [0] * (bam.get_reference_length(chrom) + 2)
        coverage_minus = [0] * (bam.get_reference_length(chrom) + 2)
        
        # Count per-base coverage on each strand
        for read in bam.fetch(chrom):
            if read.is_reverse:
                for pos in range(read.reference_start, read.reference_end):
                    coverage_minus[pos] += 1
            else:
                for pos in range(read.reference_start, read.reference_end):
                    coverage_plus[pos] += 1
        
        # Detect peaks on each strand separately
        for strand, cov_array in [("+", coverage_plus), ("-", coverage_minus)]:
            # Step 1: Find local maxima (>= min_max_cov)
            local_maxima = []
            for pos in range(1, len(cov_array) - 1):
                cov = cov_array[pos]
                if (cov >= min_max_cov and
                    cov > cov_array[pos - 1] and
                    cov >= cov_array[pos + 1]):  # Allows plateaus
                    local_maxima.append((pos, cov))
            
            # Step 2: Extend peaks around each local maximum
            for max_pos, max_cov in local_maxima:
                # Find left boundary (where coverage drops below 50% of max)
                left = max_pos
                while left > 0 and cov_array[left - 1] >= min_relative_cov * max_cov:
                    left -= 1
                
                # Find right boundary (where coverage drops below 50% of max)
                right = max_pos
                while right < len(cov_array) - 1 and cov_array[right + 1] >= min_relative_cov * max_cov:
                    right += 1
                
                peak_length = right - left + 1
                if min_len <= peak_length <= max_len:
                    peaks.append((
                        chrom,
                        left,
                        right + 1,  # BED uses 0-based start, 1-based end
                        f"peak_{peak_counter}",
                        max_cov,
                        strand
                    ))
                    peak_counter += 1
    
    # Write to BED file
    with open(output_bed, "w") as bed_file:
        for peak in peaks:
            bed_file.write(f"{peak[0]}\t{peak[1]}\t{peak[2]}\t{peak[3]}\t{peak[4]}\t{peak[5]}\n")
    
    print(f"Detected {len(peaks)} peaks. Saved to {output_bed}")

# Example usage
if __name__ == "__main__":
    find_local_max_peaks(
        bam_path="/media/maxpauel/6372EF8A3B0879C1/miRNA/bam_genome_PC_1/merged.sorted.bam",
        output_bed="/media/maxpauel/6372EF8A3B0879C1/miRNA/SEACR/peaks.bed",
        min_max_cov=144,       # Local max must have coverage >= 144
        min_relative_cov=0.5,   # Include bases with coverage >= 50% of local max
        min_len=16,             # Minimum peak width
        max_len=1000              # Maximum peak width
    )

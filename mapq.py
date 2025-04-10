import os
import sys
import pysam
from collections import Counter

def clean_sam_file(input_sam, output_sam):
    """
    Removes duplicate headers and writes a cleaned SAM file.
    """
    with open(input_sam, 'r') as infile:
        headers = []
        body = []
        seen = set()

        for line in infile:
            if line.startswith('@'):
                if line not in seen:
                    headers.append(line)
                    seen.add(line)
            else:
                body.append(line)

    with open(output_sam, 'w') as outfile:
        outfile.writelines(headers + body)

def compute_cumulative_mapq(sam_path):
    """
    Reads a SAM file, computes cumulative MAPQ distribution.
    Returns a dict {MAPQ: cumulative_count}
    """
    samfile = pysam.AlignmentFile(sam_path, "r")
    mapq_values = [read.mapping_quality for read in samfile if not read.is_unmapped]
    samfile.close()

    mapq_dist = Counter(mapq_values)

    # Cumulative count in descending MAPQ
    cumulative_counts = {}
    total = 0
    for mapq in sorted(mapq_dist, reverse=True):
        total += mapq_dist[mapq]
        cumulative_counts[mapq] = total

    return cumulative_counts

def process_sam_file(sam_file_path, output_summary_file):
    """
    Processes a given SAM file and writes cleaned version and MAPQ summary.
    """
    if not os.path.isfile(sam_file_path):
        print(f"? Error: '{sam_file_path}' does not exist or is not a file.")
        sys.exit(1)

    folder_path = os.path.dirname(sam_file_path)
    cleaned_sam = os.path.join(folder_path, "aligned_output_unique.sam")

    print(f"\n?? Processing SAM file: {sam_file_path}")
    clean_sam_file(sam_file_path, cleaned_sam)
    cumulative_mapq = compute_cumulative_mapq(cleaned_sam)

    with open(output_summary_file, "w") as f:
        f.write("MAPQ\tCumulative_Count\n")
        for mapq in sorted(cumulative_mapq):
            f.write(f"{mapq}\t{cumulative_mapq[mapq]}\n")

    print(f"? Cleaned SAM saved to: {cleaned_sam}")
    print(f"? MAPQ summary saved to: {output_summary_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 mapq.py <path_to_sam_file> <output_summary_file>")
        sys.exit(1)

    sam_file = sys.argv[1]
    output_file = sys.argv[2]
    process_sam_file(sam_file, output_file)

import sys, time

from rich.progress import Progress
from Bio import SeqIO

read_file = sys.argv[1]  # you should probably use `click`!

with open(read_file) as reads:
    # Count all records in the file so we know to track
    # overall progress (FASTQ has 4 lines per record).
    num_reads = sum(1 for _ in open(read_file)) // 4

    # rewind the handle to the beginning of the file.
    reads.seek(0)

    # process the reads (with progress!)
    with Progress() as progress:
        task = progress.add_task(
            "Loading reads", total=num_reads, unit="reads", width=40
        )
        for read in SeqIO.parse(reads, "fastq"):
            time.sleep(0.1)  # "work"
            print(read.seq)
            progress.update(task, advance=1)

print(f"Processed {num_reads} reads.")

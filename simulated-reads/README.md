# Simulated data sets

All simulated data sets have been generated from [an isolate of COVID-19 on 
NCBI].

You can find this isolate in this directory as `covid-genome.fasta`.

[an isolate of COVID-19 on NCBI]: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=genbank

## Simulated error free sets

1. `covid.fq`: This is the entire COVID genome as one giant "read".
2. `covid.fasta_simulated-allforward-errorfree.fq`: This has simulated reads;
   all are 250bp long and all are in the forward orientation (relative to the
   original genome file).
3. `covid.fasta_simulated-forwardreverse-errorfree.fq`: This has simulated
   reads; all are 250bp long and are in both the forward and reverse strand.

## Simulated with errors sets

1. `covid.fasta_simulated-errors-tips.fq`: errors that *should* appear as tip
   structures with $k=37$ (all errors are within 35bp of the ends of a read).
2. `covid.fasta_simulated-errors-bubbles.fq`: errors that *should* appear as
   bubble structures with $k=37$ (all errors are between 35bp from the start and
   35bp from the end of a read).
3. `covid.fasta_simulated-errors-tips-bubbles.fq`: errors that appear anywhere
   in the read.

<p align="center">
  <src="https://github.com/albertozeni/gelato/blob/main/media/gelato.png">
</p>

## Introduction
<p align="justify">
Short-read WGS is a critical tool for clinical applications and basic research.
Indeed, given their high level of precision and relatively low cost, short read technologies are still the foundation of the variant calling and assembly processes.
However, characterizing the variants from a reference genome using short sequences only reduces the performance of the process in both recall and precision, given their repetitive and fragmented nature.
Variant calling can also be performed by exploiting assembled genomes in order to better represent the discovered variants.
Assembling genomes with short reads often leads to low quality results, requiring sequences attained using other technologies (e.g., ONT and HiFi) or a reference genome in order to support the process in challenging regions.
In this context, we present Gelato, a graph-based genome assembler that exploits multiple-references together abd a graph k-mer structure to rapidly and accurately assemble genomes.
To benchmark short-read variant calling, we employed 141 diverse clinical Mycobacterium tuberculosis (Mtb) isolates
sequenced with Illumina short-reads, PacBio/ONT long-reads.
We systematically evaluated short-read variant calling accuracy comparing it to multiple Gelato generated genomes employing support from similar isolates or the Mtb reference genome only, and optional trusted assemblies.
</p>

## Results
<p align="justify">
We evaluate Gelato building a total of 141 tuberculosis genomes, and systematically evaluate their variant calling results against short-read based methodologies.
Reference-based Illumina variant calling demonstrated a maximum recall of 89.4% and minimum precision
of 98.5% across the evaluated datasets, while Gelato based genomes attained a maximum recall of 98.30% with minimum precision of 96.6%.
To reach short-read levels of precision we masked called variants for Gelato genomes only to regions with a high number of mapping reads, attaining a minimum precision of 98.7% while still maintaining a significantly elevated level of recall of 97.3%.
Our benchmarking results have broad implications for the use of \ac{WGS} and Hifi ONT sequences in variant calling, showcasing how short-reads based variant calling still remains a practical methodology for the foreseeable future.
We present a method that, in combination with high-quality reference assemblies, offers a powerful approach for genotyping that can effectively complement state-of-the-art variant calling pipelines. 
</p>

## Usage

### Compilation

Gelato requires C++17, libz and pthread for compilation.
It also exploits robinhood for hash functionalities. To compile Gelato simply type:
```
make -j4
```
This will create a bin/ folder with the gelato executable.

### Functionalities

Gelato exploits multiple reference datasets and short reads to assemble a genome that produces high accuracy and sensitivity when doing the variant calling process.

To explore the various inputs and outputs gelato can produce simply type:
```
gelato -h
```
This will output the following menu:
```
Usage: gelato [options]
Options:
  -h               help
  -v               verbose prints
  -N               count kmers containing the N character (default do not count)
  -s STR           save reads kmer counting info
  -S STR           save references kmer counting info
  -k INT           kmer len used for counting and graph building (max 32)
  -r STR           single reference input file
  -R STR           list of reference input files (on per line)
  -i STR           single reads input file
  -I STR           list of reads input files (on per line)
  -a STR           additional assembly input file
  -A STR           list of additional assembly input files (on per line)
  -n INT           num of times a kmer has to appear to not be considered error
  -p INT           set n. of threads to run gelato with (default 32)
  -c INT           set read coverage (default 1)
  -g STR           output filename for graph (gfa format and extension are automatically applied)
  -o STR           output filename for genome (fasta format and extension are automatically applied)
  -O STR           output filename for genome (fastq format and extension are automatically applied)
Note that gelato always counts for canonical k-mers (the count is the combined number of occurrences of both a k-mer and it reverse complement)
```

Gelato requires at least one reference genome file, and one file of reads to assemble a genome.
It than counts the kmers from all references and reads and computes the similarity values of the references, choosing the reference which is closer to the reads.
It than creates a graph using the reads' kmers and, following the path of the chosen reference genome, creates a new genome using both the information of the newly created reads kmer graph and the chosen reference.

One typical command to run gelato would be something like this:
```
gelato -I list-reads.txt -R list-refs.txt -A list-assemblies.txt -p 16 -k 31 -o assembly_best
```
This will use a list of reads datasets (paired- or single-end) called ```list-reads.txt```, this contains the paths to a fasta or fastq file (```.gz``` files are supported as well), one per line.
In case only one file is present, the user could use the ```-i``` option and directly state the path of the reads file.
The command will also use a set of references, called ```list-refs.txt``` which is a set of datasets containing the reference files wanted to be used as base by the assembler, the content of the file follows the logic of the reads file (the same goes for the single reference usage by indicating it with ```-r```).
The above command exploits also an additional feature of Gelato, by also adding the support of trusted assemblies from the reads we are selecting with the ```-A``` option. This improves the performance of the outputted Gelato dataset if the indicated assemblies are of high quality.
One use case might be to build an assembly using ```Spades```, for example, and using the outputted assembly as additional info for Gelato.
This command also specifies the number of threads and kmer length and outputs a genome called ```assembly_best```, in the ```fasta``` format.

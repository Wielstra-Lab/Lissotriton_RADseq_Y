_Notes on identification of male-linked presence/absence markers in Lissotriton vulgaris_

### Software versions: ###

**Stacks:** v2.54  
**Perl:** v5.36.1  
**Samtools:** v1.16.1  
**R:** v4.3.2  
**BLAST+:** v2.13.0  
**Primer 3:** v2.6.1  

### Scripts used: ###

RADcov.pl  
RADcoverage.sh  
RADsort_sexed_1.R  
RADsort_sexed_2.R  
fasta_to_primer_IO_2.sh 

### Files used: ###

PopMap_Liss_sexed_1.txt  
IO_header_2  
Liss_Y_1_table.txt   
combined_liss_Y.fa  

SECTION 1: Upstream RADseq bioinformatics
-----------------------------------------
> The denovo_map.pl program from the stacks package was run on the sexed data set 3 times with different values for parameter M (2, 6 and 10)

```sh
denovo_map.pl -T 16 --samples ~/data/Liss_RADseq/Sexed_samples --popmap PopMap_Liss_sexed_1.txt -o ~/data/Liss_RADseq/Sexed_stacks_1 -M 2 --paired
denovo_map.pl -T 16 --samples ~/data/Liss_RADseq/Sexed_samples --popmap PopMap_Liss_sexed_1.txt -o ~/data/Liss_RADseq/Sexed_stacks_2 -M 6 --paired
denovo_map.pl -T 16 --samples ~/data/Liss_RADseq/Sexed_samples --popmap PopMap_Liss_sexed_1.txt -o ~/data/Liss_RADseq/Sexed_stacks_3 -M 10 --paired
```

SECTION 2: Identification of male-linked markers based on coverage
------------------------------------------------------------------

### Generate initial candidate markers ###

> First the coverage of every marker in each sample was calculated using RADcov.pl (which calls RADcoverage.sh).

```sh
perl RADcov.pl -d ~/data/Liss_RADseq/Sexed_stacks_1
perl RADcov.pl -d ~/data/Liss_RADseq/Sexed_stacks_2
perl RADcov.pl -d ~/data/Liss_RADseq/Sexed_stacks_3
```
> Coverage was then filtered to produce an inital set of candiate markers using the RADsort_sexed_1.R script

```sh
Rscript RADsort_sexed_1.R PopMap_Liss_sexed_1.txt ~/data/Liss_RADseq/Sexed_stacks_1/Coverage Liss_Y_1 
Rscript RADsort_sexed_1.R PopMap_Liss_sexed_1.txt ~/data/Liss_RADseq/Sexed_stacks_2/Coverage Liss_Y_2 
Rscript RADsort_sexed_1.R PopMap_Liss_sexed_1.txt ~/data/Liss_RADseq/Sexed_stacks_3/Coverage Liss_Y_3 
```

> The lists of candidate markers names were then used to build fasta files from the catalog.fa files created by stacks in each run

 ```sh
cat Liss_Y_1_markers.txt | while read marker do; grep -A 1 ">"$marker" " ~/data/Liss_RADseq/Sexed_stacks_1/catalog.fa; done > Y_markers_Liss_1.fa
cat Liss_Y_2_markers.txt | while read marker do; grep -A 1 ">"$marker" " ~/data/Liss_RADseq/Sexed_stacks_2/catalog.fa; done > Y_markers_Liss_2.fa
cat Liss_Y_3_markers.txt | while read marker do; grep -A 1 ">"$marker" " ~/data/Liss_RADseq/Sexed_stacks_3/catalog.fa; done > Y_markers_Liss_3.fa
```
### Filter for paralogs and select top 10 markers ###

> The catalog.fa file created by stacks was then converted into a BLAST database

 ```sh
cp ~/data/Liss_RADseq/Sexed_stacks_1/catalog.fa ~/data/Liss_RADseq/Sexed_catalogs/catalog_1/catalog_1.fa
makeblastdb -dbtype nucl -in catalog_1.fa
cp ~/data/Liss_RADseq/Sexed_stacks_2/catalog.fa ~/data/Liss_RADseq/Sexed_catalogs/catalog_2/catalog_2.fa
makeblastdb -dbtype nucl -in catalog_2.fa
cp ~/data/Liss_RADseq/Sexed_stacks_3/catalog.fa ~/data/Liss_RADseq/Sexed_catalogs/catalog_3/catalog_3.fa
makeblastdb -dbtype nucl -in catalog_3.fa
```

> The markers are then BLASTed against their own catalog.fa to detect paralogs

```sh
blastn -query Y_markers_Liss_1.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_1/catalog_1.fa -outfmt 6 -perc_identity 80 -qcov_hsp_perc 25 > ~/data/Liss_RADseq/Sexed_catalogs/1_1_matches.txt
blastn -query Y_markers_Liss_2.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_2/catalog_2.fa -outfmt 6 -perc_identity 80 -qcov_hsp_perc 25 > ~/data/Liss_RADseq/Sexed_catalogs/2_2_matches.txt
blastn -query Y_markers_Liss_3.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_3/catalog_3.fa -outfmt 6 -perc_identity 80 -qcov_hsp_perc 25 > ~/data/Liss_RADseq/Sexed_catalogs/3_3_matches.txt
```

> The RADsort_sexed_2.R script is then used to select a top 10 markers after filtering by paralogs

```sh
Rscript RADsort_sexed_2.R  Liss_Y_1_table.txt ~/data/Liss_RADseq/Sexed_catalogs/1_1_matches.txt Liss_Y_1_final
Rscript RADsort_sexed_2.R  Liss_Y_2_table.txt ~/data/Liss_RADseq/Sexed_catalogs/2_2_matches.txt Liss_Y_2_final
Rscript RADsort_sexed_2.R  Liss_Y_3_table.txt ~/data/Liss_RADseq/Sexed_catalogs/3_3_matches.txt Liss_Y_3_final
```
### Intergrate final set from all 3 runs ###

> The markers are also BLASTed against the catalog.fa files from the other stacks runs, to allow for integration of the results (only the top hit is required)

```sh
blastn -query Y_markers_Liss_1.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_2/catalog_2.fa -outfmt 6 -perc_identity 90 -qcov_hsp_perc 50 | sort -k1,1 -u 
> ~/data/Liss_RADseq/Sexed_catalogs/1_2_matches.txt
blastn -query Y_markers_Liss_1.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_3/catalog_3.fa -outfmt 6 -perc_identity 90 -qcov_hsp_perc 50 | sort -k1,1 -u
> ~/data/Liss_RADseq/Sexed_catalogs/1_3_matches.txt
blastn -query Y_markers_Liss_2.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_1/catalog_1.fa -outfmt 6 -perc_identity 90 -qcov_hsp_perc 50 | sort -k1,1 -u 
> ~/data/Liss_RADseq/Sexed_catalogs/2_1_matches.txt
blastn -query Y_markers_Liss_2.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_3/catalog_3.fa -outfmt 6 -perc_identity 90 -qcov_hsp_perc 50 | sort -k1,1 -u 
> ~/data/Liss_RADseq/Sexed_catalogs/2_3_matches.txt
blastn -query Y_markers_Liss_3.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_1/catalog_1.fa -outfmt 6 -perc_identity 90 -qcov_hsp_perc 50 | sort -k1,1 -u 
> ~/data/Liss_RADseq/Sexed_catalogs/3_1_matches.txt
blastn -query Y_markers_Liss_3.fa -db ~/data/Liss_RADseq/Sexed_catalogs/catalog_2/catalog_2.fa -outfmt 6 -perc_identity 90 -qcov_hsp_perc 50 | sort -k1,1 -u 
> ~/data/Liss_RADseq/Sexed_catalogs/3_2_matches.txt
```

SECTION 3: Primer design
------------------------
> A combined fasta file was made of all selected candidate makers, where the same marker was found in multiple runs the sequence from the run with the highest value of M was chosen
> This fasta file was then used to design primers with Primer3, two sets of primers were desinged for each marker (amplifying a short and long primer each)
> The input file for primer 3 was created by running the fasta_to_primer_IO_2.sh script (the IO_header_2 contains the common settings for all primer pairs)

```sh
sh fasta_to_primer_IO_2.sh  Liss_Y_1_table.txt combined_liss_Y.fa IO_header_2 Liss_IO_2
```
> The resultant Liss_IO_2 file was then used as the input for primer 3

```sh  
~/primer3-2.6.1/src/primer3_core Liss_IO_2 > Liss_IO_2.out
```
> The primer designs from this output file were then ordered from IDT and tested for sex specificity via PCR

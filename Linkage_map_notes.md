_Notes on construction of L. vulgaris RADseq linkage map_


### Software versions: ###

**Stacks:** v2.54  
**Perl:** v5.36.1  
**Samtools:** v1.16.1  
**R:** v4.3.2  
**BLAST+:** v2.13.0   
**Lep-MAP3:** v0.5  

SECTION 1: Upstream RADseq bioinformatics
-----------------------------------------

> The denovo_map.pl program from the Stacks package is run on the fastq.gz files from all linkage maps samples.  
> The populations program from Stacks was then run with the --vcf flag to generate a joint vcf file for the linkage family.
```sh
denovo_map.pl -T 16 --samples ~/data/Liss_RADseq/Linkage_samples --popmap PopMap_Liss_Linkage_1.txt -o ~/data/Liss_RADseq/Linkage_stacks --paired
populations -P ~/data/Liss_RADseq/Linkage_stacks --popmap PopMap_Liss_Linkage_1.txt -t 16 --vcf
````

SECTION 2: Marker filtering
---------------------------

### Filter the vcf file for SNPs to be placed on the linkage map ###

> The joint vcf was strictly filtered for minor allele frequency, mean depth and missing data. Indels were removed and 1 SNP per marker was selected

```sh
vcftools --vcf ~/data/Liss_RADseq/Linkage_stacks/lissotriton.snps.vcf --recode --recode-INFO-all --out Lissotriton_joint.filtered --maf 0.2 --min-meanDP 10 --max-missing 0.95 --remove-indels --thin 500
````

### Create genotype calls for candidate Y-linked presence/absence markers ###

> A BLAST database was created from the catalog.fa file created by denovo_map.pl

````sh
cp ~/data/Liss_RADseq/Linkage_stacks/catalog.fa Blast/Liss_linkage_catalog.fa
makeblastdb -dbtype nucl -in Liss_linkage_catalog.fa
````

> The candidate Y-linked markers identified in the know-sex adults were blasted against the linkage map catalog

```sh
blastn -qcov_hsp_perc 30 -outfmt 6 -perc_identity 95 -query Liss_sexed_Y_candidates_1.fa -db Blast/Liss_linkage_catalog.fa | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > Y_candidates_in_linkage_map.blast
````

> The coverage of every marker in each sample in the linkage family was calculated using RADcov.pl (which calls RADcoverage.sh)

```sh
perl RADcov.pl -d ~/data1/Liss_RADseq/Linkage_stacks
````

> Coverage_sort_1.R was then used to produce a table of coverage in each sample for the candidate Y-linked markers

```sh
Rscript Coverage_sort_1.R PopMap_Liss_Linkage_1.txt ~/data1/Liss_RADseq/Linkage_stacks Y_candidates_in_linkage_map.blast sexed_Y_candidates_raw.cov
```

> This was then filtered with Coverage_filter_1.R, to eliminate markers with very low coverage, or coverage in the mother

```sh
Rscript Coverage_filter_1.R sexed_Y_candidates_raw.cov sexed_Y_candidates_filtered.cov 1 17011c
```

> Add_to_call_table_sexed.R was used to transform the coverage table into Psuedo-SNP calls in a format compatable with LepMAP 3

```sh
Rscript Add_to_call_table_sexed.R sexed_Y_candidates_raw.cov Y_marker_calls_1.txt Y
```

SECTION 3: Linkage map construction
-----------------------------------

### Identify Linkage Groups ###

> The LepMAP 3 package was used to create a linkage map - starting with the ParentCall 2 program

```sh
java -cp ~/LepMap/bin/ParentCall2 data = Liss_RAD.ped vcfFile = Lissotriton_joint.filtered.recode.vcf > Lissotriton_RAD_1_parent.call
```

> The resulting .call file was then concatenated with the output of the Add_to_call_table_sexed.R script to include the candiate Y-linked presence/absence markers

```sh
cat Lissotriton_RAD_1_parent.call Y_marker_calls_1.txt > Lissotriton_RAD_1_Y_marker_parent.call
```

> The next 2 steps of the the LepMAP 3 pipeline (SeperateChromosomes 2 and JoinSingles 2) were run

```sh
java -cp ~/LepMap/bin/ SeparateChromosomes2  data = Lissotriton_RAD_1_Y_marker_parent.call lodLimit = 20 distortionLod = 1 > Lissotriton_RAD_Y_map_1.txt
java -cp ~/LepMap/bin/ JoinSingles2All  data = Lissotriton_RAD_1_Y_marker_parent.call map = Lissotriton_RAD_Y_map_1.txt lodLimit = 15 > Lissotriton_RAD_Y_map_with_singles_1.txt
```

> Marker distribution across linkage groups was checked with the following command:
```sh
cut -f 1 Lissotriton_RAD_Y_map_with_singles_1.txt|sort|uniq -c|sort -n
```
> Giving the following result, the first column is number of markers per group, the second is the group number, group 0 represents unmapped markers:

        2 14
        2 15
        3 13
      386 12
      553 11
      561 10
      581 9
      763 8
      896 7
      912 6
     1023 5
     1162 4
     1199 3
     1309 2
     1349 1
     6000 0

> Groups 13-15 represent artifacts which were not analysed further

### Order markers and build maps ###

> Translate the marker designations from ParentCall 2 back to the input names

```sh
cat Lissotriton_RAD_1_Y_marker_parent.call|cut -f 1,2|awk '(NR>=7)' > snps.txt
```

> The OrderMarkers 2 program from LepMAP 3 was used to order the markers assigned to each linkage group, the groups were then grouped into final maps
> This script was run three times, creating both a sex averaged map and maternal and paternal maps

For sex averaged map:

```sh
input=Lissotriton_RAD_Y_map_with_singles_1.txt
mkdir Ordered_maps_sex_averaged

for i in {1..12}
do
  output="Ordered_maps_sex_averaged/Lissotriton_Y_map_1_group_"$i"_SA.txt"
  java -cp ~/LepMap/bin/ OrderMarkers2 map = $input data = Lissotriton_RAD_1_Y_marker_parent.call chromosome = $i  numMergeIterations = 20 numPolishIterations = 8 minError=0.02 scale = M/N 2 numThreads = 8 sexAveraged = 1 > $output
  tail -n +4 $output | while read marker position line; do echo -e $marker"\t"$position"\t"$i; done >> Ordered_maps_sex_averaged/Lissotriton_Y_map_1_SA_ordered.txt
done

awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt Ordered_maps_sex_averaged/Lissotriton_Y_map_1_SA_ordered.txt | cut -f 1,2,4,5 > Ordered_maps_sex_averaged/Lissotriton_Y_map_1_SA_mapped.txt
```

For paternal and maternal maps:

```sh
input=Lissotriton_RAD_Y_map_with_singles_1.txt
output_stem=Lissotriton_Y_map_1
outdir=Ordered_maps_sexed

mkdir $outdir

for i in {1..12}

do
  output_M=$outdir"/"$output_stem"_"$i"_M.txt"
  output_F=$outdir"/"$output_stem"_"$i"_F.txt"

  echo "working on "$output_M

  java -cp ~/LepMap/bin/ OrderMarkers2 map = $input data = Lissotriton_RAD_1_Y_marker_parent.call chromosome = $i  numMergeIterations = 20 numPolishIterations = 8 minError=0.02 scale = M/N 2 numThreads = 8 informativeMask=23 > $output_M
  tail -n +4 $output_M | while read marker pat_position mat_position line; do echo -e $marker"\t"$pat_position"\t"$mat_position"\t"$i; done >> $outdir"/"$output_stem"_male_ordered.txt"


  java -cp ~/LepMap/bin/ OrderMarkers2 map = $input data = Lissotriton_RAD_1_Y_marker_parent.call chromosome = $i  numMergeIterations = 20 numPolishIterations = 8 minError=0.02 scale = M/N 2 numThreads = 8 informativeMask=13 > $output_F
  tail -n +4 $output_F | while read marker pat_position mat_position line; do echo -e $marker"\t"$pat_position"\t"$mat_position"\t"$i; done >> $outdir"/"$output_stem"_female_ordered.txt"

done

awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt $outdir"/"$output_stem"_male_ordered.txt" | cut -f 1,2,4,5 > $outdir"/"$output_stem"_paternal_mapped.txt"

awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt $outdir"/"$output_stem"_female_ordered.txt" | cut -f 1,3,4,5 > $outdir"/"$output_stem"_maternal_mapped.txt"
```

> The sex-averaged map was re-ordered (and re-numbered) with linkage groups arranged by length (longest to shortest)

```sh
Rscript Trim_map_2.R Lissotriton_Y_map_1_SA_mapped.txt Lissotriton_Y_map_1_SA_mapped_reordered.txt NA 100
```

> The numbering and orientation of the sex-averaged map was then applied to the maternal and paternal maps

```sh
Rscript arrange_map.R Lissotriton_Y_map_1_SA_mapped_reordered.txt Lissotriton_Y_map_1_paternal_mapped.txt Lissotriton_Y_map_1_Pat_mapped_reordered.txt
Rscript arrange_map.R Lissotriton_Y_map_1_SA_mapped_reordered.txt Lissotriton_Y_map_1_maternal_mapped.txt Lissotriton_Y_map_1_Mat_mapped_reordered.txt
```

SECTION 4: Comparative genomics with _P. waltl_
-----------------------------------------------
> A BLAST database was created for the _P. waltl_ genome

```sh
makeblastdb -dbtype nucl -in aPleWal1.pri.20220803.fasta
```

> The the sequences placed on the linkage map were BLASTed against the _P. waltl_ database

```sh
blastn -query Lissotriton_Y_map_1_markers.fa  -db aPleWal1.pri.20220803.fasta -outfmt 6 -evalue 1e-20 -word_size 11 -num_threads 8 > Liss_Y_map_Blast_raw.txt

sort -k1,1 -k12,12nr Liss_Y_map_Blast_raw.txt > Liss_Y_map_Blast_sorted.txt
```

> The sorted results from the BLAST were then filtered with Filter_blast.R
```sh
Rscript Filter_blast.R Liss_Y_map_Blast_sorted.txt Liss_Y_map_Blast_filtered.txt
```

> The filtered results, along with the Linkage map and file containing information about the _P. waltl_ genome structure were used as a input for the construction of an Oxford plot
```sh
Rscript Plot_map_v_genome_1.R Lissotriton_Y_map_1_SA_mapped_reordered.txt Liss_Y_map_Blast_filtered.txt Pluro_genome_struct_chr.txt Lissotriton_Y_map_1_SA_v_Pluro
```

SECTION 5: Analysis
-------------------

> Visulisation of the linkage map was performed in R, highlighting the Y-linked markers identified in known_sex population

```sh
Rscript Draw_groups.R Lissotriton_Y_map_1_SA_mapped_reordered.txt Liss_Y_markers.txt Liss_Y_map_SA.pdf "_Lissotriton vulgaris_ Linkage Groups"
```

> A marker density plot was created for the paternal and maternal maps respectively

```sh
Rscript Plot_marker_density.R Lissotriton_Y_map_1_Pat_mapped_reordered.txt Liss_Y_markers.txt Lissotriton_Y_map_1_Mat_mapped_reordered.txt Liss_density_plot.pdf
```
> The number of paternal markers and SNPs was plotted per linkage group

```sh
Rscript Plot_paternal_markers.R Lissotriton_Y_map_1_SA_mapped_reordered.txt father_hetero_markers_filtered_by_mother.txt Paternal_specific_snps_marker_list.txt Liss_group_plots.pdf
```

> The occurence of paternal markers and SNPs was visulised in a manhatten style plot

```sh
Rscript Plot_paternal_probability.R Lissotriton_Y_map_1_SA_mapped_reordered.txt father_hetero_markers_filtered_by_mother.txt Paternal_specific_snps_marker_list.txt Liss_probability_plots.pdf
```

 

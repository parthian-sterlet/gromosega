# Gromosega
The Gromosega (GROups of MOtifs SElection by Genetic Algorithm) tool searches for groups of several enriched transcription factor (TF) binding site (TFBS) motifs respecting cooperative action of TFs. The annotated binding site motifs for known TFs are selected from species- or taxon-speific collections from the public databases HOCOMOCO for mammals ([Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077)) and JASPAR for vertebrates, insects and plants ([Ovek Baydar et al., 2026](https://doi.org/10.1093/nar/gkaf1209)). The enrichment of each motif group is deduced from the comparison of the input sets of the foreground / background DNA sequences, the promoters of differentially expressed genes (DEGs) / not DEGs for RNA-seq data.  The metrics [partial Area under Precision-Recall curve (pAUPRC)](https://github.com/parthian-sterlet/metarea) ([Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90)) is used to transit from the enrichment of single TFBS motifs to the enrichment of their groups. The approach of Genetic Algorithm (GA) is used to avoid the exhaustive search in a combinatorically large space of possible TFBS motif content of groups. The default tool application implies independent testing of certain series of motif group sizes, e.g. the groups of 1, 3, 5, ... 99 motifs are tested (from 1 to 99 motifs with the step of 2). The tool is primarely designed for the analysis of the whole transcriptome data, such as RNA-seq data, but it can be also be applied for any other DNA sequence data presuming the specific cooperative action of TFs in transcription regulation. 

# Description
Gromosega takes the basic computation principle from [MetArea tool](https://github.com/parthian-sterlet/metarea). For a motif representing a TFBS pattern, the recognition accuracy is calculated as the recognition performance metrics pAUPRC ([Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90); [Davis, Goadrich, 2006](https://doi.org/10.1145/1143844.1143874)). Consider several TFBS motifs compiled in a group. Let for this group and for given DNA sequence certain member motif has the best (most conservative) predicted site among all other motifs of this group with their predicted sites. Hence, this motif and the quality of its best site represent the whole group, i.e. different member motifs can stand for this group in various DNA sequences. Thus, for a group of motifs and two given sets of foreground and background sequences the one PR curve and its recognition performance metric pAUPRC value are computed as described earlier ([Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90). 

The genetic algorithm (GA) can help to define the most important TFBS motif groups. The GA defines its terms of individuals, the population of individuals, and the fitness function of each individual as follows. An individual is a group of several motifs; a set of about a couple hundred groups of the same number of motifs represents the population; the fitness function is the recognition performance metrics pAUPRC. Testing and optimization of the motif content for the population of many groups by GA allows to find and rank the best groups of motifs reflecting differences between two input DNA sequence sets. The participant motifs of a group can represent structurally different binding sites for the same TF and binding sites of different TFs acting together comsistently in gene transcription regulation (Levitsky et al., [2014](https://doi.org/10.1186/1471-2164-15-80), [2024](https://doi.org/10.18699/vjgb-24-90)), otherwise these motifs are combinations alternatively implementing distinct pathways of gene regulatory networks. Since the default option of TFs action is mutual cooperative interactions ([Morgunova & Taipale, 2017](https://doi.org/10.1016/j.sbi.2017.03.006)), the protein complexes of multiple TFs are reflected in in RNA-seq data where various DEGs contain certain more or less similar combinations of TFBSs for involved TFs. Hence, the functional relationship of several distinct motifs of the same group is recognized through increased recognition performance. For the fixed number of motifs, what is the condition of each individual GA run, the groups with the greater pAUPRC values are the better representatives of the functional TFs compared to the groups with the smaller pAUPRC values. Thus, Gromosega predicts TFBS motifs respecting TFs participating in distinct multiprotein complexes functioning in gene transcription regulation. Importanly, this approach identifies motifs of each group simultaneously and together, but not each of the motifs separately. Consequently, we can get the three main distinct options for TF binding _in vivo_. Namely, a target TF _in vivo_ can bind genomic DNA: (a) directly, there is a binding site of this target TF in DNA; (b) using another TF, sites of the target and another partner TF are near in DNA; (c) not directly, only the partner TF has a binding site in DNA. 

# Requirements
The source code is written in C++ language. To compile execubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Repository structure
Folder [**src**](https://github.com/parthian-sterlet/gromosega/tree/main/src) contains the [major](https://github.com/parthian-sterlet/gromosega/blob/main/src/minimax.cpp) and some supporting C++ source code files.  

Folder [**run**](https://github.com/parthian-sterlet/antinoise/tree/main/run) contains two command lines, implementing the tool pipeline for example RNA-seq data.

Folder [**examples**](https://github.com/parthian-sterlet/gromosega/tree/main/bin/examples) contains files required as the functional examples of the tool.

Folder [**genomes**](https://github.com/parthian-sterlet/antinoise/tree/main/genomes) contains whole genome annotations required for mentioned above functional examples. 

# How to compile
* In Linux system: 
```
git clone https://github.com/parthian-sterlet/gromosega
cd gromosega/run
chmod a+x build.sh
./build.sh
cd ../genomes
tar -xvzf protcod_dm6_m1kb_p100_tss_272tbp.tab.tar.gz
cd ../run

```
* In Windows system:

separate compilation of all source files in VC++

# Algorithm
Gromosega algorithm considers a pair of foreground/background sequence sets derived from RNA-seq/Microchip data. For these data the foreground and background sequences are promoters of DEGs and non-DEGs, that are defined by the default criteria for adjusted p-value (padj) & Fold Change (FC) thresholds padj<sub>CRIT</sub> & FC<sub>UP</sub> / FC<sub>DOWN</sub>: 
- for DEGs: padj < padj<sub>CRIT</sub>, & FC >= FC<sub>UP</sub> and FC <= FC<sub>DOWN</sub> for up- and down-regulated DEGs; 
- for non-DEGs: padj > padj<sub>CRIT</sub> &  FC<sub>DOWN</sub> <= FC <= FC<sub>UP</sub> .
The default value for the threshold padj<sub>CRIT</sub> is 0.05, threshold values FC<sub>UP</sub> and FC<sub>DOWN</sub> are selected to compile fixed number of DEGs/non-DEGs, e.g. 500 / 2000. Genes with the maximal /minimal values of |FC| are compiled in DEGs / non-DEGs.

## GA input data:
- a pair of foreground / background sets, containing N<sub>POS</sub> / N<sub>NEG</sub> DNA sequences;
- a collection of M<sub>TOT</sub> TFBS motifs derived from public databases like [Hocomoco v 14](https://hocomoco14.autosome.org/) ([Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077)) or [Jaspar2026](https://jaspar.elixir.no/) ([Ovek Baydar et al., 2026](https://doi.org/10.1093/nar/gkaf1209)), this collection can be complemented by motifs derived by _de novo_ motif search for an appropriate dataset, for each motif it is required its name (e.g. TF name and motif designation) and short description of its TF class/family according to the hierarchical classification of mammalian ([Wingender et al., 2018](https://doi.org/10.1093/nar/gkx987)) and plant ([Blanc-Mathieu et al., 2024](https://doi.org/10.1016/j.tplants.2023.06.023)) TFs by DNA-binding domain (DBD) structure;
- two matrices of the sizes M<sub>TOT</sub> x N<sub>POS</sub> and M<sub>TOT</sub> x N<sub>NEG</sub>. The rows of these matrices represent sequences (1 ≤ n ≤ N<sub>POS</sub>, 1 ≤ n ≤ N<sub>NEG</sub>). The m-th column (1 ≤ m ≤ M<sub>TOT</sub>) of these matrices contain -Log<sub>10</sub>(ERR) values (Expected Recognition Rate, ERR) for best predicted hits of m-th motif in the respective sequence set ([Tsukanov et al., 2022](https://doi.org/10.3389/fpls.2022.938545)). 

GA obtains for these input data the list of motif groups ranked by pAUPRC accuracy, each group includes exactly M motifs, M < M<sub>TOT</sub>.

## GA output data:
- a list of elite (top-scored) M motif groups ranked in the descending order of their pAUPRC recognition accuracy. For the group of one motifs (single motifs), the elite includes all motifs of the input collection, otherwise, the default size of the elite is 100 motifs.
- a list of PR curves for M motif groups, the list is also ranked in the descending order of the pAUPRC values;
- internal structure of motif groups, a list of triangle matrices for M motif groups, computed separately for foreground and background sequence sets, each matrix contains M × (M - 1) / 2 Pearson's correlation coefficients for various pairs of -Log<sub>10</sub>(ERR) vectors representing separate motifs of the same group, the list of matrices is ranked in the descending order of the pAUPRC values;
- external structure of motif groups, two triangle matrices for all elite motifs computed separately for foreground and background sequence sets, each matrix contains M × (M - 1) / 2 Pearson's correlation coefficients for various pairs of -Log<sub>10</sub>(ERR) vectors representing the elite groups. 

# Source code and command line arguments
## TFBS motif recognition for all genes
Preliminary computed data are the results of TFBS motif recognition for promoters of all genes, they represent a table of WG (rows, number all genes in genome) × M<sub>TOT</sub> (columns, number of all motif in the input collection) of -Log<sub>10</sub>(ERR) values, these are best scores of motifs for promoters of all WG genes of genome. This step prepare recognition result for whole genome, this results fits any possible RNA-seq data containing gene lists, only one option may be here 5' and 3' ends of promoters relative to the transcription start sites, recognition for different borders of promoter gene regions can resolve this issues. The next preliminary analysis performs two steps. 

## Preparation of the lists of DEGs and non-DEGs according to rhe results of RNA-seq
[table_rnaseq_adapted.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/table_rnaseq_adapted.cpp) selects the lists of up-/down-regulated DEGs and non-DEGs from the RNA-seq data.
1. input file - table from RNA-seq experiment with a list of gene IDs and log2Fold (Logarithm of the FoldChange value to a base of 2) and padj (adjusted p-value).
2. integer value - column number of gene IDs in the RNA-seq table (argument #1). Currently, for _H. sapiens_ / _M. musculus_, _A. thaliana_ and _D. melanogaster_ Ensembl gene IDs, TAIR AGI codes and FyBase gene IDs are supported, e.g. ENSG00000160072 / ENSMUSG00000033813,
AT1G01200 and FBgn0000008.
3. integer value - column number of log2Fold values in the RNA-seq table (argument #1).
4. integer value - column number of padj values in the RNA-seq table (argument #1).
5. input file - the table for whole genome containing gene IDs (presumed these are all WG protein coding genes of genome). Currently for _H. sapiens_ / _M. musculus_, _A. thaliana_ and _D. melanogaster_ protein coding genes for the recent genome releases hg38, mm38, TAIR10 and dm6 are considered.
6. integer value - the column number of gene IDs in the table for whole genome input file (argument #5).
7. integer value - minimal number of DEGs, selected (either up or down), default values 500 or 250.
8. integer value - minimal number of non-DEGs, selected, default value 1000.
9. double threshold for adjusted p-value (padj), default value for up-/down-regulated DEGs and not DEGs: padj = 0.05.
10. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on up-regulated DEGs, e.g. default: padj < 0.05 & FC >= FCup. 
11. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on down-regulated DEGs, e.g. default: padj < 0.05 & FC <= FCdown.
12. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on non-DEGs, e.g. default: padj > 0.05 &  FCdown <= FC <= FCup.
13. integer value - marks 0 / 1 indicate absence / presence of a swap between up- and down-regultated DEGs in output data

This option forms three files marking for the whole genome list of WG genes up-/down-regulated DEGs and non-DEGs.

[table_rnaseq_knife.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/table_rnaseq_knife.cpp) selects from the RNA-seq data sets of non-overlapped portions of up-/down-regulated DEGs, and for each of these portions the set of the contrast non-DEGs.
1. input file - table from RNA-seq experiment with a list of gene IDs and log2Fold (Logarithm of the FoldChange value to a base of 2) and padj (adjusted p-value).
2. integer value - column number of gene IDs in the RNA-seq table (argument #1). Currently, for _H. sapiens_ / _M. musculus_, _A. thaliana_ and _D. melanogaster_ Ensembl gene IDs, TAIR AGI codes and FyBase gene IDs are supported, e.g. ENSG00000160072 / ENSMUSG00000033813,
AT1G01200 and FBgn0000008.
3. integer value - column number of log2Fold values in the RNA-seq table (argument #1).
4. integer value - column number of padj values in the RNA-seq table (argument #1).
5. input file - the table for whole genome containing gene IDs (presumed these are all WG protein coding genes of genome). Currently for _H. sapiens_ / _M. musculus_, _A. thaliana_ and _D. melanogaster_ protein coding genes for the recent genome releases hg38, mm38, TAIR10 and dm6 are considered.
6. integer value - the column number of gene IDs in the table for whole genome input file (argument #5).
7. integer value - number of DEGs per one portion, selected (either up or down), default values 50 or 25.
8. integer value - number of DEGs portions, default 10.
9. integer value - minimal number of non-DEGs, selected, default value 1000.
10. double threshold for adjusted p-value (padj), default value for up-/down-regulated DEGs and not DEGs: padj = 0.05.
11. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on up-regulated DEGs, e.g. default: padj < 0.05 & FC >= FCup. 
12. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on down-regulated DEGs, e.g. default: padj < 0.05 & FC <= FCdown.
13. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on non-DEGs, e.g. default: padj > 0.05 &  FCdown <= FC <= FCup.
14. integer value - marks 0 / 1 indicate absence / presence of a swap between up- and down-regultated DEGs in output data
15. integer value - marks 1 / -1 indicate taking DEGs with the highest/lowest values of |log2fold|, default values 1. This option is required to study fold-specific distribution of enriched TFBS motifs in DEGs with distinct |log2fold| values.

This step forms three set of files marking gene portions of up-/down-regulated DEGs and sets of WG genes (non-DEGs).

## Extraction TFBS recognition data for DEGs and non-DEGs
[select_lines01.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/select_lines01.cpp) selects the lines of pre-computed TFBS motif recognition data for all up-/down-regulated DEGs and non-DEGs from the RNA-seq data.
1. input file - a table of WG rows (input table).
2. input file - file with WG rows, in each row only one symbol 0 or 1 (input list), so that only N<sub>POS</sub>/Neg rows for up- or down-regulated DEGs / not-DEG contain values of 1.
3. output file - filtered input table containing only rows respecting 1 values in the input list (argument #2).

Now everything is ready for the search of the motif groups.

## Search of TFBS motif groups by Genetic Algorithm
[minimax.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/minimax.cpp) implements the GA search of motif groups.
1. input file - motif recognition table of -Log<sub>10</sub>(ERR) values for up- or down-regulated DEGs (they are required two separate runs), the table has sizes N<sub>POS</sub> (rows, number of up- or down-regulated DEGs) × M<sub>TOT</sub> (columns, number of all motif in the input collection).
2. input file - motif recognition table of -Log<sub>10</sub>(ERR) values for non-DEGs, the table has sizes N<sub>NEG</sub> (rows, number of non-DEGs) × M<sub>TOT</sub> (columns, number of all motif in the input collection).
3. input file - input table representing the collection of all M<sub>TOT</sub> motifs. The table containing four columns: (a) motif IDs from [Hocomoco](https://hocomoco14.autosome.org/) or TF names from [Jaspar2026](https://jaspar.elixir.no/); (b) motif class/family names for these M<sub>TOT</sub> motifs, one unique short description represents each family, depending in the classification of TFs by the DBD structure avalable; (c) indices of these class/family names; (d) notations for the unique family / class names, marks the names encountered in the list for the first time with the number 1, and the rest with the number 0. For RNA-seq data analysis theses input tables are available for several options of the input TFBS motif collection: [_H. sapiens_](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/hocomoco14_tfclass_hs1595), [_M. musculus_](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/hocomoco14_tfclass_mm1245), [_A. thaliana_](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/jaspar2026_tfclass_plants876) and [_D. melanogaster_](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/jaspar2026_tfclass_insects272tbp). _H. sapiens_ / _M. musculus_, _A. thaliana_ and _D. melanogaster_ collections includes 1595/1245 motifs from ([Hocomoco v14](https://hocomoco14.autosome.org/)), 742 motifs for plant TFs from ([Jaspar2026 Plants](https://jaspar.elixir.no/), filtered for -Log<sub>10</sub>(ERR) > 3.4) and 239 motifs for insect TFs, among them 238 motifs are from ([Jaspar2026 Insects](https://jaspar.elixir.no/), filtered for -Log<sub>10</sub>(ERR) > 3.4), and one motif for human TF TBP ([Hocomoco v14](https://hocomoco14.autosome.org/)) is added. 
4. integer value - M<sub>SEL</sub> value, number of motifs in each group (group size), values from 1 to several tens are recommended.
5. integer value - M<sub>SELMAX</sub> value, number of motifs in each group, the maximum value for a series of independent runs, e.g. if these runs imply the group sizes M<sub>SEL</sub> = {1, 2, 3}, then M<sub>SELMAX</sub> = 3.
6. integer value - 0 / 1 marks single / batch run modes, batch mode means concatenatation of data from multiple runs, hence some headers in output files are missed, single mode implies printing of all headers
7. double value - -Log<sub>10</sub>(ERR<sub>MAX</sub>) threshold, it is required to filter out low affinity predictions respecting the left tail of the PR curve, the default value -Log<sub>10</sub>(ERR<sub>MAX</sub>) = 4 respects the probability to find sites (ERR<sub>MAX</sub> = 1E-4, one recognized site per 10 kbp.
8. char path_out - path for output data folder
9. output PR curve file - PR curves for found groups of motifs
10. output file - PR curve for the motif group with the first rank is written in one line: Precision values are marked for Recall values of 0.05, 0.1, 0.15, etc. up to 1. This output file will be concatanated with corresponding files from other runs with other values of the group size M<sub>SEL</sub>. The final output file shows dynamics of the PR curve for the motif group with the first rank as a function of the number of motifs (group size).
11. output external correlation file - correlations of -Log<sub>10</sub>(ERR) scores between found groups of motifs. Two separate correleation matrices are computed: for the foreground and background DNA sequence sets.
12. output internal correlation file - correlations of -Log<sub>10</sub>(ERR) scores between motifs within found groups. For each motif group two separate correleation matrices are computed: for the foreground and background DNA sequence sets.
13. output motif class/family file (all groups) - the list of found groups of motifs, for each groups are marked the recognition accuracy pAUPRC, motif names and class/family names for all participant motifs.
14. output motif class/family file (1st group) - the group of motifs with the first rank, for this group are marked the recognition accuracy pAUPRC, motif names and class/family names for all participant motifs. 
15. output file - distribution of motifs from the group with the first rank by classes/families. This output file will be concatanated with corresponding files from other runs with other values of the group size M<sub>SEL</sub>. The final output file with concatenated results shows dynamics of the family content for the motif group with the first rank as a function of the number of motifs.
16. output file - distribution of motifs from the group with the first rank by present/absent motifs. This output file will be concatanated with corresponding files from other runs with other values of the group size M<sub>SEL</sub>. The final output file with concatenated results shows dynamics of the motif content for the motif group with the first rank as a function of the number of motifs.
17. output log file - GA evolution showing numbers of mutations and recombinations in each iteration of GA.

## Summary table of recognition performance and heatmaps showing emergence ranks and fold enrichments for TFBS motifs for TFs from different families 
[minimax_pars.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/minimax_pars.cpp) creates (a) summary table of recognition accuracy of the groups of motifs with highest recognition performance and (b) heatmaps of emergence ranks, fold enrichments and relative impact to recognition performance for TFBS motifs and TFBS motifs clusterd by TF families. For each run (e.g. up-regulted DEGs vs. not DEGs or down-regulted DEGs vs. not DEGs) totally eight file with text data for heatmaps are prepared; these files differ in the threshold of the fraction of the recognition performance pAUPRC (AUPRC<sub>i</sub> - pAUPRC<sub>1</sub>)/(AUPRC<sub>MAX</sub> - pAUPRC<sub>1</sub>) considered, here pAUPRC<sub>i</sub>, pAUPRC<sub>1</sub> and pAUPRC<sub>MAX</sub> are pAUPRC values for the groups of one motif, several (i) motifs, and for the selected group with the maximal pAUPRC. This approach allows selection of the most important amount of motifs if their number is relatively large.
1. char path_work - path of input and output data. Here input data are output data of GA [minimax.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/minimax.cpp)
2. input PR curve file - output data from the GA search, this file respects argument #9 of [minimax.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/minimax.cpp)
3. input class/family distribution file - output data from the GA search, this file respects argument #14 of [minimax.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/minimax.cpp)
4. input motif distribution file - output data from the GA search, this file respects argument #15 of [minimax.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/minimax.cpp)
5. input file - input table representing the collection of all M<sub>TOT</sub> motifs, this file is the same as described as the argument #3 of [minimax.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/minimax.cpp). Here all names of families/classes, TFs and motifs are required.
6. integer value - number of paralel runs of GA with different group size, default value 50
7. integer value - minimal group size, default value 1
8. integer value - step of group size, default value 1, e.g. if step = 1 then group sizes of 1, 2, ... 50 are tested; if step = 10 then four group sizes of 1, 11, ... 41 are tested.
9. integer value - value 0 or 1 implies that a header of a summary table is printed or not printed.
10. output file - source data for performance chart (default file dynamics), dynamics of the recognition performance for the total range of motif group sizes (default range from 1 to 50 implies 50 preliminary runs of GA for 50 distinct sizes of groups).
11. output file - source data for heatmaps (default files fold_*), fold enrichments for TF families. These folds are defined as follows. If (a) the collection of all motifs has M<sub>TOT</sub> motifs, among them the family FAM has M<sub>FAM</sub> motifs, and (b) in the group of total N motifs only N<sub>FAM</sub> motifs belong to the FAM family, then the fold enrichment is equal to {N<sub>FAM</sub> / N} / M<sub>FAM</sub> / {M<sub>TOT</sub>}. The concatanated heatmap contains values for all group sizes from the minimal group of 1 motif to the maximal group size respecting the highest recognition performance, the same true for the arguments #12, #14, #16 and #18 below.
12. output file - source data for heatmaps (default files fold_best), fold enrichments for families only for one selected group size respecting the highest recognition performance.
13. output file - source data for heatmaps (default files rank_dbd*), ranks for families, if (a) in the group of size {N - 1} the motifs from certain TF family are absent, and (b) at least one motif from this family emerges in the group of size {N} then the rank of this family is N. 
14. output file - source data for heatmaps (default file rank_dbd_best), ranks for families only for one selected group size respecting the highest recognition performance.
15. output file - source data for heatmaps, ranks for motifs, all groups (default files rank_mot*). This is a mask for eight distinct heatmaps.
16. output file - source data for heatmaps, ranks for families only for one selected group size respecting the highest recognition performance (file rank_mot_best).
17. output file - source data for heatmaps, impacts for families. These relative impacts to recofnition performance are defined as follows. The impact for the group of one motif is  pAUPRC<sub>1</sub> /  pAUPRC<sub>1</sub> = 1 (default files ratio_mot*). If (a) in the group of size {N - 1} the motifs from certain TF family are absent, and (b) at least one motif from this family emerges in the group of size {N} then the impact of this family is {pAUPRC<sub>N</sub> - pAUPRC<sub>N-1</sub>} / pAUPRC<sub>1</sub>, here index means the the number of motif in the group. 
18. output file - source data for heatmaps, impacts for families, only for one selected group size respecting the highest recognition performance (default file ratio_mot_best).
19. output file - source data for heatmaps, impacts for motifs, all groups. This is a mask for eight distinct heatmaps.
20. output file - source data for heatmaps impacts for motifs only for one selected group size respecting the highest recognition performance (default file ratio_mot_best).

# Command line examples
## RNA-seq data 
[RNA-seq data](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/E-MTAB-5516-query-results.tsv) [ID E-MTAB-5516](https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5516/resources/ExperimentDownloadSupplier.BulkDifferential/tsv) ([Ragheb et al., 2017](https://doi.org/10.1038/s41598-017-14969-7)) on whole female drosophila wounded and/or infected with _Pseudomonas entomophila_. Data *'wounding' vs 'none' in 'Pseudomonas entomophila' at '3 hour'* were used as a functional example of the tool applicatiion
Two commamd line examples [gmsga_rnaseq_total](https://github.com/parthian-sterlet/minimax/blob/main/run/gmsga_rnaseq_total) and [gmsga_rnaseq_slice](https://github.com/parthian-sterlet/minimax/blob/main/run/gmsga_rnaseq_slice) illustrate two slightly different modes to the tool: 
- The [first](https://github.com/parthian-sterlet/minimax/blob/main/run/gmsga_rnaseq_total) takes a whole set of DEGs, e.g. 250 top-scored by |log2fold| values up- or down-regulated DEGs vs. 1000 non-DEGs.
- The [second](https://github.com/parthian-sterlet/minimax/blob/main/run/gmsga_rnaseq_slice) takes the same whole of DEGs subdivided into certain number of consequent non-overlapping portions, e.g. 5 x 50 top-scored by |log2fold| values up- or down-regulated DEGs vs. 1000 non-DEGs. These separate portion reveal TF regulators potentially involved in the transcription of genes with high or low intensity estimated by the high or low  |log2fold| values.

Gromosega was applied for the motif collection [Jaspar2026 Insecta](https://jaspar.elixir.no/), up-/down-regulated DEGs and non-DEGs were selected by criteria described above for RNA-seq data analysis among [(-1000; +100) regions of promoters of protein coding genes](https://github.com/parthian-sterlet/gromosega/blob/main/genomes/protcod_dm6_m1kb_p100_tss.bed). [Archive](https://github.com/parthian-sterlet/gromosega/blob/main/genomes/protcod_dm6_m1kb_p100_tss_272tbp.tab.tar.gz) represents the best scores for motif recognition of the collection of motifs for {272 insect TFs from Jaspar2026 and human TBP} described above for this whole-genome set of Drosophila promoters. Folder [out](https://github.com/parthian-sterlet/gromosega/tree/main/examples/rnaseq/out) contains the example of main output files, their names respecing four RNA-seq experiments from [E-MTAB-5516](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5516/Results), the suffice *3h* respects the data *'wounding' vs 'none' in 'Pseudomonas entomophila' at '3 hour'*. These files have prefixes *per*, *prb*, *dbd* and *mot*. The prefix *pre* marks output tables for {pAUPRC values, TF names, motif IDs and TF family names}, the prefix *prb* denotes {pAUPRC values and PR curves as Precision values for the standard Recall value set of (5%, 10%, etc. up to 100%)}, and the prefices *dbd* and *mot* mean the distribution for the first-ranked group by TF families and TFBS motif names. Each file contains concatenated output data for all runs for up- and down regulated genes.

# Next applications
[Hocomoco](https://hocomoco.autosome.org/) database (version 14) provided the lists of [1595 human](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/hocomoco14_tfclass_hs1595) / [1245 mouse](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/hocomoco14_tfclass_mm1245) TFBS motifs classified by [TFClass database](https://doi.org/10.1093/nar/gkx987) into families, [Jaspar](https://jaspar.elixir.no/) database (version JASPAR2026) classified TFBS motifs according the [TFClass](https://doi.org/10.1093/nar/gkx987) and [Plant-TFClass](https://doi.org/10.1016/j.tplants.2023.06.023) databases and provided the lists of motifs for plants and insects, the filtration of these lists by [the approach descriped for MCOT package](https://doi.org/10.1093/nar/gkz800) with the threshold of -Log<sub>10</sub>(ERR) = 3.4 brought [876 plant](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/jaspar2024_tfclass_plants742)/[272 insect](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/jaspar2026_tfclass_insects272tbp) TFBS motifs, we added to the list of insects the motif for [human TBP motif](https://hocomoco14.autosome.org/motif/TBP.H14CORE.0.P.B) since we proved earlier ([Zhimulev et al., 2024](https://doi.org/10.3390/ijms25074068)) its enrichement for Drosophila developmental genes (the state 'Ruby' of closed chromatin). 

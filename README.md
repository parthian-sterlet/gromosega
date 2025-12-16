# GroMoSeGA
The GroMoSeGA (GROups of MOtifs SElection by Genetic Algotitm) tool searches for groups of several enriched transcription factor (TF) binding site (TFBS) motifs respecting cooperatively functioning TFs. The motifs are selected from a large collection of annotated TFBS motifs for known TFs derived from the public databases of motifs ([Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077); [Ovek Baydar et al., 2026](https://doi.org/10.1093/nar/gkaf1209)), this collection can be complemented by other motifs obtained by _de novo_ motif search. The enrichment of motif groups is estimated in the comparison of two input sets of positive and negative DNA sequences. The metrics [partial Area under Precision-Recall curve (pAUPRC)](https://github.com/parthian-sterlet/metarea) ([Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90)) is used to transit from single motifs to motif groups. The approach of Genetic Algorithm (GA) is used to avoid the exhaustive search in the combinatorially large space of possible motif content of groups. The tool is designed for the analysis of RNA-seq and ChIP-seq data, and may be applied for any other DNA sequence sets presuming action of different transcription regulatory machines. Hence, the positive/negative sequences are the promoters of differentially expressed genes (DEGs) / not DEGs for RNA-seq data, and ChIP-seq peaks and randomly selected genomic loci for ChIP-seq data. 

# Description
 GroMoSeGA took the basic computation principle from [MetArea tool](https://github.com/parthian-sterlet/metarea). For a motif representing TFBS pattern, the recognition accuracy is calculated as the recogntion performance pAUPRC ([Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90); [Davis, Goadrich, 2006](https://doi.org/10.1145/1143844.1143874)). Several motifs compiles a group of motifs. For this group for given DNA sequence a motif of this group with the best hit among all other motifs of this group defines the hit of the whole group, i.e. different motifs may represent group in various sequences. Thus, for the group of motifs the one PR curve and its pAUPRC value is computed. Hence, we can search the groups that better in terms of the pAUPRC measure distinguishes the positive and negative DNA sequence sets. Testing and optimization of the motif content for the population of many groups by the approach of genetic algorithm (GA) allows to detect and rank groups of motifs that reflect the difference between two sequence sets. The motifs of the same group can represent structurally different binding sites for the same TF and binding sites of different TFs acting together in single multiprotein complexes (Levitsky et al., [2014](https://doi.org/10.1186/1471-2164-15-80), [2024](https://doi.org/10.18699/vjgb-24-90)). Since the default option of TFs action is mutual cooperative interactions ([Morgunova & Taipale, 2017](https://doi.org/10.1016/j.sbi.2017.03.006)), the single protein complex of multiple TFs is reflected in various portions of peaks as different enriched motifs of various TFs, including TFs participating in transcription regulation of DEGs from the RNA-seq data, and the target and partner TFs related to TFBS motif enrichment for ChIP-seq data. Hence, the functional relationship of several distinct motifs of the same group is recognized through increased recognition performance. For the fixed number of motifs, the groups with the greater pAUPRC values are the better representatives of the functional TFs compared to other groups with the smaller pAUPRC values. Thus, GroMoSeGA predicts TFBS motifs respecting TFs participating in the same multiprotein complexes functioning in gene transcription regulation.

# Requirements
The source code is written in C++ language. To compile execubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Repository structure
Folder [**src**](https://github.com/parthian-sterlet/gromosega/tree/main/src) contains the [major](https://github.com/parthian-sterlet/gromosega/blob/main/src/minimax.cpp) and some supporting C++ source code files.  

Folder [**run**](https://github.com/parthian-sterlet/antinoise/tree/main/run) contains two command lines, implementing the tool pipeline for example RNA-seq and ChIP-seq data.

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
tar -xvzf protcod_dm6_m1kb_p100_tss_247tbp.tab.tar.gz
cd ..

```
* In Windows system:

separate compilation of all source files in VC++

# Algorithm
GroMoSeGA algorithm considers a pair of positive/negative sequence sets derived either from RNA-seq or ChIP-seq data. peaks. For RNA-seq data the positive and negative sequences are promoters of DEGs and not-DEGs, that are defined by the default criteria: 
- {adjusted p-value < 0.05 & log2(FoldChange) > 1} / {adjusted p-value < 0.05 & log2(FoldChange) < -1}  for up-/down-regulated DEGs, and 
- {adjusted p-value > 0.05 & 0.8 < FoldChange < 1.25} for not DEGs.

For ChIP-seq data the positive / negative sequences are ChIP-seq peaks / randomly selected genomic loci, correspondingly. In the case if RNA-seq data, the , respectively. For ChIP-seq/ATAC-seq data the negative set contains [randomly selected genomic loci, adopted by G/C-content selected by the AntiNoise tool](https://github.com/parthian-sterlet/antinoise/). 

## GA input data:
- a pair of positive / negative sets, containng N<sub>POS</sub> / N<sub>NEG</sub> DNA sequences;
- a collection of M<sub>TOT</sub> TFBS motifs derived from public databases like [Hocomoco v 14](https://hocomoco14.autosome.org/) ([Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077)) or [Jaspar2026](https://jaspar.elixir.no/) ([Ovek Baydar et al., 2026](https://doi.org/10.1093/nar/gkaf1209)), this collection can be complemented by motifs derived by _de novo_ motif search for an appropriate dataset, for each motif it is required its name (e.g. TF name and motif designation) and short description of its TF class/family according to the hierarchical classification of mammalian ([Wingender et al., 2018](https://doi.org/10.1093/nar/gkx987)) and plant ([Blanc-Mathieu et al., 2024](https://doi.org/10.1016/j.tplants.2023.06.023)) TFs by DNA-binding domain (DBD) structure;
- two matrices of the sizes M<sub>TOT</sub> x N<sub>POS</sub> and M<sub>TOT</sub> x N<sub>NEG</sub>. The rows of these matrices represent sequences (1 ≤ n ≤ N<sub>POS</sub>, 1 ≤ n ≤ N<sub>NEG</sub>). The m-th column (1 ≤ m ≤ M<sub>TOT</sub>) of these matrices contain -Log<sub>10</sub>(ERR) values (Expected Recognition Rate, ERR) for best predicted hits of m-th motif in the respective sequence set ([Tsukanov et al., 2022](https://doi.org/10.3389/fpls.2022.938545)). 

GA obtains for these input data the list of motif groups ranked by pAUPRC accuracy, each group includes exactly M motifs, M < M<sub>TOT</sub>.

## GA output data:
- a list of elite (top-scored) ME motif groups ranked in the descending order of their pAUPRC recognition accuracy. For the group of one motifs (single motifs), the elite includes all motifs of the input collection, otherwise, the default size of the elite is 100 motifs.
- a list of PR curves for ME motif groups, the list is also ranked in the descending order of the pAUPRC values;
- internal structure of motif groups, a list of triangle matrices for ME motif groups, computed separately for positive and negative sequence sets, each matrix contains M × (M - 1) / 2 Pearson's correlation coefficients for various pairs of -Log<sub>10</sub>(ERR) vectors representing separate motifs of the same group, the list of matrices is ranked in the descending order of the pAUPRC values;
- external structure of motif groups, two triangle matrices for all elite motifs computed separately for positive and negative sequence sets, each matrix contains ME × (ME - 1) / 2 Pearson's correlation coefficients for various pairs of -Log<sub>10</sub>(ERR) vectors representing the elite groups. 

# Source code and command line arguments

The analysis of ChIP-seq data is more simple, the positive and negative DNA sequence sets are almost ready. Below a pipeline for the more complicate case of RNA-seq data analysis is described.

## Zero step
Preliminary computed data are the results of TFBS motif recognition for promoters of all genes, they represent a table of WG (rows, number all genes in genome) × M<sub>TOT</sub> (columns, number of all motif in the input collection) of -Log<sub>10</sub>(ERR) values, these are best scores of motifs for promoters of all WG genes of genome. This step prepare recognition result for whole genome, this results fits any possible RNA-seq data containing gene lists, only one option may be here 5' and 3' ends of promoters relative to the transcription start sites, recognition for different borders of promoter gene regions can resolve this issues. The next preliminary analysis performs two steps. 

## First step
[table_rnaseq_filter.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/table_rnaseq_filter.cpp) selects the lists of up-/down-regulated DEGs and not-DEGs from the RNA-seq data.
1. input file - table from RNA-seq experiment with a list of gene IDs and log2Fold (Logarithm of the FoldChange value to a base of 2) and padj (adjusted p-value).
2. integer value - column number of gene IDs in the RNA-seq table (argument #1). Currently, for _H. sapiens_ / _M. musculus_, _A. thaliana_ and _D. melanogaster_ Ensembl gene IDs, TAIR AGI codes and FyBase gene IDs are supported, e.g. ENSG00000160072 / ENSMUSG00000033813,
AT1G01200 and FBgn0000008.
3. integer value - column number of log2Fold values in the RNA-seq table (argument #1).
4. integer value - column number of padj values in the RNA-seq table (argument #1).
5. input file - the table for whole genome containing gene IDs (presumed these are all WG protein coding genes of genome). Currently for _H. sapiens_ / _M. musculus_, _A. thaliana_ and _D. melanogaster_ protein coding genes for the recent genome releases hg38, mm39, TAIR10 and dm6 are considered.
6. integer value - the column number of gene IDs in the table for whole genome (argument #5).
7. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on up-regulated DEGs, {adjusted p-value < 0.05 & log2(FoldChange) > 1.
8. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on down-regulated DEGs, {adjusted p-value < 0.05 & log2(FoldChange) < -1.
9. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on not-DEGs, {adjusted p-value > 0.05 &  0.8 < FoldChange) < 1.25.

This first step forms three files marking for the whole genome list of WG genes up-/down-regulated DEGs and not-DEGs.

## Second step 
[select_lines01.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/select_lines01.cpp) selects the lines of pre-computed TFBS motif recognition data for all up-/down-regulated DEGs and not-DEGs from the RNA-seq data.
1. input file - a table of WG rows (input table).
2. input file - file with WG rows, in each row only one symbol 0 or 1 (input list), so that only N<sub>POS</sub>/Neg rows for up- or down-regulated DEGs / not-DEG contain values of 1.
3. output file - filtered input table containing only rows respecting 1 values in the input list (argument #2).

Now everything is ready for the search of the motif groups.

## Main analysis step
[minimax.cpp](https://github.com/parthian-sterlet/gromosega/blob/main/cpp/minimax.cpp) implements the GA search of motif groups.
1. input file - motif recognition table of -Log<sub>10</sub>(ERR) values for up- or down-regulated DEGs (they are required two separate runs), the table has sizes N<sub>POS</sub> (rows, number of up- or down-regulated DEGs) × M<sub>TOT</sub> (columns, number of all motif in the input collection).
2. input file - motif recognition table of -Log<sub>10</sub>(ERR) values for not-DEGs, the table has sizes N<sub>NEG</sub> (rows, number of not-DEGs) × M<sub>TOT</sub> (columns, number of all motif in the input collection).
3. input file - input table representing the collection of all M<sub>TOT</sub> motifs. The table containing four columns: (a) motif IDs from [Hocomoco](https://hocomoco14.autosome.org/) or TF names from [Jaspar2026](https://jaspar.elixir.no/); (b) motif class/family names for these M<sub>TOT</sub> motifs, one unique short description represents each family/class, depending in the classification of TFs by the DBD structure avalable; (c) indices of these class/family names; (d) notations for the unique family / class names, marks the names encountered in the list for the first time with the number 1, and the rest with the number 0. For RNA-seq data analysis theses input tables are available for several options of the input TFBS motif collection: [_H. sapiens_](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/hocomoco14_tfclass_hs1595), [_M. musculus_](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/hocomoco14_tfclass_mm1245), [_A. thaliana_](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/jaspar2026_tfclass_plants859) and [_D. melanogaster_](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/jaspar2026_tfclass_insects247tbp). _H. sapiens_ / _M. musculus_, _A. thaliana_ and _D. melanogaster_ collections includes 1595/1245 motifs from ([Hocomoco v14](https://hocomoco14.autosome.org/)), 742 motifs for plant TFs from ([Jaspar2026 Plants](https://jaspar.elixir.no/), filtered for -Log<sub>10</sub>(ERR) > 3.6) and 239 motifs for insect TFs, among them 238 motifs are from ([Jaspar2026 Insects](https://jaspar.elixir.no/), filtered for -Log<sub>10</sub>(ERR) > 3.6), and one motif for human TF TBP ([Hocomoco v14](https://hocomoco14.autosome.org/))} is added. For ChIP-seq data analysis this input table should be prepared separately for each set of input TFBS motifs, [XLSX file example](https://github.com/parthian-sterlet/gromosega/blob/main/examples/chipseq/Gromosega_table_for_TFClass_annotation.xlsx) provides the explanation of this input table for given ChIP-seq data [GSM2845664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2845664) and respective [list](https://github.com/parthian-sterlet/gromosega/tree/main/examples/chipseq) of motifs from _de novo_ motif search.
4. integer value - M<sub>SEL</sub> value, number of motifs in each group (group size), values from 1 to several tens are recommended.
5. integer value - M<sub>SELMAX</sub> value, number of motifs in each group, the maximum value for a series of independent runs, e.g. if these runs imply the group sizes M<sub>SEL</sub> = {1, 2, 3}, then M<sub>SELMAX</sub> = 3.
6. double value - -Log<sub>10</sub>(ERR<sub>MAX</sub>) threshold, ot is required to filter out a left tail of the PR curve where the potential sites of the lowest affinity are expected, typically this value should be in the range from the more stringent maximal -Log<sub>10</sub>(ERR<sub>MAX</sub>) = 3.3 (ERR<sub>MAX</sub> = 5E-4, ~ 1 site is recognized per 2kb) to the more mild minimal -Log<sub>10</sub>(ERR<sub>MAX</sub>) = 2.69 (ERR<sub>MAX</sub> = 2E-3, 1 site per 500 bp).
7. output PR curve file - PR curves for found groups of motifs
8. output file - PR curve for the motif group with the first rank is written in one line: Precision values are marked for Recall values of 0.05, 0.1, 0.15, etc. up to 1. This output file will be concatanated with corresponding files from other runs with other values of the group size M<sub>SEL</sub>. The final output file shows dynamics of the PR curve for the motif group with the first rank as a function of the number of motifs (group size).
9. output external correlation file - correlations of -Log<sub>10</sub>(ERR) scores between found groups of motifs. Two separate correleation matrices are computed: for the positive and negative DNA sequence sets.
10. output internal correlation file - correlations of -Log<sub>10</sub>(ERR) scores between motifs within found groups. For each motif group two separate correleation matrices are computed: for the positive and negative DNA sequence sets.
11. output motif class/family file (all groups) - the list of found groups of motifs, for each groups are marked the recognition accuracy pAUPRC, motif names and class/family names for all participant motifs.
12. output motif class/family file (1st group) - the group of motifs with the first rank, for this group are marked the recognition accuracy pAUPRC, motif names and class/family names for all participant motifs.
13. output file - distribution of motifs from the group with the first rank by classes/families. This output file will be concatanated with corresponding files from other runs with other values of the group size M<sub>SEL</sub>. The final output file shows dynamics of the family/class content for the motif group with the first rank as a function of the number of motifs (group size).
14. output log file - GA evolution showing numbers of mutations and recombinations in each iteration of GA.
    
# Command line examples
## ChIP-seq data 
[commamd_line example](https://github.com/parthian-sterlet/gromosega/blob/main/run/chipseq_com_line) for ChIP-seq data ([Czimmerer et al., 2018](https://doi.org/10.1016/j.immuni.2017.12.010)) [GSM2845664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2845664) for _M. musculus_ target TF STAT6 in bone marrow-derived macrophage (BMDM) cells differentiated in the presence of MCSF (Macrophage Colony-Stimulating Factor) and treated (1h) with interleukin-4 (IL-4) cytokine. _De novo_ motif search tool [STREME](https://doi.org/10.1093/bioinformatics/btab203) applied for the foreground set of [ChIP-seq peaks](https://github.com/parthian-sterlet/gromosega/blob/main/examples/chipseq/GSM2845664_BMDM_1hIL4_STAT6_mm10.fa) and background set of [random genomic loci derived by AntiNoise tool](https://github.com/parthian-sterlet/gromosega/blob/main/examples/chipseq/GSM2845664_BMDM_1hIL4_STAT6_mm10_gb.fa) revealed 14 enriched motifs respecting [proven TFBS motifs for known TFs](https://github.com/parthian-sterlet/gromosega/blob/main/examples/chipseq/mnames). GroMoseGA was applied for these motifs and the same foreground/background sequence sets.
## RNA-seq data 
[commamd_line example](https://github.com/parthian-sterlet/minimax/blob/main/run/rnaseq_com_line) illustrates RNA-seq data ([Jafari et al., 2021](https://doi.org/10.1371/journal.pbio.3001101)) [E-MTAB-9805](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-9805/Results) on stress and neuronal activity on _D. melanogaster_ antennae, mutant genotype Orco-/- vs. wild type, 1 day old (newly hatched flies). Mutants with Orco-/- genotype have a severe loss of olfactory sensory neurons since [Orco protein](https://flybase.org/reports/FBgn0037324) is essential for the function and development of odorant receptor complexes. Hence, mutants are unable to respond to most odors. GroMoseGA was applied for the motif collection [Jaspar2026 Insects](https://jaspar.elixir.no/), down-regulated DEGs and not-DEGs were selected by criteria described above for RNA-seq data analysis among [(-1000; +100) regions of promoters of protein coding genes](https://github.com/parthian-sterlet/gromosega/blob/main/genomes/protcod_dm6_m1kb_p100_tss.bed). [Archive](https://github.com/parthian-sterlet/gromosega/blob/main/genomes/protcod_dm6_m1kb_p100_tss_247tbp.tab.tar.gz) represents the best scores for motif recognition of the collection of motifs for {247 insect TFs from Jaspar2026 and human TBP} described above for this whole-genome set of Drosophila promoters.

# Next applications
[Hocomoco](https://hocomoco.autosome.org/) database (version 14) provided the lists of [1595 human](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/hocomoco14_tfclass_hs1595) / [1245 mouse](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/hocomoco14_tfclass_mm1245) TFBS motifs classified by [TFClass database](https://doi.org/10.1093/nar/gkx987) into families, [Jaspar](https://jaspar.elixir.no/) database (version JASPAR2026) classified TFBS motifs according the [TFClass](https://doi.org/10.1093/nar/gkx987) and [Plant-TFClass](https://doi.org/10.1016/j.tplants.2023.06.023) databases and provided the lists of motifs for plants and insects, the filtration of these lists by [the approach descriped for MCOT package](https://doi.org/10.1093/nar/gkz800) with the threshold of -Log<sub>10</sub>(ERR) = 3.6 brought [742 plant](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/jaspar2024_tfclass_plants742)/[239 insect](https://github.com/parthian-sterlet/gromosega/blob/main/examples/rnaseq/jaspar2026_tfclass_insects247tbp) TFBS motifs, we added to the list of insects the motif for [human TFBS motif](https://hocomoco14.autosome.org/motif/TBP.H14CORE.0.P.B) since we proved earlier ([Zhimulev et al., 2024](https://doi.org/10.3390/ijms25074068)) it enrichement for a substantial part of all Drosophila genes (developmental genes). 

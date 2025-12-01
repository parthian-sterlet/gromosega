# GroMoSeGA
The GroMoSeGA (GROups of MOtifs SElection by Genetic Algotitm) tool searches for groups of several enriched transcription factor (TF) binding site (TFBS) motifs respecting cooperatively functioning TFs. The motifs are selected from the large library of TFBS motifs for known TFs derived from public databases of motifs, this library can be complemented by other motifs obtained by _de novo_ motif search. The motif enrichment is estimated in the comparison of the given positive and negative DNA sequence set. The tool is designed for the analysis of RNA-seq and ChIP-seq data. Hence, the positive/negative sequences are the promoters of differentially expressed genes (DEGs) / not DEGs for RNA-seq data, and ChIP-seq peaks and randomly selected genomic loci for ChIP-seq data. 

# Description
Basic computation principle of the GroMoSeGA tool comes from the [MetArea tool](https://github.com/parthian-sterlet/metarea). For a motif representing TFBSs, the recognition accuracy is calculated as the partial area under the PR (Precision-Recall) curve (pAUPRC, [Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90); [Davis, Goadrich, 2006](https://doi.org/10.1145/1143844.1143874)). Several motifs compiles a group of motifs. For this group for any DNA sequence certain motif of this group with the best hit among other motifs of this group defines the hit of the whole group, i.e. different motifs may represent group in various sequences. Thus, for the group of motifs the one PR curve and its pAUPRC value is computed. Hence, we can search the groups that better in terms of the pAUPRC measure distinguishes the positive and negative DNA sequence sets. Testing and optimizarion of the motif content for the population of many groups by the approach of genetic algorithm (GA) allows to detect and rank groups of motifs that reflect the diference between two sequence sets. The motifs of the same group can represent structurally different binding sites for the same TF and binding sites of different TFs acting together in single multiprotein complexes (Levitsky et al., [2014](https://doi.org/10.1186/1471-2164-15-80), [2024](https://doi.org/10.18699/vjgb-24-90)). Since the default option of TFs action is mutual cooperative interactions ([Morgunova & Taipale, 2017](https://doi.org/10.1016/j.sbi.2017.03.006)), the single protein complex of multiple TFs is reflected in various portions of peaks as different enriched motifs of various TFs, including TFs participating in transcription regulation of DEGs from the RNA-seq data, and the target and partner TFs related to TFBS motif enrichment for ChIP-seq data. Hence, the functional relationship of several distinct motifs of the same group is recognized through increased recognition performance. For the fixed number of motifs, the groups with the greater pAUPRC values are the better representatives of the functional TFs compared to other groups with the smaller pAUPRC values. Thus, GroMoSeGA predicts TFBS motifs respecting TFs paricipating in the same multiprotein complexes functioning in gene transcription regulation.

# Requirements
The source code is written in C++ language. To compile exetubables from the source code you need:

* In Linux system, C++ compiler, e.g. [GCC](https://gcc.gnu.org/) compiler 
* In Windows system any VC++ package, e.g. [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/vs/community/)

# Repository structure
Folder [**src**](https://github.com/parthian-sterlet/GroMoSeGA/tree/main/src) contains the [major](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/src/GroMoSeGA.cpp) and some supporting C++ source code files.  

Folder [**run**](https://github.com/parthian-sterlet/antinoise/tree/main/run) contains perl scripts command line examples, implementing the tool pipeline for RNA-seq and ChIP-seq data.

Folder [**examples**](https://github.com/parthian-sterlet/GroMoSeGA/tree/main/bin/examples) contains files required as the functional examples of the tool.

Folder [**genomes**](https://github.com/parthian-sterlet/antinoise/tree/main/genomes) contains whole genome anntations required for mentioned above functional examples. 

# How to compile
* In Linux system: 
```
git clone https://github.com/parthian-sterlet/GroMoSeGA
cd GroMoSeGA/run
chmod a+x build.sh
./build.sh
cd ../genomes
tar -xvzf m1kb_p1_dm6_err238tbp.tab.tar.gz
cd ..

```
* In Windows system:

separate compilation of all source files in VC++

# Algorithm
GroMoSeGA algorithm considers a pair of positive/negative sequence sets derived either from from RNA-seq or ChIP-seq data. peaks. For RNA-seq data the positive and negative sequences are promoters of DEGs and not-DEGs, that are defined by the default criteria {adjusted p-value < 0.05 & log2(FoldChange) > 1 / log2(FoldChange) < -1  for up-/down-regulated DEGs} and {adjusted p-value > 0.05 & 0.8 < FoldChange) < 1.25}, correspondingly. For ChIP-seq data the positive/negative sequences are ChIP-seq peaksor randomly selected genomic loci, correspondingly. In the case if RNA-seq data, the , respectively. For ChIP-seq/ATAC-seq data the negative set contains [randomly selected genomic loci, adopted by G/C-content selected by the AntiNoise tool](https://github.com/parthian-sterlet/antinoise/). 

## GA input data:
- a pair of positive/negative sets, containg Npos/Nneg sequences;
- a library of Mtot TFBS motifs derived from public databases like [Hocomoco](https://hocomoco14.autosome.org/) ([Vorontsov et al., 2024](https://doi.org/10.1093/nar/gkad1077)) or [Jaspar](https://jaspar.elixir.no/) ([Rauluseviciute et al., 2024](https://doi.org/10.1093/nar/gkad1059)), this library can be complemented by motifs derived by _de novo_ motif search for an apropriate dataset, for each motif it is required its name (e.g. TF name and motif designation) and short description of its TF class/family according to the hierarchical classification of mammalian ([Wingender et al., 2018](https://doi.org/10.1093/nar/gkx987)) and plant ([Blanc-Mathieu et al., 2024](https://doi.org/10.1016/j.tplants.2023.06.023)) TFs by DNA-binding domain (DBD) structure;
- two matrices of the sizes Mtot x Npos and Mtot x Nneg. The rows of these matrices represent sequences (1 ≤ n ≤ Npos, 1 ≤ n ≤ Nneg). The m-th column (1 ≤ m ≤ Mtot) of these matrices contain -Log<sub>10</sub>(ERR) values (Expected Recognition Rate, ERR) for best predicted hits of m-th motif in the respective sequence set ([Tsukanov et al., 2022](https://doi.org/10.3389/fpls.2022.938545)). 

GA obtains for these input data the list of motif groups ranked by pAUPRC accuracy, each group includes exactly M motifs, M < Mtot.

## GA output data:
- a list of elite (top-scored) ME motif groups ranked in the descending order of their pAUPRC recognition accuracy. For the group of one motifs (single motifs), the elite includes all motifs of the input library, otherwise, the default size of the elite is 100 motifs.
- a list of PR curves for ME motif groups, the list is also ranked in the descending order of the pAUPRC values;
- internal structure of motif groups, a list of triangle matrices for ME motif groups, computed separately for positive and negative sequence sets, each matrix contains M × (M - 1) / 2 Pearson's correlation coefficients for various pairs of -Log<sub>10</sub>(ERR) vectors representing separate motifs of the same group, the list of matrices is ranked in the descending order of the pAUPRC values;
- external structure of motif groups, two triangle matrices for all elite motifs computed separately for positive and negative sequence sets, each matrix contains ME × (ME - 1) / 2 Pearson's correlation coefficients for various pairs of -Log<sub>10</sub>(ERR) vectors representing the elite groups. 

# Source code and command line arguments

The analysis of ChIP-seq data is more simple, the positive and negative DNA sequence sets are almost ready. Below a pipeline for the more complicate case of RNA-seq data analysis is described.

## Zero step
Preliminary computed data are the results of TFBS motif recognition for promoters of all genes, they represent a table of WG (rows, number all genes in genome) × Mtot (columns, number of all motif in the input library) of -Log<sub>10</sub>(ERR) values, these are best scores of motifs for promoters of all WG genes of genome. This step prepare recognition result for whole genome, this results fits any possible RNA-seq data containing gene lists, only one option may be here 5' and 3' ends of promoters relative to the transcription start sites, recognition for diiferent edges of promoters can resolve this issues. The next preliminary analysis performs two steps. 

## First step
[table_rnaseq_filter.cpp](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/cpp/table_rnaseq_filter.cpp) selects the lists of up-/down-regulated DEGs and not-DEGs from the RNA-seq data.
1. input file - table from RNA-seq experiment with a list of gene IDs and log2Fold (Logarithm of the FoldChange value to a base of 2) and padj (adjusted p-value).
2. integer value - column number of gene IDs in the RNA-seq table (argument #1). Currently, for _H. sapiens_ / _M.musculus_, _A. thaliana_ and _D. melanogaster_ Ensembl gen ID, TAIR AGI codes and FyBase gene ID are sipported, e.g. ENSG00000160072, AT1G01200 and FBgn0000008
3. integer value - column number of log2Fold values in the RNA-seq table (argument #1).
4. integer value - column number of padj values in the RNA-seq table (argument #1).
5. input file - the table for whole genome containing gene IDs (presumed these are all WG protein coding genes of genome). Currently for  _H. sapiens_ / _M.musculus_, _A. thaliana_ and _D. melanogaster_ 19795/19991 (hg38/mm10), 27202 (TAIR10) and 13773 (dm6) genes are considered.
6. integer value - the column number of gene IDs in the table for whole genome (argument #5).
7. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on up-regulated DEGs, {adjusted p-value < 0.05 & log2(FoldChange) > 1.
8. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on down-regulated DEGs, {adjusted p-value < 0.05 & log2(FoldChange) < -1.
9. output file -list of all WG integer values (0 or 1) marking gene satisfying the default criterion on not-DEGs, {adjusted p-value > 0.05 &  0.8 < FoldChange) < 1.25.

This first step forms three files marking for the whole genome list of WG genes up-/down-regulated DEGs and not-DEGs.

## Second step 
[select_lines01.cpp](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/cpp/select_lines01.cpp) selects the lines of pre-computed TFBS motif recognition data for all up-/down-regulated DEGs and not-DEGs from the RNA-seq data.
1. input file - a table of WG rows (input table).
2. input file - file with WG rows, in each row only one symbol 0 or 1 (input list), so that only Npos/Neg rows for up- or down-regulated DEGs / not-DEG contain values of 1.
3. output file - filtered input table containing only rows respecting 1 values in the input list (argument #2).

Now everithing is ready for the search of motif groups.

## Main analysis step
[minimax.cpp](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/cpp/minimax.cpp) implements the GA search of motif groups.
1. input file - motif recognition table of -Log<sub>10</sub>(ERR) values for up- or down-regulated DEGs (they are required two separate runs), the table has sizes Npos (rows, number of up- or down-regulated DEGs) × Mtot (columns, number of all motif in the input library).
2. input file - motif recognition table of -Log<sub>10</sub>(ERR) values for not-DEGs, the table has sizes Nneg (rows, number of not-DEGs) × Mtot (columns, number of all motif in the input library).
3. input file - list of all Mtot motif names for the input library (for [Jaspar](https://jaspar.elixir.no/) these are TF names, for [Hocomoco(https://hocomoco.autosome.org/) - motif IDs), this list includes Mtot motifs, Mtot is total number of motifs in the input library. Currently, for  _H. sapiens_ / _M.musculus_, _A. thaliana_ and _D. melanogaster_ these numbers are 1595/1245 ([Hocomoco v14](https://hocomoco14.autosome.org/)), 742 ([Jaspar Plants](https://jaspar.elixir.no/), filtered for -Log<sub>10</sub>(ERR) > 3.6) and 239 {238 ([Jaspar Insects](https://jaspar.elixir.no/), filtered for -Log<sub>10</sub>(ERR) > 3.6) + human motif for TBP TF ([Hocomoco v14](https://hocomoco14.autosome.org/))}. See examples for [_H. sapiens_](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/rnaseq/namemot_hcmc1595), [_A. thaliana_](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/rnaseq/name_jaspar742) and [_D. melanogaster_](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/rnaseq/name_jaspar238_tbp), this list includes Mtot motifs.
4. input file - list of all Mtot motif class/family names for the input library. One unique short description is required. See examples for [_H. sapiens_](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/rnaseq/nameclass_hcmc1595), [_A. thaliana_](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/rnaseq/name_class742) and [_D. melanogaster_](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/rnaseq/name_class238_tbp), this list includes Mtot motifs.
6. integer value - M value, number of motifs in each group, values from 1 to several tens are recommended.
7. double value - -Log<sub>10</sub>(ERR<sub>MAX</sub>) threshold, ot is required to filter out a left tail of the PR curve where the potential sites of the lowest affinity are expected, typically this value should be in the range from the more stringent maximal -Log<sub>10</sub>(ERR<sub>MAX</sub>) = 3.3 (ERR<sub>MAX</sub> = 5E-4, ~ 1 site is recognized per 2kb) to the more mild minimal -Log<sub>10</sub>(ERR<sub>MAX</sub>) = 2.69 (ERR<sub>MAX</sub> = 2E-3, 1 site per 500 bp).
8. mask for output PR curve file - PR curves for found groups of motifs, output file name is "mask"_M.
9. mask for output internal correlation file - correlations of -Log<sub>10</sub>(ERR) scores between found groups of motifs, output file name is "mask"_M.
10. mask for output external correlation file - correlations of -Log<sub>10</sub>(ERR) scores between motifs for found groups, output file name is "mask"_M.
11. mask for main output file - the list of found groups of motifs, for each groups are marked the recognition accuracy pAUPRC (fitness function from GA), motif names and class/family names for all M motifs, output file name is "mask"_M.
12. mask for output file - log of GA evolution showing numbers of mutations and recombinations in iterations of GA, output file name is "mask"_M.
    
# Command line examples
- ChIP-seq data: [comnamd_line example](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/run/chipseq_com_line) for ChIP-seq data [GSM2845664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2845664) for _M.musculus_ target TF STAT6 in bone marrow-derived macrophage (BMDM) cells differentiated in the presence of MCSF (Macrophage Colony-Stimulating Factor) and treated (1h) with interleukin-4 (IL-4) cytokine. _De novo_ motif search tool [STREME](https://doi.org/10.1093/bioinformatics/btab203) applied for the foreground set of [ChIP-seq peaks](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/chipseq/GSM2845664_BMDM_1hIL4_STAT6_mm10.fa) and background set of [random genome loci derived by AntiNoise tool](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/chipseq/GSM2845664_BMDM_1hIL4_STAT6_mm10_gb.fa) revealed 14 enriched motifs respecting [proven TFBS motifs for known TFs](https://github.com/parthian-sterlet/GroMoSeGA/blob/main/examples/chipseq/mnames). GroMoseGA was applied for these motifs and the same forground/background sequence sets.
-  RNA-seq data: [comnamd_line example](https://github.com/parthian-sterlet/minimax/blob/main/run/rnaseq_com_line) illustrates RNA-seq data [E-MTAB-9805](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-9805/Results) on stress and neuronal activity in _D. melanigaster_ antenna (Orco-/- vs. wild type genotype, 1 day). GroMoseGA was applied for applied for the motif library [Jaspar Insects](https://jaspar.elixir.no/), down-regulated DEGs and not-DEGs were selected by criteria described above for RNA-seq data analysis.

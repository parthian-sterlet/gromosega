#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <cmath>

#define MMAX 1700  //max total number of motifs
#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
// parsing of gromosega
struct grpa // for all runs
{
	double rat;//(allg[i].auprc - auprc_min) / auprc_max - auprc_min)
	double auprc;
	double ratio; ////0 -> 1, i -> (allg[i].auprc - allg[i - 1].auprc) / allg[0].auprc)
	int gsi;
	int irez;
};
struct grps //for selected 8 points
{
	double rat;
	double auprc;
	int gsi;
	int irez;	
};
int StrNStr(char* str, char c, int n)
{
	if (n == 0)return -1;
	int i, len = (int)strlen(str);
	int k = 1;
	for (i = 0; i < len; i++)
	{
		if (str[i] == c)
		{
			if (k == n)return i;
			k++;
		}
	}
	return -10000;
}
void DelHole(char* str)
{
	char* hole;
	hole = strstr(str, "\n");
	if (hole != NULL) *hole = 0;
}
double UnderStolDouble(char* str, int nstol, char* ret, size_t size, char razd)
{
	memset(ret, 0, size);
	if (nstol == 0)return atof(str);
	int p1 = StrNStr(str, razd, nstol);
	int p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -10000)
	{
		p2 = (int)strlen(str);
	}
	if (p1 == -10000 || p2 == -10000) return -10000;
	int len = p2 - p1 - 1;
	strncpy(ret, &str[p1 + 1], len);
	ret[len] = '\0';
	int cd = (int)ret[0];
	if (isdigit(cd) || cd == 45)return atof(ret);//0123456789 or -
	else return -10000;
}
int UnderStolInt(char* str, int nstol, char* ret, size_t size, char razd)
{
	memset(ret, 0, size);
	if (nstol == 0)return atoi(str);
	int p1 = StrNStr(str, razd, nstol);
	int p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -10000)
	{
		p2 = (int)strlen(str);
	}
	if (p1 == -10000 || p2 == -10000) return -10000;
	int len = p2 - p1 - 1;
	strncpy(ret, &str[p1 + 1], len);
	ret[len] = '\0';
	int cd = (int)ret[0];
	if (isdigit(cd) || cd == 45)return atoi(ret);//0123456789 or -
	else return -10000;
}
char* UnderStolStr(char* str, char* ret, size_t size, int nstol, char razd)
{
	memset(ret, 0, size);
	int p1, p2;
	if (nstol == 0)p1 = -1;
	else p1 = StrNStr(str, razd, nstol);
	p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -10000)
	{
		p2 = (int)strlen(str);
	}
	if (p1 == -10000 || p2 == -10000) return NULL;
	int len = p2 - p1 - 1;
	strncpy(ret, &str[p1 + 1], len);
	ret[len] = '\0';
	return ret;
}
// delete symbol 'c' from input string
void DelChar(char* str, char c)
{
	int i, lens, size;

	size = 0;
	lens = (int)strlen(str);
	for (i = 0; i < lens; i++)
	{
		if (str[i] != c)str[size++] = str[i];
	}
	str[size] = '\0';
}
int main(int argc, char* argv[])
{
	int i, j, k;
	char d[50000], filei_prb[300], filei_dbd[300], filei_mot[300], filei_motifs_tfclass[300], file_sta[300];
	char fileo_fold[300], fileo_fold_best[300], fileo_polka[300];
	char fileo_rank_dbd[300], fileo_rank_motif_base[300], fileo_rank_motif[300], fileo_rank_dbd_best[300], fileo_rank_motif_best[300];
	char fileo_ratio_dbd[300], fileo_ratio_motif_base[300], fileo_ratio_motif[300], fileo_ratio_dbd_best[300], fileo_ratio_motif_best[300];
	
	FILE* in_prb, * in_dbd, * in_motif, * in_motif_tfcass;
	FILE* out_fold, * out_fold_best, * out_polka, * out_sta;
	FILE* out_rank_dbd, * out_rank_dbd_best, * out_rank_motif_best;
	FILE* out_ratio_dbd, * out_ratio_dbd_best, * out_ratio_motif_best;
	FILE** out_rank_motif, ** out_ratio_motif;

	if (argc != 20)
	{
		puts("Sintax: 1file in_table_prb 2file in_table_dbd file 3in_table_mot 4filei_TFClass_table"); //file in_table_per,
		puts("5int nrun(default 50) 6int min_group_size 7int step_group_size 8int print_headers");
		puts("9file_out_dynamics 10file_out_fold 11file_out_fold_best" );
		puts("12file_out_rank_dbd 13file_out_rank_dbd_best 14file_out_rank_motif_base 15file_out_rank_motif_best ");
		puts("16file_out_ratio_dbd 17file_out_ratio_dbd_best 18file_out_ratio_motif_base 19file_out_ratio_motif_best ");
		return -1;
	}
	//strcpy(filei_per, argv[1]);//in_file
	strcpy(filei_prb, argv[1]);//in_file
	strcpy(filei_dbd, argv[2]);//in_file
	strcpy(filei_mot, argv[3]);//in_file
	strcpy(filei_motifs_tfclass, argv[4]);
	int nrun = atoi(argv[5]); // no. of ga runs
	int min_gros = atoi(argv[6]);//min group size
	int step_gros = atoi(argv[7]);//step group size
	int print_headers = atoi(argv[8]);// 0 / 1 == print or not headers to file_table*best files 
	//                1    2    3    4    5     6    7     8
	char ratc[8][4] = { "50", "60", "70", "80", "85", "90", "95", "100" };
	double rat[8];	
	for (j = 0; j < 8; j++)rat[j] = 0;
	//dynamics & fold
	strcpy(fileo_polka, argv[9]);//out_file
	strcpy(fileo_fold, argv[10]);//out_file
	strcpy(fileo_fold_best, argv[11]);//out_file
	//rank
	strcpy(fileo_rank_dbd, argv[12]);//out_file
	strcpy(fileo_rank_dbd_best, argv[13]);//out_file
	strcpy(fileo_rank_motif_base, argv[14]);//out_file		
	strcpy(fileo_rank_motif_best, argv[15]);//out_file	
	//ratio
	strcpy(fileo_ratio_dbd, argv[16]);//out_file
	strcpy(fileo_ratio_dbd_best, argv[17]);//out_file
	strcpy(fileo_ratio_motif_base, argv[18]);//out_file		
	strcpy(fileo_ratio_motif_best, argv[19]);//out_file	

	strcpy(file_sta, "minimax_mot_rank.txt");//out_file
	int nrun1 = nrun - 1;
	int max_gros = 1 + nrun1 * step_gros;

	if ((in_motif_tfcass = fopen(filei_motifs_tfclass, "rt")) == NULL)
	{
		printf("Input file %s can't be opened\n", filei_motifs_tfclass);
		return -1;
	}
	int mtot = 0;
	fgets(d, sizeof(d), in_motif_tfcass);//header
	while (fgets(d, sizeof(d), in_motif_tfcass) != NULL)
	{
		if (*d != '\n' && *d != '\t')mtot++;
	}
	rewind(in_motif_tfcass);
	fgets(d, sizeof(d), in_motif_tfcass);//header	
	if (mtot > MMAX)
	{
		printf("Number of motifs %d above the upper limit %d", mtot, MMAX);
		exit(1);
	}
	grpa* allg;
	allg = new grpa[nrun];
	if (allg == NULL) { puts("Out of memory..."); exit(1); }
	grps selg[8];

	char** tf_names;
	char** class_names;
	char** motif_names;
	{
		const size_t lens = 50;
		char sep = '\t';
		tf_names = new char* [mtot];
		if (tf_names == NULL) { printf("Out of memory..."); return -1; };
		for (i = 0; i < mtot; i++)
		{
			tf_names[i] = new char[lens];
			if (tf_names[i] == NULL) { puts("Out of memory..."); exit(1); }
		}
		motif_names = new char* [mtot];
		if (motif_names == NULL) { printf("Out of memory..."); return -1; };
		for (i = 0; i < mtot; i++)
		{
			motif_names[i] = new char[lens];
			if (motif_names[i] == NULL) { puts("Out of memory..."); exit(1); }
		}
		class_names = new char* [mtot];
		if (class_names == NULL) { printf("Out of memory..."); return -1; };
		for (i = 0; i < mtot; i++)
		{
			class_names[i] = new char[lens];
			if (class_names[i] == NULL) { puts("Out of memory..."); exit(1); }
		}
		size_t sizemot = lens * sizeof(tf_names[0][0]);
		for (i = 0; i < mtot; i++)
		{
			memset(tf_names[i], '\0', sizeof(sizemot));
			memset(class_names[i], '\0', sizeof(sizemot));
			fgets(d, sizeof(d), in_motif_tfcass);
			DelChar(d, '\n');
			if (UnderStolStr(d, tf_names[i], sizemot, 0, sep) == NULL) { printf("Wrong format %s\n", filei_motifs_tfclass); return(-1); }
			if (UnderStolStr(d, class_names[i], sizemot, 1, sep) == NULL) { printf("Wrong format %s\n", filei_motifs_tfclass); return(-1); }
			if (UnderStolStr(d, motif_names[i], sizemot, 4, sep) == NULL) { printf("Wrong format %s\n", filei_motifs_tfclass); return(-1); }
		}
	}
	fclose(in_motif_tfcass);

	if ((in_prb = fopen(filei_prb, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", filei_prb);
		return -1;
	}
	if ((in_dbd = fopen(filei_dbd, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", filei_dbd);
		return -1;
	}
	if ((in_motif = fopen(filei_mot, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", filei_mot);
		return -1;
	}
	if ((out_fold = fopen(fileo_fold, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_fold);
		return -1;
	}
	int num_thr;
	if (nrun > 1)num_thr = 8;
	else num_thr = 1;
	int num_thr1 = num_thr - 1;
	out_polka = out_sta = NULL;
	out_rank_dbd = out_rank_dbd_best = out_rank_motif_best = NULL;
	out_ratio_dbd = out_ratio_dbd_best = out_ratio_motif_best = NULL;
	if (num_thr > 1)
	{
		if ((out_rank_dbd = fopen(fileo_rank_dbd, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_dbd);
			return -1;
		}
		if ((out_ratio_dbd = fopen(fileo_ratio_dbd, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_ratio_dbd);
			return -1;
		}
	}
	char omode[3];
	memset(omode, '\0', sizeof(omode));
	if (print_headers == 1)strcpy(omode, "wt");
	else strcpy(omode, "at");
	
	if (num_thr > 1)
	{
		if ((out_rank_dbd_best = fopen(fileo_rank_dbd_best, omode)) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_dbd_best);
			return -1;
		}
		if ((out_rank_motif_best = fopen(fileo_rank_motif_best, omode)) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_motif_best);
			return -1;
		}
		if ((out_ratio_dbd_best = fopen(fileo_ratio_dbd_best, omode)) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_ratio_dbd_best);
			return -1;
		}
		if ((out_ratio_motif_best = fopen(fileo_ratio_motif_best, omode)) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_ratio_motif_best);
			return -1;
		}
		if ((out_sta = fopen(file_sta, omode)) == NULL)
		{
			printf("Input file %s can't be opened!", file_sta);
			return -1;
		}
		if ((out_polka = fopen(fileo_polka, omode)) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_polka);
			return -1;
		}
	}
	if ((out_fold_best = fopen(fileo_fold_best, omode)) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_fold_best);
		return -1;
	}
	allg[0].gsi = min_gros;
	allg[0].irez = 1;
	//printf("AllRuns\n");
	//printf("%d %d\t", allg[0].irez, allg[0].gsi);	
	{
		if (num_thr > 1)
		{
			if (print_headers == 1)fprintf(out_polka, "\t%d", allg[0].gsi);
			for (i = 1; i < nrun; i++)
			{
				int i1 = i - 1;
				allg[i].gsi = allg[i1].gsi + step_gros;
				allg[i].irez = allg[i1].irez + 1;
				if (print_headers == 1)fprintf(out_polka, "\t%d", allg[i].gsi);
		//		printf("%d %d\t", allg[i].irez, allg[i].gsi);
			}
			if (print_headers == 1)fprintf(out_polka, "\n");
		}
	}
	//printf("\n");
	out_rank_motif = out_ratio_motif = NULL;
	if (num_thr > 1)
	{
		out_rank_motif = new FILE * [num_thr];
		if (out_rank_motif == NULL) { printf("Not  enough memory!"); return -1; }
		for (i = 0; i < num_thr; i++)
		{
			memset(fileo_rank_motif, '\0', sizeof(fileo_rank_motif));
			strcpy(fileo_rank_motif, fileo_rank_motif_base);
			strcat(fileo_rank_motif, "_");
			strcat(fileo_rank_motif, ratc[i]);			
			rat[i] = atof(ratc[i]);
			rat[i] /= 100;
			if ((out_rank_motif[i] = fopen(fileo_rank_motif, "wt")) == NULL)
			{
				printf("Out file %s can't be opened!\n", fileo_rank_motif);
				return -1;
			}
		}
		out_ratio_motif = new FILE * [num_thr];
		if (out_ratio_motif == NULL) { printf("Not  enough memory!"); return -1; }
		for (i = 0; i < num_thr; i++)
		{
			memset(fileo_ratio_motif, '\0', sizeof(fileo_ratio_motif));
			strcpy(fileo_ratio_motif, fileo_ratio_motif_base);
			strcat(fileo_ratio_motif, "_");
			strcat(fileo_ratio_motif, ratc[i]);
			rat[i] = atof(ratc[i]);
			rat[i] /= 100;
			if ((out_ratio_motif[i] = fopen(fileo_ratio_motif, "wt")) == NULL)
			{
				printf("Out file %s can't be opened!\n", fileo_ratio_motif);
				return -1;
			}
		}
	}
	char tab = '\t';
	char buf[50];
	//int col_prb_gros = max_gros + 4;	
	int col_prb_gros = 3;
	int col_prb_pr = col_prb_gros + 1;
	fgets(d, sizeof(d), in_prb);//header	
	char track_name[300];
	char track_name0[300];
	const size_t lens = 100;
	int irez_max = 0;
	int gsi_max = 0;
	double auprc_max = -1;
	double auprc_max_zero = -1;
	int irez_max_zero = 0;
	int gsi_max_zero = 0;
	double precision_thr = 0.5;
	int best_index = 0;
	int best_index_zero = 0;
	int decil_tested = 2;
	for (i = 0; i < nrun; i++)
	{
		if (fgets(d, sizeof(d), in_prb) == NULL)
		{
			printf("Input file %s reading error!", filei_prb);
			exit(1);
		}
		DelChar(d, '\n');
		if (i == 0)
		{
			size_t sizemot = lens * sizeof(track_name[0]);
			memset(track_name, '\0', sizeof(sizemot));
			memset(track_name0, '\0', sizeof(sizemot));
			if (UnderStolStr(d, track_name, sizemot, 0, tab) == NULL)
			{
				printf("Wrong format in file %s\n", filei_dbd);
				return(-1);
			}
			if (num_thr > 1)fprintf(out_polka, "%s", track_name);
			strcpy(track_name0, track_name);
		}
		double prec_sta = 0;
		int gom = 0;
		if (num_thr > 1)
		{
			for (j = 0; j < decil_tested; j++)
			{
				prec_sta = UnderStolDouble(d, col_prb_pr + j, buf, sizeof(buf), tab);
				if (prec_sta <= precision_thr)
				{
					gom++;
				}
				else break;
			}
		}
		double auprc1 = UnderStolDouble(d, 1, buf, sizeof(buf), tab);
		if (num_thr > 1)fprintf(out_polka, "\t%f", auprc1);
		allg[i].auprc = auprc1;
		if (gom < 2)
		{
			//int cur_gros = UnderStolInt(d, col_prb_gros, buf, sizeof(buf), tab);
			if (auprc1 > auprc_max)
			{
				auprc_max = auprc1;
				best_index = i;
				irez_max = allg[i].irez;
				gsi_max = allg[i].gsi;
			}
		}
		else
		{
			if (auprc1 > auprc_max_zero)
			{
				auprc_max_zero = auprc1;
				best_index_zero = i;
				irez_max_zero = allg[i].irez;
				gsi_max_zero = allg[i].gsi;
			}
		}
	}	
	allg[0].ratio = 1;
	for (i = 1; i < nrun; i++)
	{
		double ratio = (allg[i].auprc - allg[i - 1].auprc) / allg[0].auprc;
		if (ratio > 0)allg[i].ratio = ratio;
		else allg[i].ratio = 0;
	}
	double auprc_min = allg[0].auprc;
	double auprc_raz;
	if (auprc_max > -1)
	{
		auprc_raz = auprc_max - auprc_min;
		selg[num_thr1].auprc = auprc_max;
		selg[num_thr1].gsi = gsi_max;
		selg[num_thr1].irez = irez_max;
	}
	else
	{
		auprc_raz = auprc_max_zero - auprc_min;
		selg[num_thr1].auprc = auprc_max_zero;
		selg[num_thr1].gsi = gsi_max_zero;
		selg[num_thr1].irez = irez_max_zero;
		auprc_max = auprc_max_zero;
		best_index = best_index_zero;
	}
	selg[num_thr1].rat = 1;
	allg[0].rat = 0;
	for (i = 1; i < nrun; i++)
	{
		int i1 = i - 1;
		allg[i].rat = allg[i1].rat;
		double ratio = (allg[i].auprc - auprc_min) / auprc_raz;
		if (ratio > allg[i].rat)allg[i].rat = ratio;
	}
	/*printf("AUPRC Ratio 8\n");
	for (i = 0; i < num_thr; i++)
	{
		printf("%.2f\t", rat[i]);
	}
	printf("\n");
	printf("AUPRC Runs\n");
	for (i = 0; i < nrun; i++)
	{
		int i1 = i - 1;
		if (i1 % 10 == 0)printf("\n");
		printf("%.3f\t", allg[i].rat);
	}
	printf("\n");*/
	
	int nmots[8]; //no. of motifs, from 8 various ratios
	for (j = 0; j < num_thr; j++)nmots[j] = 0;
	for (j = 0; j < num_thr1; j++)selg[j].rat = 0;
	{
		int ix = 0;
		for (j = 0; j < num_thr1; j++)
		{
			for (i = ix; i < nrun; i++)
			{
				if (selg[j].rat == 0)
				{
					if (allg[i].rat >= rat[j] && (i == 0 || allg[i - 1].rat < rat[j]))
					{
						selg[j].rat = allg[i].rat;
						selg[j].gsi = allg[i].gsi;
						selg[j].irez = allg[i].irez;
						selg[j].auprc = allg[i].auprc;
				//		printf("Ratio %.2f\t%.2f\t%d\t%d\n", rat[j], selg[j].rat, selg[j].irez, selg[j].gsi);
						ix = i;
						break;
					}
				}								
			}			
		}
	}
	/*printf("AUPRC Selected\n");
	for (j = 0; j < num_thr; j++)printf("\t%.2f", selg[j].rat); printf("\n");
	printf("GroupIndex\n");
	for (j = 0; j < num_thr; j++)printf("\t%d", selg[j].irez); printf("\n");
	printf("GroupSize\n");
	for (j = 0; j < num_thr; j++)printf("\t%d", selg[j].gsi); printf("\n");
	*/
	if (num_thr > 1)
	{
		fprintf(out_polka, "\n");		
	}
	fclose(out_polka);
	fclose(in_prb);
	/*if (auprc_max == -1)
	{
		printf("Input file %s - maximal pAUPRC is not found!", filei_prb);
		exit(1);
	}*/
	fgets(d, sizeof(d), in_motif);
	if (print_headers == 1)
	{
		if (num_thr > 1)
		{
			fprintf(out_rank_motif_best, "%s", d);//header 1 2 3 4 5
			fprintf(out_ratio_motif_best, "%s", d);//header 1 2 3 4 5
		}
	}
	int nfamilies = 0;
	fgets(d, sizeof(d), in_dbd);//1st header numbers count of motifs in families	
	if (num_thr > 1)
	{
		fprintf(out_rank_dbd, "%s", d);
		fprintf(out_ratio_dbd, "%s", d);
	}
	fprintf(out_fold, "%s", d);
	if (print_headers == 1)
	{
		if (num_thr > 1)
		{
			fprintf(out_rank_dbd_best, "%s", d);
			fprintf(out_ratio_dbd_best, "%s", d);
		}
		fprintf(out_fold_best, "%s", d);
	}
	DelChar(d, '\n');
	{
		int hlen = (int)strlen(d);
		hlen--;
		for (i = 1; i < hlen; i++)
		{
			if (d[i] == tab)
			{
				int c = (int)d[i + 1];
				if (isdigit(c))nfamilies++;
			}
		}
	}
	int best_index1 = best_index + 1;
	int* motif_rank_min;
	motif_rank_min = new int[mtot];
	if (motif_rank_min == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < mtot; i++)motif_rank_min[i] = 0;
	int** motif_rank;
	motif_rank = new int* [best_index1];
	if (motif_rank == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < best_index1; i++)
	{
		motif_rank[i] = new int[mtot];
		if (motif_rank == NULL) { printf("Out of memory..."); return -1; };
	}
	for (j = 0; j < best_index1; j++)for (i = 0; i < mtot; i++)motif_rank[j][i] = 0;
	double** motif_ratio;
	motif_ratio = new double* [best_index1];
	if (motif_ratio == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < best_index1; i++)
	{
		motif_ratio[i] = new double[mtot];
		if (motif_ratio == NULL) { printf("Out of memory..."); return -1; };
	}
	for (j = 0; j < best_index1; j++)for (i = 0; i < mtot; i++)motif_ratio[j][i] = 0;
	double* family_ratio;
	family_ratio = new double[nfamilies];
	if (family_ratio == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < nfamilies; i++)family_ratio[i] = 0;
	int* family_rank;
	family_rank = new int[nfamilies];
	if (family_rank == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < nfamilies; i++)family_rank[i] = 0;
	int* family_exp;
	family_exp = new int[nfamilies];
	if (family_exp == NULL) { printf("Out of memory..."); return -1; };
	int family_sum = 0;
	{
		char buf[50];
		for (i = 0; i < nfamilies; i++)
		{
			family_exp[i] = UnderStolInt(d, 2 + i, buf, sizeof(buf), tab);
			family_sum += family_exp[i];
		}
	}
	fgets(d, sizeof(d), in_dbd);// 2nd header	
	if (num_thr > 1)
	{
		fprintf(out_rank_dbd, "%s", d);
		fprintf(out_ratio_dbd, "%s", d);
	}
	fprintf(out_fold, "%s", d);
	if (print_headers == 1)
	{
		if (num_thr > 1)
		{
			fprintf(out_rank_dbd_best, "%s", d);
			fprintf(out_ratio_dbd_best, "%s", d);
		}
		fprintf(out_fold_best, "%s", d);
	}	
	size_t sizemot = lens * sizeof(track_name[0]);
	for (j = 0; j < best_index1; j++)
	{
		fgets(d, sizeof(d), in_motif);
		DelChar(d, '\n');
		memset(track_name, '\0', sizeof(sizemot));
		if (UnderStolStr(d, track_name, sizemot, 0, tab) == NULL)
		{
			printf("Wrong format in file %s\n", filei_mot);
			return(-1);
		}
		int cur_gros = UnderStolInt(d, 1, buf, sizeof(buf), tab);
		//fprintf(out_rank_motif, "%s", track_name);		
		//fprintf(out_rank_motif, "\t%d", cur_gros);
		if (j == best_index)
		{
			if (num_thr > 1)
			{
				fprintf(out_rank_motif_best, "%s\t%d", track_name, cur_gros);
				fprintf(out_ratio_motif_best, "%s\t%d", track_name, cur_gros);
			}
		}
		for (i = 0; i < mtot; i++)
		{
			if (UnderStolStr(d, track_name, sizemot, 2 + i, tab) == NULL)
			{
				printf("Wrong format in file %s\n", filei_mot);
				return(-1);
			}
			int str_tr = (int)strlen(track_name);
			if (str_tr > 0)
			{
				if (j == 0)
				{
					motif_rank[j][i] = cur_gros;
					motif_ratio[j][i] = allg[j].ratio;
				}
				else
				{
					int j1 = j - 1;
					if (motif_rank[j1][i] > 0)
					{
						motif_rank[j][i] = motif_rank[j1][i];
						motif_ratio[j][i] = motif_ratio[j1][i];
					}
					else
					{
						motif_rank[j][i] = cur_gros;
						motif_ratio[j][i] = allg[j].ratio;
					}
				}
				if (motif_rank_min[i] == 0)motif_rank_min[i] = motif_rank[j][i];
			}
			else motif_rank[j][i] = 0;
			//fprintf(out_rank_motif, "\t");			
			/*if (motif_rank[i] > 0)
			{
				fprintf(out_rank_motif, "%d", motif_rank[i]);
			}*/
		}
		fgets(d, sizeof(d), in_dbd);
		DelChar(d, '\n');
		//size_t sizemot = lens * sizeof(track_name[0]);
		memset(track_name, '\0', sizeof(sizemot));
		if (UnderStolStr(d, track_name, sizemot, 0, tab) == NULL)
		{
			printf("Wrong format in file %s\n", filei_dbd);
			return(-1);
		}
		cur_gros = UnderStolInt(d, 1, buf, sizeof(buf), tab);
		if (num_thr > 1)
		{
			fprintf(out_rank_dbd, "%s\t%d", track_name, cur_gros);
			fprintf(out_ratio_dbd, "%s\t%d", track_name, cur_gros);
		}
		fprintf(out_fold, "%s\t%d", track_name, cur_gros);
		if (j == best_index)
		{
			if (num_thr > 1)
			{
				fprintf(out_rank_dbd_best, "%s\t%d", track_name, cur_gros);
				fprintf(out_ratio_dbd_best, "%s\t%d", track_name, cur_gros);
			}
			fprintf(out_fold_best, "%s\t%d", track_name, cur_gros);			
		}
		int family_obs;
		for (i = 0; i < nfamilies; i++)
		{
			if (UnderStolStr(d, track_name, sizemot, 2 + i, tab) == NULL)
			{
				printf("Wrong format in file %s\n", filei_dbd);
				return(-1);
			}
			double ratio = 0;
			int str_tr = (int)strlen(track_name);
			if (str_tr > 0)
			{
				family_obs = atoi(track_name);
				if (family_rank[i] == 0)
				{
					family_rank[i] = cur_gros;
					family_ratio[i] = allg[j].ratio;
				}
				ratio = ((double)family_obs / family_exp[i]) / ((double)cur_gros / family_sum);
			}
			else
			{
				family_rank[i] = 0;
				family_ratio[i] = 0;
			}
			if (num_thr > 1)
			{
				fprintf(out_rank_dbd, "\t");
				fprintf(out_ratio_dbd, "\t");
			}
			fprintf(out_fold, "\t");
			if (family_rank[i] > 0)
			{
				if (num_thr > 1)
				{
					fprintf(out_rank_dbd, "%d", family_rank[i]);
					fprintf(out_ratio_dbd, "%f", family_ratio[i]);
				}
				fprintf(out_fold, "%f", ratio);
			}
			if (j == best_index)
			{
				if (num_thr > 1)
				{
					fprintf(out_rank_dbd_best, "\t");
					fprintf(out_ratio_dbd_best, "\t");
				}
				fprintf(out_fold_best, "\t");
				if (family_rank[i] > 0)
				{
					if (num_thr > 1)
					{
						fprintf(out_rank_dbd_best, "%d", family_rank[i]);
						fprintf(out_ratio_dbd_best, "%f", family_ratio[i]);
					}
					fprintf(out_fold_best, "%f", ratio);
				}
			}
		}
		if (num_thr > 1)
		{
			fprintf(out_rank_dbd, "\n");
			fprintf(out_ratio_dbd, "\n");
		}
		fprintf(out_fold, "\n");
		if (j == best_index)
		{
			if (num_thr > 1)
			{
				fprintf(out_rank_dbd_best, "\n");
				fprintf(out_ratio_dbd_best, "\n");
			}
			fprintf(out_fold_best, "\n");
		}
	}
	if (num_thr > 1)
	{
		j = best_index;
		{
			for (i = 0; i < mtot; i++)
			{
				fprintf(out_rank_motif_best, "\t");
				fprintf(out_ratio_motif_best, "\t");
				if (motif_rank[j][i] > 0)
				{
					fprintf(out_rank_motif_best, "%d", motif_rank[j][i]);
					fprintf(out_ratio_motif_best, "%f", motif_ratio[j][i]);
				}
			}
			fprintf(out_rank_motif_best, "\n");
			fprintf(out_ratio_motif_best, "\n");
		}
		for (k = 0; k < num_thr; k++)
		{
			fprintf(out_rank_motif[k], "%s\n\t\t", track_name0);
			fprintf(out_ratio_motif[k], "%s\n\t\t", track_name0);
			{
				int x = min_gros;
				for (j = 0; j < selg[k].irez; j++)
				{
					fprintf(out_rank_motif[k], "\t%d", x);
					fprintf(out_ratio_motif[k], "\t%d", x);
					x += step_gros;
				}
			}
			fprintf(out_rank_motif[k], "\n");
			fprintf(out_ratio_motif[k], "\n");
		}
		{			
			for (i = 0; i < mtot; i++)
			{
				if (motif_rank_min[i] != 0)
				{
					//printf("Motif %d\tThr %d\n", i + 1, motif_rank_min[i]);
					for (k = num_thr1; k >=0; k--)
					{
						if (motif_rank_min[i] <= selg[k].gsi)
						{
							//printf("Rez %d\tThr %d\n", k + 1, selg[k].gsi);
							nmots[k]++;
							fprintf(out_rank_motif[k], "%s\t%s\t%s", class_names[i], tf_names[i], motif_names[i]);							
							fprintf(out_ratio_motif[k], "%s\t%s\t%s", class_names[i], tf_names[i], motif_names[i]);
							for (j = 0; j < selg[k].irez; j++)
							{
								fprintf(out_rank_motif[k], "\t");
								fprintf(out_ratio_motif[k], "\t");
								if (motif_rank[j][i] > 0)
								{
									fprintf(out_rank_motif[k], "%d", motif_rank[j][i]);
									fprintf(out_ratio_motif[k], "%f", motif_ratio[j][i]);
									//printf("%d %d\t", j + 1, motif_rank[j][i]);
								}
							}
						//	printf("\n");
							fprintf(out_rank_motif[k], "\n");
							fprintf(out_ratio_motif[k], "\n");
						}
						else break;
					}
				}
			}
		}
		for (k = 0; k < num_thr; k++)fclose(out_rank_motif[k]);
		for (k = 0; k < num_thr; k++)fclose(out_ratio_motif[k]);
		for (k = num_thr1; k >= 0; k--)
		{
			fprintf(out_sta, "%s\t%s\tAUPRC\t%f\tRatioAUPRC\t%f\tNMotifsMax\t%d\tNMotifsTotal\t%d\n", fileo_rank_motif_base, ratc[k], selg[k].auprc, selg[k].rat, selg[k].gsi, nmots[k]);
		}
		fclose(out_sta);
		fclose(out_rank_dbd);
		fclose(out_rank_dbd_best);
		fclose(out_rank_motif_best);
		fclose(out_ratio_dbd);
		fclose(out_ratio_dbd_best);
		fclose(out_ratio_motif_best);
	}
	fclose(out_fold);
	fclose(out_fold_best);
	fclose(in_dbd);
	delete[] family_exp;
	delete[] family_rank;
	delete[] family_ratio;
	delete[] motif_rank_min;
	delete[] out_rank_motif;
	for (i = 0; i < mtot; i++)
	{
		delete[] tf_names[i];
		delete[] motif_names[i];
	}
	for (i = 0; i < best_index1; i++)delete[] motif_ratio[i];
	delete[] motif_ratio;
	for (i = 0; i < best_index1; i++)delete[] motif_rank[i];
	delete[] motif_rank;
	delete[] tf_names;
	delete[] motif_names;
	delete[] allg;
	return 1;
}

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
// rnaseq table filtered for crireria in two columns log2Fold & p-adj to extract list of up- & down regulated DEGs & not DEGs
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
	int i, j;
	char d[50000], filei_prb[300], filei_dbd[300], filei_mot[300], filei_motifs_tfclass[300]; //filei_per[300],
	char fileo_rank_dbd[300], fileo_rank_motif[300], fileo_fold[300], fileo_rank_dbd_best[300], fileo_rank_motif_best[300], fileo_fold_best[300], fileo_polka[300];
	double precision_thr = 0.5;
	int best_index = 0, best_size = 1, decil_tested = 1;
	FILE* in_prb, * in_dbd, * in_motif, * in_motif_tfcass;
	FILE *out_rank_dbd, * out_rank_motif, * out_rank_dbd_best, * out_rank_motif_best, * out_fold, * out_fold_best, * out_polka;//* in_per,

	if (argc != 16)
	{
		puts("Sintax: 1file in_table_prb 2file in_table_dbd file 3in_table_mot 4filei_TFClass_table"); //file in_table_per,
		puts("5int nrun(default 50) 6int min_group_size 7int step_group_size 8int print_headers ");
		puts("9file_out_rank_dbd 10file_out_rank_motif 11file_out_fold");		
		puts("12file_out_rank_dbd_best 13file_out_rank_motif_best 14file_out_fold_best 15file_out_dynamics");
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
	//3
	strcpy(fileo_rank_dbd, argv[9]);//out_file
	strcpy(fileo_rank_motif, argv[10]);//out_file
	strcpy(fileo_fold, argv[11]);//out_file
	//3
	strcpy(fileo_rank_dbd_best, argv[12]);//out_file
	strcpy(fileo_rank_motif_best, argv[13]);//out_file
	strcpy(fileo_fold_best, argv[14]);//out_file
	//1
	strcpy(fileo_polka, argv[15]);//out_file
	int max_gros = nrun * step_gros;
	
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
	if ((out_rank_dbd = fopen(fileo_rank_dbd, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_rank_dbd);
		return -1;
	}
	if ((out_rank_motif = fopen(fileo_rank_motif, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_rank_motif);
		return -1;
	}
	if (print_headers == 1)
	{
		if ((out_rank_dbd_best = fopen(fileo_rank_dbd_best, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_dbd_best);
			return -1;
		}
		if ((out_rank_motif_best = fopen(fileo_rank_motif_best, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_motif_best);
			return -1;
		}
		if ((out_fold_best = fopen(fileo_fold_best, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_fold_best);
			return -1;
		}
		if ((out_polka = fopen(fileo_polka, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_polka);
			return -1;
		}
	}
	else
	{
		if ((out_rank_dbd_best = fopen(fileo_rank_dbd_best, "at")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_dbd_best);
			return -1;
		}
		if ((out_rank_motif_best = fopen(fileo_rank_motif_best, "at")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_motif_best);
			return -1;
		}
		if ((out_fold_best = fopen(fileo_fold_best, "at")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_fold_best);
			return -1;
		}
		if ((out_polka = fopen(fileo_polka, "at")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_polka);
			return -1;
		}
	}
	if (print_headers == 1)
	{		
		int cur_gros = min_gros;
		for (i = 0; i < nrun; i++)
		{
			fprintf(out_polka, "\t%d", cur_gros);
			cur_gros += step_gros;
		}
		fprintf(out_polka, "\n");
	}
	char tab = '\t';	
	char buf[50];
	//int col_prb_gros = max_gros + 4;	
	int col_prb_gros = 3;
	int col_prb_pr = col_prb_gros + 1;
	fgets(d, sizeof(d), in_prb);//header
	double auprc_max = -1;
	char track_name[300];
	char track_name0[300];	
	const size_t lens = 100;
	for (i = 0; i < nrun;i++)
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
			fprintf(out_polka, "%s", track_name);
			strcpy(track_name0, track_name);
		}
		double prec_sta = 0;
		int gom = 0;
		for (j = 0; j < decil_tested; j++)
		{
			prec_sta = UnderStolDouble(d, col_prb_pr + j, buf, sizeof(buf), tab);			
			if (prec_sta <= precision_thr)
			{
				gom++;				
			}
			else break;
		}
		double auprc = UnderStolDouble(d, 1, buf, sizeof(buf), tab);
		fprintf(out_polka, "\t%f", auprc);
		if (gom < 2)
		{
			int cur_gros = UnderStolInt(d, col_prb_gros, buf, sizeof(buf), tab);						
			if (auprc > auprc_max)
			{
				auprc_max = auprc;
				best_index = i + 1;
				best_size = cur_gros;
			}
		}
	}
	int best_index1 = best_index - 1;
	fprintf(out_polka, "\n");
	fclose(out_polka);
	fclose(in_prb);
	if(auprc_max == -1)
	{
		printf("Input file %s - maximal pAUPRC is not found!", filei_prb);
		exit(1);
	}		
	fgets(d, sizeof(d), in_motif);	
	if (print_headers == 1)
	{
		fprintf(out_rank_motif_best, "%s", d);//header 1 2 3 4 5
	}	
	int nfamilies = 0;
	fgets(d, sizeof(d), in_dbd);//1st header numbers count of motifs in families
	fprintf(out_rank_dbd, "%s", d);
	fprintf(out_fold, "%s", d);
	if (print_headers == 1)
	{
		fprintf(out_rank_dbd_best, "%s", d);
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
	int* motif_occ;
	motif_occ = new int[mtot];
	if (motif_occ == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < mtot; i++)motif_occ[i] = 0;
	int** motif_rank;
	motif_rank = new int*[best_index];
	if (motif_rank == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < best_index; i++)
	{
		motif_rank[i] = new int[mtot];
		if (motif_rank == NULL) { printf("Out of memory..."); return -1; };
	}	
	for (j = 0; j < best_index; j++)for (i = 0; i < mtot; i++)motif_rank[j][i] = 0;
	int* family_rank;
	family_rank = new int[nfamilies];
	if (family_rank == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < nfamilies; i++)family_rank[i] = 0;
	int *family_exp;
	family_exp = new int[nfamilies];	
	if (family_exp == NULL) { printf("Out of memory..."); return -1; };
	int family_sum = 0;
	{		
		char buf[50];
		for (i = 0; i < nfamilies; i++)
		{
			family_exp[i] = UnderStolInt(d, 2+i, buf, sizeof(buf), tab);
			family_sum += family_exp[i];
		}		
	}
	fgets(d, sizeof(d), in_dbd);// 2nd header	
	fprintf(out_rank_dbd, "%s", d);
	fprintf(out_fold, "%s", d);
	if (print_headers == 1)
	{
		fprintf(out_rank_dbd_best, "%s", d);
		fprintf(out_fold_best, "%s", d);
	}
	size_t sizemot = lens * sizeof(track_name[0]);
	for (j = 0; j < best_index; j++)
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
		if (j == best_index1)
		{			
			fprintf(out_rank_motif_best, "%s\t%d", track_name,cur_gros);
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
				if (j == 0)motif_rank[j][i] = cur_gros;
				else
				{
					if (motif_rank[j - 1][i] > 0)motif_rank[j][i] = motif_rank[j - 1][i];
					else motif_rank[j][i] = cur_gros;
				}
				motif_occ[i] = 1;
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
		fprintf(out_rank_dbd, "%s", track_name);
		fprintf(out_fold, "%s", track_name);		
		fprintf(out_rank_dbd, "\t%d", cur_gros);
		fprintf(out_fold, "\t%d", cur_gros);
		if (j == best_index1)
		{
			fprintf(out_rank_dbd_best, "%s", track_name);
			fprintf(out_fold_best, "%s", track_name);
			fprintf(out_rank_dbd_best, "\t%d", cur_gros);
			fprintf(out_fold_best, "\t%d", cur_gros);
		}
		int family_obs;
		for (i = 0; i < nfamilies; i++)
		{
			if (UnderStolStr(d, track_name, sizemot, 2+i, tab) == NULL)
			{
				printf("Wrong format in file %s\n", filei_dbd);
				return(-1);
			}
			double ratio = 0;
			int str_tr = (int)strlen(track_name);
			if (str_tr > 0)
			{
				family_obs = atoi(track_name);
				if (family_rank[i] == 0)family_rank[i] = cur_gros;
				ratio = ((double)family_obs / family_exp[i]) / ((double)cur_gros / family_sum);
			}
			else
			{
				family_rank[i] = 0;
			}			
			fprintf(out_rank_dbd, "\t");			
			fprintf(out_fold, "\t");
			if (family_rank[i] > 0)
			{
				fprintf(out_rank_dbd, "%d", family_rank[i]);
				fprintf(out_fold, "%f", ratio);
			}
			if (j == best_index1)
			{
				fprintf(out_rank_dbd_best, "\t");
				fprintf(out_fold_best, "\t");
				if (family_rank[i] > 0)
				{
					fprintf(out_rank_dbd_best, "%d", family_rank[i]);
					fprintf(out_fold_best, "%f", ratio);
				}
			}
		}
		fprintf(out_rank_dbd, "\n");
		fprintf(out_fold, "\n");
		if (j == best_index1)
		{
			fprintf(out_rank_dbd_best, "\n");
			fprintf(out_fold_best, "\n");
		}
	}
	j = best_index1;
	{
		for (i = 0; i < mtot; i++)
		{
			fprintf(out_rank_motif_best, "\t");
			if (motif_rank[j][i] > 0)
			{
				fprintf(out_rank_motif_best, "%d", motif_rank[j][i]);
			}			
		}
		fprintf(out_rank_motif_best, "\n");
	}
	fprintf(out_rank_motif, "%s\n\t\t",track_name0);
	{
		int x = min_gros;
		for (j = 0; j < best_index; j++)
		{
			fprintf(out_rank_motif, "\t%d", x);
			x += step_gros;
		}
	}
	fprintf(out_rank_motif, "\n");
	for (i = 0; i < mtot; i++)
	{
		if (motif_occ[i] == 1)
		{
			fprintf(out_rank_motif, "%s\t%s\t%s",class_names[i],tf_names[i],motif_names[i]);
			for (j = 0; j < best_index; j++)
			{
				fprintf(out_rank_motif, "\t");
				if(motif_rank[j][i] > 0)fprintf(out_rank_motif, "%d", motif_rank[j][i]);
			}
			fprintf(out_rank_motif, "\n");
		}
	}
	fclose(out_rank_dbd);
	fclose(out_rank_motif);
	fclose(out_fold);
	fclose(out_rank_dbd_best);
	fclose(out_fold_best);	
	//fclose(in_per);
	fclose(in_dbd);	
	delete[] family_exp;
	delete[] family_rank;
	delete[] motif_occ;
	for (i = 0; i < mtot; i++)
	{			
		delete[] tf_names[i];
		delete[] motif_names[i];
	}
	for (i = 0; i < best_index; i++)delete[] motif_rank[i];
	delete[] tf_names;
	delete[] motif_names;	
	delete[] motif_rank; 
	return 1;
}

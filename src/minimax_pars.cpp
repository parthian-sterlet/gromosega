#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <cmath>

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
	char d[20000], filei_prb[300], filei_dbd[300]; //filei_per[300],
	char fileo_rank_evo[300], fileo_fold_evo[300], fileo_rank_best[300], fileo_fold_best[300], fileo_polka[300];	
	double precision_thr = 0.5;
	int best_index = 0, best_size = 1, decil_tested = 2;
	FILE * in_prb, * in_dbd, * out_rank_evo, * out_rank_best, * out_fold_evo, * out_fold_best, *out_polka;//* in_per,

	if (argc != 12)
	{
		puts("Sintax: 1file in_table_prb 2file in_table_dbd "); //file in_table_per,
		puts("3int nrun(default 50) 4int min_group_size 5int step_group_size 6int print_headers ");
		puts("7file_out_table_rank_evo 8file_out_table_fold_evo 9file_out_table_rank_best 10file_out_table_fold_best 11file_out_polka");		
		return -1;
	}
	//strcpy(filei_per, argv[1]);//in_file
	strcpy(filei_prb, argv[1]);//in_file
	strcpy(filei_dbd, argv[2]);//in_file
	int nrun = atoi(argv[3]); // no. of ga runs
	int min_gros = atoi(argv[4]);//min group size
	int step_gros = atoi(argv[5]);//step group size
	int print_headers = atoi(argv[6]);// 0 / 1 == print or not headers to file_table*best files 
	strcpy(fileo_rank_evo, argv[7]);//out_file
	strcpy(fileo_fold_evo, argv[8]);//out_file
	strcpy(fileo_rank_best, argv[9]);//out_file
	strcpy(fileo_fold_best, argv[10]);//out_file
	strcpy(fileo_polka, argv[11]);//out_file
	int max_gros = nrun * step_gros;

	/*if ((in_per = fopen(filei_per, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", filei_per);
		return -1;
	}*/
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
	if ((out_fold_evo = fopen(fileo_fold_evo, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_fold_evo);
		return -1;
	}
	if ((out_rank_evo = fopen(fileo_rank_evo, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_rank_evo);
		return -1;
	}
	if (print_headers == 1)
	{
		if ((out_rank_best = fopen(fileo_rank_best, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_best);
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
		if ((out_rank_best = fopen(fileo_rank_best, "at")) == NULL)
		{
			printf("Input file %s can't be opened!", fileo_rank_best);
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
	int col_prb_gros = max_gros + 4;	
	int col_prb_pr = col_prb_gros + 1;
	fgets(d, sizeof(d), in_prb);//header
	double auprc_max = -1;
	char track_name[300];
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
			if (UnderStolStr(d, track_name, sizemot, 0, tab) == NULL)
			{
				printf("Wrong format in file %s\n", filei_dbd);
				return(-1);
			}
			fprintf(out_polka, "%s", track_name);
		}
		double prec_sta = 0;
		int gom = 0;
		for (j = 0; j < decil_tested; j++)
		{
			prec_sta = UnderStolDouble(d, col_prb_pr + j, buf, sizeof(buf), tab);			
			if (prec_sta <= precision_thr)
			{
				gom = 1;
				break;
			}
		}
		double auprc = UnderStolDouble(d, 1, buf, sizeof(buf), tab);
		fprintf(out_polka, "\t%f", auprc);
		if (gom == 0)
		{
			int cur_gros = UnderStolInt(d, col_prb_gros, buf, sizeof(buf), tab);						
			if (auprc > auprc_max)
			{
				auprc_max = auprc;
				best_index = i;
				best_size = cur_gros;
			}
		}
	}
	fprintf(out_polka, "\n");
	fclose(out_polka);
	fclose(in_prb);
	if(auprc_max == -1)
	{
		printf("Input file %s - maximal pAUPRC is not found!", filei_prb);
		exit(1);
	}
	int nfamilies = 0;
	fgets(d, sizeof(d), in_dbd);//header
	fprintf(out_rank_evo, "%s", d);
	fprintf(out_fold_evo, "%s", d);
	if (print_headers == 1)
	{
		fprintf(out_rank_best, "%s", d);
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
	fgets(d, sizeof(d), in_dbd);//header
	fprintf(out_rank_evo, "%s", d);
	fprintf(out_fold_evo, "%s", d);
	if (print_headers == 1)
	{
		fprintf(out_rank_best, "%s", d);
		fprintf(out_fold_best, "%s", d);
	}
	DelChar(d, '\n');	
	/*char** family_names;
	{			
		family_names = new char* [nfamilies];
		if (family_names == NULL) { printf("Out of memory..."); return -1; };
		for (i = 0; i < nfamilies; i++)
		{
			family_names[i] = new char[lens];
			if (family_names[i] == NULL) { puts("Out of memory..."); exit(1); }
		}
		size_t sizemot = lens * sizeof(family_names[0][0]);
		for (i = 0; i < nfamilies; i++)
		{			
			memset(family_names[i], '\0', sizeof(sizemot));			
			if (UnderStolStr(d, family_names[i], sizemot, 2+i, tab) == NULL) 
			{ 
				printf("Wrong format in file %s\n", filei_dbd); 
				return(-1); 
			}			
		}
	}*/		
	for (j = 0; j <= best_index; j++)
	{		
		fgets(d, sizeof(d), in_dbd);
		DelChar(d, '\n');
		size_t sizemot = lens * sizeof(track_name[0]);		
		memset(track_name, '\0', sizeof(sizemot));
		if (UnderStolStr(d, track_name, sizemot, 0, tab) == NULL)
		{
			printf("Wrong format in file %s\n", filei_dbd);
			return(-1);
		}			
		int cur_gros = UnderStolInt(d, 1, buf, sizeof(buf), tab);			
		fprintf(out_rank_evo, "%s", track_name);
		fprintf(out_fold_evo, "%s", track_name);		
		fprintf(out_rank_evo, "\t%d", cur_gros);
		fprintf(out_fold_evo, "\t%d", cur_gros);
		if (j == best_index)
		{
			fprintf(out_rank_best, "%s", track_name);
			fprintf(out_fold_best, "%s", track_name);
			fprintf(out_rank_best, "\t%d", cur_gros);
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
			fprintf(out_rank_evo, "\t");			
			fprintf(out_fold_evo, "\t");
			if (family_rank[i] > 0)
			{
				fprintf(out_rank_evo, "%d", family_rank[i]);
				fprintf(out_fold_evo, "%f", ratio);
			}
			if (j == best_index)
			{
				fprintf(out_rank_best, "\t");
				fprintf(out_fold_best, "\t");
				if (family_rank[i] > 0)
				{
					fprintf(out_rank_best, "%d", family_rank[i]);
					fprintf(out_fold_best, "%f", ratio);
				}
			}
		}
		fprintf(out_rank_evo, "\n");
		fprintf(out_fold_evo, "\n");
		if (j == best_index)
		{
			fprintf(out_rank_best, "\n");
			fprintf(out_fold_best, "\n");
		}
	}
	fclose(out_rank_evo);
	fclose(out_fold_evo);
	fclose(out_rank_best);
	fclose(out_fold_best);	
	//fclose(in_per);
	fclose(in_dbd);	
	delete[] family_exp;
	delete[] family_rank;
	/*for (i = 0; i < nfamilies; i++)
	{
		delete[] family_names[i];		
	}
	delete[] family_names; */
	return 1;
}

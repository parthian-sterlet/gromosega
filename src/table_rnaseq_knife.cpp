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
struct table_rna
{
	//	char id[50];
	int num;
	double padj;
	double log2fc;
};
int compare_padj(const void* X1, const void* X2)//descending
{
	struct table_rna* S1 = (struct table_rna*)X1;
	struct table_rna* S2 = (struct table_rna*)X2;
	double dif = S1->padj - S2->padj;
	if (dif > 0)return 1;
	else return -1;
	return 0;
}
int compare_fc(const void* X1, const void* X2)//ascending
{
	struct table_rna* S1 = (struct table_rna*)X1;
	struct table_rna* S2 = (struct table_rna*)X2;
	double abs1 = fabs(S1->log2fc);
	double abs2 = fabs(S2->log2fc);
	double dif = abs1 - abs2;
	if (dif > 0)return -1;
	else return 1;
	return 0;
}
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
double UnderStol(char* str, int nstol, char* ret, size_t size, char razd)
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
	char d[15000], filei_genelist[300], filei_rnaseq[300], fileo_up_base[300], fileo_do_base[300], fileo_no_base[300], fileo_up[300], fileo_do[300], fileo_no[300], fileo_sta[80];
	//	double *val;
	FILE* in_genelist, * in_rnaseq, * out_sta;
	FILE ** out_up, ** out_do, ** out_no;

	if (argc != 15)
	{
		puts("Sintax: 1file rnaseq table, 2,3,4int columns gene_id,log2Fold,padj 5file all_gene's_ID_list 6int columns gene_id ");
		puts("7int #DEGcount per segment (default 50) 8int number of segments (default 5) 9int int #notDEGmin (default 1000) 10double threshold padj (default 0.05)");
		puts("11file_out_upDEG(0,1) 12file_out_downDEG(0,1) 13file_out_noDEG(0,1) 14int reverse up/down (1yes 0no)");
		return -1;
	}
	strcpy(filei_rnaseq, argv[1]);//out_file
	int col_gene = atoi(argv[2]);
	int col_log2fold = atoi(argv[3]);
	int col_padj = atoi(argv[4]);
	strcpy(filei_genelist, argv[5]);//in_file	
	int col_genome = atoi(argv[6]);
	int size_deg_seg = atoi(argv[7]);
	int nseg_deg = atoi(argv[8]);
	int size_nedeg_min = atoi(argv[9]);
	double padj_thr = atof(argv[10]); //0.05;	
	strcpy(fileo_up_base, argv[11]);//out_file
	strcpy(fileo_do_base, argv[12]);//out_file	
	strcpy(fileo_no_base, argv[13]);//out_file
	int rev = atoi(argv[14]);

	int size_deg_min = size_deg_seg * nseg_deg;
	col_gene--;
	col_log2fold--;
	col_padj--;
	col_genome--;
	
	out_up = new FILE * [nseg_deg];
	if (out_up == NULL) { printf("Not  enough memory!"); return -1; }
	out_do = new FILE * [nseg_deg];
	if (out_do == NULL) { printf("Not  enough memory!"); return -1; }
	out_no = new FILE * [nseg_deg];
	if (out_no == NULL) { printf("Not  enough memory!"); return -1; }

	if ((in_rnaseq = fopen(filei_rnaseq, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", filei_rnaseq);
		return -1;
	}
	if ((in_genelist = fopen(filei_genelist, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!", filei_genelist);
		return -1;
	}	
	if (rev == 0)
	{
		for (j = 0; j < nseg_deg; j++)
		{
			char buf[3];
			memset(buf, '\0', sizeof(buf));
			snprintf(buf, sizeof(buf), "%d", j + 1);			
			strcpy(fileo_up, fileo_up_base);
			strcpy(fileo_do, fileo_do_base);
			strcpy(fileo_no, fileo_no_base);
			strcat(fileo_up, "_");
			strcat(fileo_do, "_");
			strcat(fileo_no, "_");
			strcat(fileo_up, buf);
			strcat(fileo_do, buf);
			strcat(fileo_no, buf);
			if ((out_up[j] = fopen(fileo_up, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!", fileo_up);
				return -1;
			}
			if ((out_do[j] = fopen(fileo_do, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!", fileo_do);
				return -1;
			}
			if ((out_no[j] = fopen(fileo_no, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!", fileo_no);
				return -1;
			}
		}
	}
	else
	{
		for (j = 0; j < nseg_deg; j++)
		{
			char buf[3];
			memset(buf, '\0', sizeof(buf));			
			snprintf(buf, sizeof(buf), "%d", j + 1);
			strcpy(fileo_up, fileo_up_base);
			strcpy(fileo_do, fileo_do_base);
			strcpy(fileo_no, fileo_no_base);
			strcat(fileo_up, "_");
			strcat(fileo_do, "_"); 
			strcat(fileo_no, "_");
			strcat(fileo_up, buf);
			strcat(fileo_do, buf);
			strcat(fileo_no, buf);
			if ((out_up[j] = fopen(fileo_do, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!", fileo_do);
				return -1;
			}
			if ((out_do[j] = fopen(fileo_up, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!", fileo_up);
				return -1;
			}
			if ((out_no[j] = fopen(fileo_no, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!", fileo_no);
				return -1;
			}
		}
	}
	strcpy(fileo_sta, "table_rnaseq_knife.txt");
	if ((out_sta = fopen(fileo_sta, "at")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_sta);
		return -1;
	}
	char tab = '\t';
	int n_genes = 0;
	//fgets(d, sizeof(d), in_genelist);//header
	while (fgets(d, sizeof(d), in_genelist) != NULL)
	{
		DelChar(d, '\n');
		char gene_id[100];
		memset(gene_id, 0, sizeof(gene_id));
		if (UnderStolStr(d, gene_id, sizeof(gene_id), col_genome, tab) == NULL)break;
		int cd = (int)gene_id[0];
		if (cd >= 65 && cd <= 90) // ASCD..XYZ
		{
			n_genes++;
			continue;
		}
		else break;
	}
	if (n_genes == 0)
	{
		printf("Genes not found! Input file %s reading error!\n%s", filei_genelist, d);
		return -1;
	}
	rewind(in_genelist);
	int* gen_type; //-1 0 1 10 down no up neither
	gen_type = new int[n_genes];
	if (gen_type == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < n_genes; i++)gen_type[i] = 1000;
	char** genes;
	{
		const size_t lens = 100;
		genes = new char* [n_genes];
		if (genes == NULL) { printf("Out of memory..."); return -1; };
		for (i = 0; i < n_genes; i++)
		{
			genes[i] = new char[lens];
			if (genes[i] == NULL) { puts("Out of memory..."); exit(1); }
		}
		size_t sizemot = lens * sizeof(genes[0][0]);
		for (i = 0; i < n_genes; i++)memset(genes[i], '\0', sizeof(sizemot));
		i = 0;
		while (fgets(d, sizeof(d), in_genelist) != NULL)
		{
			DelChar(d, '\n');
			char gene_id[100];
			memset(gene_id, 0, sizeof(gene_id));
			if (UnderStolStr(d, gene_id, sizeof(gene_id), col_genome, tab) == NULL)break;
			int cd = (int)gene_id[0];
			if (cd >= 65 && cd <= 90) // ASCD..XYZ
			{
				int dlen = (int)strlen(gene_id);
				strncpy(genes[i], gene_id, dlen);
				genes[i][dlen] = '\0';
				i++;
				continue;
			}
			else break;
		}
	}
	fclose(in_genelist);

	int n_str = 0, n_str_cor = 0;
	int ndegup_tot = 0, ndegdo_tot = 0, nnedeg_tot = 0;
	char resh = '#';
	//	double log2fold_thr2 = 0.321928094887362, log2fold_thr1 = -log2fold_thr2;//log2(1.25) = -log2(0.8)
	while (fgets(d, sizeof(d), in_rnaseq) != NULL)
	{
		if (*d != resh)break;
	}
	double padj_min = 0, padj_max = 1, log2fc_min = -1000, log2fc_max = 1000;
	while (fgets(d, sizeof(d), in_rnaseq) != NULL)
	{
		DelChar(d, '\n');
		DelChar(d, ' ');
		n_str++;
		char gene_id[50];
		memset(gene_id, 0, sizeof(gene_id));
		if (UnderStolStr(d, gene_id, sizeof(gene_id), col_gene, tab) == NULL)
		{
			printf("1st time Rnaseq table parsing failed! Input file %s reading error\n %s!", filei_rnaseq, d);
			return -1;
		}
		int inx = -1;
		for (i = 0; i < n_genes; i++)
		{
			if (strcmp(gene_id, genes[i]) == 0)
			{
				inx = i;
				break;
			}
		}
		if (inx >= 0)
		{
			char buf[100];
			memset(buf, '\0', sizeof(buf));
			double log2fc = UnderStol(d, col_log2fold, buf, sizeof(buf), tab);
			memset(buf, '\0', sizeof(buf));
			double padj = UnderStol(d, col_padj, buf, sizeof(buf), tab);
			if (log2fc == -10000 || padj == -10000)continue;
			if (padj < padj_min || padj > padj_max)continue;
			if (log2fc < log2fc_min || log2fc > log2fc_max)continue;
			n_str_cor++;
			if (padj >= padj_thr)nnedeg_tot++;
			else
			{
				if (log2fc > 0)ndegup_tot++;
				else ndegdo_tot++;
			}
		}
	}
	rewind(in_rnaseq);
	int ndeg_tot = ndegup_tot + ndegdo_tot;
	table_rna* rnaseq;
	rnaseq = new table_rna[n_str_cor];
	if (rnaseq == NULL) { puts("Out of memory..."); exit(1); }
	while (fgets(d, sizeof(d), in_rnaseq) != NULL)
	{
		if (*d != resh)break;
	}
	j = 0;
	while (fgets(d, sizeof(d), in_rnaseq) != NULL)
	{
		DelChar(d, '\n');
		DelChar(d, ' ');
		char gene_id[50];
		memset(gene_id, 0, sizeof(gene_id));
		if (UnderStolStr(d, gene_id, sizeof(gene_id), col_gene, tab) == NULL)
		{
			printf("2nd time Rnaseq table parsing failed! Input file %s reading error\n %s!", filei_rnaseq, d);
			return -1;
		}
		int inx = -1;
		for (i = 0; i < n_genes; i++)
		{
			if (strcmp(gene_id, genes[i]) == 0)
			{
				inx = i;
				break;
			}
		}
		if (inx >= 0)
		{
			rnaseq[j].num = inx;
			char buf[100];
			memset(buf, '\0', sizeof(buf));
			double log2fc = UnderStol(d, col_log2fold, buf, sizeof(buf), tab);
			memset(buf, '\0', sizeof(buf));
			double padj = UnderStol(d, col_padj, buf, sizeof(buf), tab);
			if (log2fc == -10000 || padj == -10000)continue;
			if (padj < padj_min || padj > padj_max)continue;
			if (log2fc < log2fc_min || log2fc > log2fc_max)continue;
			rnaseq[j].padj = padj;
			rnaseq[j].log2fc = log2fc;
			//memset(rnaseq[j].id, '\0', sizeof(rnaseq[j].id));
		//	strcpy(rnaseq[j].id, gene_id);
			j++;
		}
	}
	fclose(in_rnaseq);
	qsort(rnaseq, n_str_cor, sizeof(rnaseq[0]), compare_padj);
	qsort(rnaseq, ndeg_tot, sizeof(rnaseq[0]), compare_fc);
	qsort(&rnaseq[ndeg_tot], nnedeg_tot, sizeof(rnaseq[0]), compare_fc);
	double* thr_log2fc_up1, * thr_log2fc_do1;
	thr_log2fc_up1 = new double[nseg_deg];
	if (thr_log2fc_up1 == NULL) { printf("Out of memory..."); return -1; };
	thr_log2fc_do1 = new double[nseg_deg];
	if (thr_log2fc_do1 == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < nseg_deg; i++)thr_log2fc_up1[i] = 1000;
	for (i = 0; i < nseg_deg; i++)thr_log2fc_do1[i] = 1000;

	int take_up = 0, take_do = 0, take_no = 0;
	int print_up = 1, print_do = 1;
	int iup_last = -1, ido_last = -1;
	int count_up = 0, count_do = 0;
	for (i = 0; i < ndeg_tot; i++)
	{
		if (rnaseq[i].log2fc > 0)
		{
			if (print_up == 1)
			{
				take_up++;
				int rest = take_up % size_deg_seg;
				gen_type[rnaseq[i].num] = 1 + count_up;
				if (rest == 0)
				{
					thr_log2fc_up1[count_up++] = rnaseq[i].log2fc;					
				}
				iup_last = i;
				if (take_up == size_deg_min)print_up = 0;
			}
		}
		else
		{
			if (print_do == 1)
			{
				take_do++;
				int rest = take_do % size_deg_seg;				
				gen_type[rnaseq[i].num] = -1 - count_do;
				ido_last = i;
				if (rest == 0)
				{
					thr_log2fc_do1[count_do++] = rnaseq[i].log2fc;
				}
				if (take_do == size_deg_min)print_do = 0;
			}
		}
		if (print_up == 0 && print_do == 0)break;
	}
	double thr_log2fc_up = 10, thr_log2fc_do = -10;
	if (take_up > 0)
	{
		thr_log2fc_up = rnaseq[iup_last].log2fc;
	}
	if (take_do > 0)
	{
		thr_log2fc_do = rnaseq[ido_last].log2fc;
	}
	double thr_log2fc_nemin = 0, thr_log2fc_nemax = 0;
	for (i = n_str_cor - 1; i >= 0; i--)
	{
		double log2fc = rnaseq[i].log2fc;
		if (log2fc > thr_log2fc_do && log2fc < thr_log2fc_up)
		{
			if (log2fc > thr_log2fc_nemax)thr_log2fc_nemax = log2fc;
			if (log2fc < thr_log2fc_nemin)thr_log2fc_nemin = log2fc;
			take_no++;
			gen_type[rnaseq[i].num] = 0;
			if (take_no == size_nedeg_min)break;
		}
	}
	if (take_up > 0)
	{
		thr_log2fc_up = pow((double)2, thr_log2fc_up);
	}
	if (take_do > 0)
	{
		thr_log2fc_do = pow((double)2, thr_log2fc_do);
	}
	if (take_no > 0)
	{
		if (thr_log2fc_nemax > 0)thr_log2fc_nemax = pow((double)2, thr_log2fc_nemax);
		if (thr_log2fc_nemin < 0)thr_log2fc_nemin = pow((double)2, thr_log2fc_nemin);
	}
	for (i = 0; i < n_genes; i++)
	{
		int ii = gen_type[i];
		if (ii == 1000)// not found
		{
			for (j = 0; j < nseg_deg; j++)
			{
				fprintf(out_up[j], "0\n");
				fprintf(out_do[j], "0\n");
				fprintf(out_no[j], "0\n");
			}			
			continue;
		}
		if (ii == 0) //not deg
		{			
			for (j = 0; j < nseg_deg; j++)
			{
				fprintf(out_up[j], "0\n");
				fprintf(out_do[j], "0\n");
				fprintf(out_no[j], "1\n");
			}			
			continue;
		}
		if (ii >= 1) //up deg
		{
			int cell = ii - 1; 
			for (j = 0; j < nseg_deg; j++)
			{
				if(j == cell)fprintf(out_up[j], "1\n");
				else fprintf(out_up[j], "0\n");
				fprintf(out_do[j], "0\n");
				fprintf(out_no[j], "0\n");
			}						
			continue;
		}
		if (ii <= -1)//down deg
		{
			int cell = - ii - 1;
			for (j = 0; j < nseg_deg; j++)
			{
				if (j == cell)fprintf(out_do[j], "1\n");
				else fprintf(out_do[j], "0\n");
				fprintf(out_up[j], "0\n");
				fprintf(out_no[j], "0\n");
			}			
			continue;
		}				
	}
	for (j = 0; j < nseg_deg; j++)
	{
		fclose(out_up[j]);
		fclose(out_do[j]);
		fclose(out_no[j]);
	}
	for (j = 0; j < nseg_deg; j++)
	{
		if (thr_log2fc_up1[j] != 1000)thr_log2fc_up1[j] = pow((double)2, thr_log2fc_up1[j]);
		if (thr_log2fc_do1[j] != 1000)thr_log2fc_do1[j] = pow((double)2, thr_log2fc_do1[j]);
	}
	if (rev == 0)
	{
		fprintf(out_sta, "%s\t%s\t%s\t%s\t%s\tTotal genes\t%d\t%s\t%d\t%s\t%d\t%s\t%d\tUpFC\t%f\tDoFC\t%f\tNoFC\t%f\t%f\t", filei_rnaseq, filei_genelist, argv[7], argv[8], argv[9], n_genes, argv[10], take_up, argv[11], take_do, argv[12], take_no, thr_log2fc_up, thr_log2fc_do, thr_log2fc_nemax, thr_log2fc_nemin);
		fprintf(out_sta, "UPs");
		for (j = 0; j < nseg_deg; j++)fprintf(out_sta, "\t%f", thr_log2fc_up1[j]);
		fprintf(out_sta, "\tDOs");
		for (j = 0; j < nseg_deg; j++)fprintf(out_sta, "\t%f", thr_log2fc_do1[j]);
		fprintf(out_sta, "\n");
	}
	else
	{
		fprintf(out_sta, "%s\t%s\t%s\t%s\t%s\tReverse_up&down\t%d\tTotal genes\t%d\t%s\t%d\t%s\t%d\t%s\t%d\tUpFC\t%f\tDoFC\t%f\tNoFC\t%f\t%f\t", filei_rnaseq, filei_genelist, argv[7], argv[8], argv[9], rev, n_genes, argv[10], take_do, argv[11], take_up, argv[12], take_no, thr_log2fc_up, thr_log2fc_do, thr_log2fc_nemax, thr_log2fc_nemin);
		fprintf(out_sta, "UPs");
		for (j = 0; j < nseg_deg; j++)fprintf(out_sta, "\t%f", thr_log2fc_do1[j]);
		fprintf(out_sta, "\tDOs");
		for (j = 0; j < nseg_deg; j++)fprintf(out_sta, "\t%f", thr_log2fc_up1[j]);
		fprintf(out_sta, "\n");
	}
	fclose(out_sta);
	for (i = 0; i < n_genes; i++)
	{
		delete[] genes[i];
	}
	delete[] genes;
	delete[] gen_type;
	delete[] out_up;
	delete[] out_do;
	delete[] thr_log2fc_up1;
	delete[] thr_log2fc_do1;
	//delete[] rnaseq;
	return 1;
}


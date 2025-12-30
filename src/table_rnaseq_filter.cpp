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
	if(isdigit(cd) || cd == 45)return atof(ret);//0123456789 or -
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
	int i;
	char d[200], filei_genelist[300], filei_rnaseq[300], fileo_up[300], fileo_do[300], fileo_no[300], fileo_sta[80];
	//	double *val;
	FILE * in_genelist, * in_rnaseq, * out_sta, *out_up, * out_do, * out_no;

	if (argc != 13)
	{
		puts("Sintax: 1file rnaseq table, 2,3,4int columns gene_id,log2Fold,padj 5file all_gene's_ID_list 6int columns gene_id ");
		puts("7double threshold log2Fold up-/downDEG (default 2) 8double threshold log2Fold notDEG (default 1.25) 9double threshold padj (default 0.05)");
		puts("10file_out_upDEG(0,1) 11file_out_downDEG(0,1) 12file_out_noDEG(0,1)");
		return -1;
	}
	strcpy(filei_rnaseq, argv[1]);//out_file
	int col_gene = atoi(argv[2]);
	int col_log2fold = atoi(argv[3]);
	int col_padj = atoi(argv[4]);
	strcpy(filei_genelist, argv[5]);//in_file	
	int col_genome = atoi(argv[6]);
	double threh_deg = atof(argv[7]);// 2 -> log2(2) means 1
	double threh_nedeg = atof(argv[8]);//1.25 -> log2(1.25) = 0.32...
	double padj_thr = atof(argv[9]); //0.05;
	strcpy(fileo_up, argv[10]);//out_file
	strcpy(fileo_do, argv[11]);//out_file
	strcpy(fileo_no, argv[12]);//out_file
	col_gene--;
	col_log2fold--;
	col_padj--;
	col_genome--;

	threh_deg = log2(threh_deg);
	//threh_deg = threh_nedeg = 0.074;
	threh_nedeg = log2(threh_nedeg);
//	printf("Thresholds: %f %f\n", threh_deg, threh_nedeg);

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
	if ((out_up = fopen(fileo_up, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_up);
		return -1;
	}
	if ((out_do = fopen(fileo_do, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_do);
		return -1;
	}
	if ((out_no = fopen(fileo_no, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_no);
		return -1;
	}
	strcpy(fileo_sta, "table_rnaseq_filter.txt");
	if ((out_sta = fopen(fileo_sta, "at")) == NULL)
	{
		printf("Input file %s can't be opened!", fileo_sta);
		return -1;
	}
	char tab = '\t';
	int n_genes = 0;
	fgets(d, sizeof(d), in_genelist);//header
	while (fgets(d, sizeof(d), in_genelist) != NULL)
	{	
		DelChar(d, '\n');
		char gene_id[50];
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
	if(n_genes == 0)
	{
		printf("Input file %s reading error!", d);
		return -1;
	}
	rewind(in_genelist);
	fgets(d, sizeof(d), in_genelist);//header
	char** genes;
	{
		const size_t lens = 50;
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
			char gene_id[50];
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
	int* iup, * ido, * ino;
	iup = new int [n_genes];
	if (iup == NULL) { printf("Out of memory..."); return -1; };
	ido = new int[n_genes];
	if (ido == NULL) { printf("Out of memory..."); return -1; };
	ino = new int[n_genes];
	if (ino == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < n_genes; i++)iup[i] = ido[i] = ino[i] = 0;

	int n_str = 0, total_up =0, total_do = 0, total_no = 0;		
//	double log2fold_thr2 = 0.321928094887362, log2fold_thr1 = -log2fold_thr2;//log2(1.25) = -log2(0.8)
	fgets(d, sizeof(d), in_rnaseq);// header
	while (fgets(d, sizeof(d), in_rnaseq) != NULL)
	{
		DelChar(d, '\n');
		DelChar(d, ' ');
		n_str++;
		int take_up = 0, take_do =0, take_no =0;		
		char gene_id[50];
		memset(gene_id, 0, sizeof(gene_id));
		if (UnderStolStr(d, gene_id, sizeof(gene_id),col_gene, tab) == NULL)
		{
			printf("Input file %s reading error!", d);
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
			char buf[50];
			double test1 = UnderStol(d, col_log2fold, buf, sizeof(buf), tab);
			double test2 = UnderStol(d, col_padj, buf, sizeof(buf), tab);
			if (test1 == -10000 || test2 == -10000)continue;
			double log2fold = log2(test1), padj = test2;
			if (padj >= padj_thr)
			{
				if (log2fold > -threh_nedeg && log2fold < threh_nedeg)
				{
					ino[inx] = 1;
				//	printf("%s\t%d\tneDEG\tpadj\t%f\tlog2fold\t%f\n", gene_id, inx, padj, log2fold);
					total_no++;
				}
			}
			if (padj < padj_thr)
			{
				if (log2fold > threh_deg)
				{
					iup[inx] = 1;
			//		printf("%s\t%d\tupDEG\tpadj\t%f\tlog2fold\t%f\n", gene_id, inx, padj, log2fold);
					total_up++;
				}
				if (log2fold < -threh_deg)
				{
					ido[inx] = 1;
				//	printf("%s\t%d\tdoDEG\tpadj\t%f\tlog2fold\t%f\n", gene_id, inx, padj, log2fold);
					total_do++;
				}
			}
		}		
	}
	fclose(in_rnaseq);			
	for (i = 0; i < n_genes; i++)
	{
		fprintf(out_up, "%d\n", iup[i]);
		fprintf(out_do, "%d\n", ido[i]);
		fprintf(out_no, "%d\n", ino[i]);
	}
	fclose(out_up);
	fclose(out_no);
	fclose(out_do);
	fprintf(out_sta, "%s\t%s\t%s\t%s\t%s\tTotal genes %d\tupDEGs\t%d\tdownDEGs\t%d\tnotDEGs\t%d\n", filei_rnaseq,filei_genelist, argv[7], argv[8], argv[9],n_genes, total_up, total_do, total_no);
	fclose(out_sta);
	for (i = 0; i < n_genes; i++)
	{
		delete[] genes[i];
	}
	delete[] genes;
	delete[] iup;
	delete[] ido;
	delete[] ino;
	return 1;
}

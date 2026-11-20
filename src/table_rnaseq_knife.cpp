#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include  <time.h>
#include <cmath>

#define SPACESIZE 11
double uw[SPACESIZE][SPACESIZE];

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
//genes without names
struct table_gene
{
	int numg;
	int num;//index in RNAseq table DEGs or non-DEG(padj,FC)
	double log2fc;
	double padj;
	int gen_type; // -1 padj or FC absent or out of range, 0 non-DEG wrong by FC, 1 non-DEG(padj,FC), 2up, 3down, -100000 default		
	double moc;// c/g
	double dic[SPACESIZE - 1]; // 10 dibucl types
	int type; // 0 non-DEG without mono- di- content, >0 / <0 up-DEG / down-DEG
	int sim; // 1 non-DEGs similar by mono- di- content, 0 all others
	double rm_up; //mahalanobis distance
	double rm_do; //mahalanobis distance
	double rm2;
};
// rnaseq table filtered for crireria in two columns log2Fold & p-adj to extract list of up- & down regulated DEGs & not DEGs
struct table_rna
{
	int gen_type; // -1 padj or FC absent or out of range, 0 non-DEG wrong by FC, 1 non-DEG(padj,FC), 2up, 3down, -100000 default
	int num;
	double padj;
	double log2fc;
};
int compare_genetype(const void* X1, const void* X2)//descending
{
	struct table_gene* S1 = (struct table_gene*)X1;
	struct table_gene* S2 = (struct table_gene*)X2;
	int dif = S1->gen_type - S2->gen_type;
	if (dif < 0)return 1;
	else return -1;
	return 0;
}
int compare_numg(const void* X1, const void* X2)//descending
{
	struct table_gene* S1 = (struct table_gene*)X1;
	struct table_gene* S2 = (struct table_gene*)X2;
	int dif = S1->numg - S2->numg;
	if (dif > 0)return 1;
	else return -1;
	return 0;
}
int compare_padj(const void* X1, const void* X2)//descending
{
	struct table_gene* S1 = (struct table_gene*)X1;
	struct table_gene* S2 = (struct table_gene*)X2;
	double dif = S1->padj - S2->padj;
	if (dif > 0)return 1;
	else return -1;
	return 0;
}
int compare_fc_top(const void* X1, const void* X2)//ascending
{
	struct table_gene* S1 = (struct table_gene*)X1;
	struct table_gene* S2 = (struct table_gene*)X2;
	double abs1 = fabs(S1->log2fc);
	double abs2 = fabs(S2->log2fc);
	double dif = abs1 - abs2;	
	if (dif > 0)return -1;
	else
	{
		if(dif < 0)return 1;
		else
		{
			double difY = S1->padj - S2->padj;
			if (difY > 0)return 1;
			else return -1;
		}				
	}
	return 0;
}
int compare_fc_bottom(const void* X1, const void* X2)//ascending
{
	struct table_gene* S1 = (struct table_gene*)X1;
	struct table_gene* S2 = (struct table_gene*)X2;
	double abs1 = fabs(S1->log2fc);
	double abs2 = fabs(S2->log2fc);
	double dif = abs1 - abs2;
	if (dif < 0)return -1;
	else
	{
		if (dif > 0)return 1;
		{
			double difY = S1->padj - S2->padj;
			if (difY < 0)return 1;
			else return -1;
		}
	}
	return 0;
}
int compare_mah(const void* X1, const void* X2)//ascending
{
	struct table_gene* S1 = (struct table_gene*)X1;
	struct table_gene* S2 = (struct table_gene*)X2;
	double rm1 = S1->rm2;
	double rm2 = S2->rm2;
	double dif = rm1 - rm2;
	if (dif < 0)return -1;
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
void MinusStr1(int size, int j, int j0, double buf, double b[SPACESIZE][SPACESIZE])
{
	int k;
	b[j][j0] = 0;
	for (k = j0 + 1; k < size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
void MinusStr2(int size, int j, int j0, double buf, double b[SPACESIZE][SPACESIZE])
{
	int k;
	b[j][j0] = 0;
	for (k = j - 1; k >= j0; k--)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
void MinusStr(int size, int j, int j0, double buf, double b[SPACESIZE][SPACESIZE])
{
	int k;
	for (k = 0; k < size; k++)
	{
		if (b[j0][k] == 0)continue;
		b[j][k] -= buf * b[j0][k];
	}
}
int BackMat(int size, FILE* out)
{
	int i, j;
	double buf, b[SPACESIZE][SPACESIZE];
	if (size == 1)
	{
		uw[0][0] = 1 / uw[0][0];
		return 1;
	}
	for (j = 0; j < size; j++)
	{
		for (i = 0; i < size; i++)
		{
			b[i][j] = 0;
		}
	}
	for (i = 0; i < size; i++)b[i][i] = 1;
	for (j = 0; j < size - 1; j++)
	{
		for (i = j + 1; i < size; i++)
		{
			if (uw[i][j] == 0)continue;
			buf = uw[i][j] / uw[j][j];
			MinusStr1(size, i, j, buf, uw);
			MinusStr(size, i, j, buf, b);
			if (fabs(buf) < 1e-035)
			{
				fprintf(out, "\n1back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//fprintf(out,"%g(%d)\n",uw[j][j],j);
				return -1;
				//
			}
			uw[i][j] = 0;
		}
	}
	for (j = size - 1; j > 0; j--)
	{
		for (i = j - 1; i >= 0; i--)
		{
			if (uw[i][j] == 0)continue;
			buf = uw[i][j] / uw[j][j];
			MinusStr2(size, i, j, buf, uw);
			MinusStr(size, i, j, buf, b);
			if (fabs(buf) < 1e-035)
			{
				fprintf(out, "\n2back %g(%d,%d) bij%g(%d,%d)\t", uw[i][j], i, j, b[i][j], i, j);
				//fprintf(out,"%g(%d)\n",uw[j][j],j);
				//
				return -1;
			}
			uw[i][j] = 0;
		}
	}
	for (i = 0; i < size; i++)
	{
		buf = uw[i][i];
		if (fabs(buf) < 1e-030)
		{
			fprintf(out, "\nback3 %g %d\t", uw[i][i], i);
			return -1;
		}
		for (j = 0; j < size; j++)
		{
			uw[i][j] /= buf;
			b[i][j] /= buf;
		}
	}
	/*
	for(i=0;i<size;i++)
	{

	for(j=0;j<size;j++)
	{
	fprintf(out,"%.2e ",uw[i][j]);
	}
	fprintf(out,"\n");
	}
	fprintf(out,"\n");
	for(i=0;i<size;i++)
	{
	for(j=0;j<size;j++)
	{
	fprintf(out,"%.1e ",b[i][j]);
	}
	fprintf(out,"\n");
	}
	fprintf(out,"\n");
	int n;
	for(i=0;i<size;i++)
	{
	for(j=0;j<size;j++)
	{
	double x=0;
	for(n=0;n<size;n++)x+=uw0[j][n]*b[n][i];
	fprintf(out,"%.1e ",x);
	}
	fprintf(out,"\n");
	}
	*/
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			uw[i][j] = b[i][j];
		}
	}
	return 1;
}
void Mix(int* a, int* b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
}
void Mix(double* a, double* b)
{
	double buf = *a;
	*a = *b;
	*b = buf;
}
double Laplace(double z)
{
	const double pi = 3.141592653589793;
	const double step = 0.0001;
	const double step2 = step / 2;
	double x, y = 0, dy;
	double prec = 1E-20;
	x = z;
	double mnoj = sqrt(2 * pi);
	do
	{
		x += step2;
		dy = step * exp(-x * x / 2);
		dy /= mnoj;
		y += dy;
		x += step2;
	} while (dy > prec);
	return(y);
}
double BackLaplace(double u)
{
	const double pi = 3.141592653589793;
	const double step = 0.0001;
	double x, y;
	if (u <= 0) return(-10000);
	if (u >= 1) return(-10000);
	x = 0;
	y = 0.5;
	double mnoj = 1 / sqrt(2 * pi);	
	if (u < 0.5)
	{
		while (y > u)
		{
			x -= step / 2;
			y -= step * mnoj * exp(-x * x / 2);
			x -= step / 2;
		}
	}
	else
	{
		while (y < u)
		{
			x += step / 2;
			y += step * mnoj * exp(-x * x / 2);
			x += step / 2;
		}
	}	
	return(x);	
}
int main(int argc, char* argv[])
{
	int i, j, k;
	char d[15000], dmi[300], ddi[300];
	char filei_genelist[300], filei_rnaseq[300], fileo_up_base[300], fileo_do_base[300], fileo_no_base[300], fileo_up[300], fileo_do[300], fileo_no[300], fileo_sta[80];
	char filei_moc[300], filei_dic[300], file_log[300];
	FILE* in_genelist, * in_rnaseq, * in_moc, * in_dic;
	FILE* out_sta;
	FILE** out_up, ** out_do, ** out_no, * outlog;

	if (argc != 20)
	{
		puts("Sintax: 1file rnaseq table, 2,3,4int columns gene_id,log2Fold,padj 5file all_gene's_ID_list 6int columns gene_id ");
		puts("7int #DEGcount per segment (default 50) 8int number of segments (default 5) 9int fraction {#notDEG (#DEGs)mono-dicontent,p-adj &FC}/{#notDEG (#DEGs)p-adj &FC}) 10double p-adj (default 0.05)");
		puts("11file_out_upDEG(0,1) 12file_out_downDEG(0,1) 13file_out_noDEG(0,1) 14int reverse up/down (1yes 0no) 15int 1 / -1 take DEGs from the top / bottom of |log2fold|");
		puts("16int randomization 0no=default 1yes 17file_in_mono_content 18file_in_di_content 19 #dinucleotide types 1CG 4AT,TA,CG,GC 10All10");
		return -1;
	}
	strcpy(filei_rnaseq, argv[1]);//out_file
	int col_gene = atoi(argv[2]);
	int col_log2fold = atoi(argv[3]);
	int col_padj = atoi(argv[4]);
	strcpy(filei_genelist, argv[5]);//in_file	
	int col_genome = atoi(argv[6]);
	int size_deg_seg = atoi(argv[7]);
	int nseq_deg = atoi(argv[8]);
	double size_nedeg_fraction = atof(argv[9]);// {no. of neDESs p-adj, FC, and adopted mono- di-content}/{no. of neDESs by p-adj & FC}
	double padj_thr = atof(argv[10]); //0.05;	
	strcpy(fileo_up_base, argv[11]);//out_file
	strcpy(fileo_do_base, argv[12]);//out_file	
	strcpy(fileo_no_base, argv[13]);//out_file
	int rev = atoi(argv[14]);
	int border = atoi(argv[15]); // 1 / -1 take DEGs with highest / lowest |log2fold| values	
	int random = atoi(argv[16]); // 0no / 1yes permute RNA-seq genes with present padj & log2fold 	
	strcpy(filei_moc, argv[17]);//in_file	
	strcpy(filei_dic, argv[18]);//in_file		
	strcpy(file_log, "backmat_err.txt");
	int ndity = atoi(argv[19]);
	double rat_nf_nb_max = 3;
	int back_err = 0;
	double z_sco_thr = BackLaplace(size_nedeg_fraction);
	if (z_sco_thr == -10000)
	{
		printf("Wrong size_nedeg_fraction = %f parameter, only values (0, 1) are allowed!\n", size_nedeg_fraction);
		return -1;
	}
	{
		int allowv[3] = { 1,4,10 };
		int ndity_test = 0;
		for (i = 0; i < 3; i++)
		{
			if (ndity == allowv[i]) {
				ndity_test = 1; break;
			}			
		}
		if (ndity_test == 0)
		{
			printf("Wrong number of dinucleotide types %d, only 1, 4 and 10 are allowed for {CG}, {AT,TA,CG,GC} and All 10 dinicleotidr complementary pairs\n", ndity);
			exit(1);
		}
	}	
	int space_size = 1 + ndity, space_size1 = ndity;

	int size_deg_min = size_deg_seg * nseq_deg;
	col_gene--;
	col_log2fold--;
	col_padj--;
	col_genome--;

	out_up = new FILE * [nseq_deg];
	if (out_up == NULL) { printf("Not  enough memory!"); return -1; }
	out_do = new FILE * [nseq_deg];
	if (out_do == NULL) { printf("Not  enough memory!"); return -1; }
	out_no = new FILE * [nseq_deg];
	if (out_no == NULL) { printf("Not  enough memory!"); return -1; }

	if ((in_rnaseq = fopen(filei_rnaseq, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", filei_rnaseq);
		return -1;
	}
	if ((in_genelist = fopen(filei_genelist, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", filei_genelist);
		return -1;
	}
	if ((in_moc = fopen(filei_moc, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", filei_moc);
		return -1;
	}
	if ((in_dic = fopen(filei_dic, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", filei_dic);
		return -1;
	}
	if (rev == 0)
	{
		for (j = 0; j < nseq_deg; j++)
		{
			strcpy(fileo_up, fileo_up_base);
			strcpy(fileo_do, fileo_do_base);
			strcpy(fileo_no, fileo_no_base);
			if (nseq_deg > 1)
			{
				char buf[3];
				memset(buf, '\0', sizeof(buf));
				snprintf(buf, sizeof(buf), "%d", j + 1);
				strcat(fileo_up, "_");
				strcat(fileo_do, "_");
				strcat(fileo_no, "_");
				strcat(fileo_up, buf);
				strcat(fileo_do, buf);
				strcat(fileo_no, buf);
			}
			if ((out_up[j] = fopen(fileo_up, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!\n", fileo_up);
				return -1;
			}
			if ((out_do[j] = fopen(fileo_do, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!\n", fileo_do);
				return -1;
			}
			if ((out_no[j] = fopen(fileo_no, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!\n", fileo_no);
				return -1;
			}
		}
	}
	else
	{
		for (j = 0; j < nseq_deg; j++)
		{
			strcpy(fileo_up, fileo_up_base);
			strcpy(fileo_do, fileo_do_base);
			strcpy(fileo_no, fileo_no_base);
			if (nseq_deg > 1)
			{
				char buf[3];
				memset(buf, '\0', sizeof(buf));
				snprintf(buf, sizeof(buf), "%d", j + 1);
				strcat(fileo_up, "_");
				strcat(fileo_do, "_");
				strcat(fileo_no, "_");
				strcat(fileo_up, buf);
				strcat(fileo_do, buf);
				strcat(fileo_no, buf);
			}
			if ((out_up[j] = fopen(fileo_do, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!\n", fileo_do);
				return -1;
			}
			if ((out_do[j] = fopen(fileo_up, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!\n", fileo_up);
				return -1;
			}
			if ((out_no[j] = fopen(fileo_no, "wt")) == NULL)
			{
				printf("Input file %s can't be opened!\n", fileo_no);
				return -1;
			}
		}
	}
	strcpy(fileo_sta, "table_rnaseq_knife.txt");
	{
		char arexp[] = "Gene ID";
		char resh = '#';
		while (fgets(d, sizeof(d), in_rnaseq) != NULL)
		{
			if (strncmp(d, arexp, 7) == 0)continue;
			if (*d != resh)break;
		}
	}
	int ngene_rnaseq = 0;
	char tab = '\t';
	do
	{
		DelChar(d, '\n');
		DelChar(d, ' ');
		char gene_id[50];
		memset(gene_id, 0, sizeof(gene_id));
		if (UnderStolStr(d, gene_id, sizeof(gene_id), col_gene, tab) == NULL)
		{
			printf("1st time Rnaseq table parsing failed! Input file %s reading error\n %s!", filei_rnaseq, d);
			return -1;
		}
		ngene_rnaseq++;
	} while (fgets(d, sizeof(d), in_rnaseq) != NULL);
	rewind(in_rnaseq);

	table_rna* rnaseq;
	rnaseq = new table_rna[ngene_rnaseq];
	if (rnaseq == NULL) { puts("Out of memory..."); exit(1); }
	int bad_gen_type = -10000;
	for (i = 0; i < ngene_rnaseq; i++)rnaseq[i].gen_type = bad_gen_type;
	int ngen_absent = 0, ngen_present = 0;
	double padj_min = 0, padj_max = 1, log2fc_min = -1000, log2fc_max = 1000, log2fc_do_max = -1000, log2fc_up_min = 1000;

	{
		char arexp[] = "Gene ID";
		char resh = '#';
		while (fgets(d, sizeof(d), in_rnaseq) != NULL)
		{
			if (strncmp(d, arexp, 7) == 0)continue;
			if (*d != resh)break;
		}
	}
	char** genes;
	{
		const size_t lens = 100;
		genes = new char* [ngene_rnaseq];
		if (genes == NULL) { printf("Out of memory..."); return -1; };
		for (i = 0; i < ngene_rnaseq; i++)
		{
			genes[i] = new char[lens];
			if (genes[i] == NULL) { puts("Out of memory..."); exit(1); }
		}
		size_t sizemot = lens * sizeof(genes[0][0]);
		for (i = 0; i < ngene_rnaseq; i++)memset(genes[i], '\0', sizeof(sizemot));
	}
	i = 0;
	do
	{
		DelChar(d, '\n');
		char gene_id[100];
		memset(gene_id, 0, sizeof(gene_id));
		if (UnderStolStr(d, gene_id, sizeof(gene_id), col_gene, tab) == NULL)break;
		int cd = (int)gene_id[0];
		if (cd >= 65 && cd <= 90) // ASCD..XYZ
		{
			int dlen = (int)strlen(gene_id);
			strncpy(genes[i], gene_id, dlen);
			genes[i][dlen] = '\0';
			//printf("%d\t%s\n", i + 1, genes[i]);
			char buf[100];
			memset(buf, '\0', sizeof(buf));
			double log2fc = UnderStolDouble(d, col_log2fold, buf, sizeof(buf), tab);
			if (log2fc != -10000)rnaseq[i].log2fc = log2fc;
			memset(buf, '\0', sizeof(buf));
			double padj = UnderStolDouble(d, col_padj, buf, sizeof(buf), tab);
			if (padj != -10000)rnaseq[i].padj = padj;
			if (log2fc == -10000 || padj == -10000)
			{
				ngen_absent++;
				rnaseq[i].padj = rnaseq[i].log2fc = bad_gen_type;
				rnaseq[i].gen_type = -1;
			}
			else
			{
				ngen_present++;
				if (padj < padj_thr)
				{
					if (log2fc > 0)
					{
						if (log2fc < log2fc_up_min)log2fc_up_min = log2fc;
						rnaseq[i].gen_type = 2;
					}
					else
					{
						if (log2fc > log2fc_do_max)log2fc_do_max = log2fc;
						rnaseq[i].gen_type = 3;
					}
				}
				else rnaseq[i].gen_type = 0;
			}
			i++;
			continue;
		}
		else break;
	} while (fgets(d, sizeof(d), in_rnaseq) != NULL);
	fclose(in_rnaseq);
	int ngen_rna_good = i;
	for (i = 0; i < ngen_rna_good; i++)
	{

		if (rnaseq[i].gen_type == 0)
		{
			if (log2fc_do_max <= rnaseq[i].log2fc && rnaseq[i].log2fc <= log2fc_up_min)
			{
				rnaseq[i].gen_type = 1;
			}
		}
	}	
	int n_genes_tot = 0;
	while (fgets(d, sizeof(d), in_genelist) != NULL)
	{
		DelChar(d, '\n');
		char gene_id[100];
		memset(gene_id, 0, sizeof(gene_id));
		if (UnderStolStr(d, gene_id, sizeof(gene_id), col_genome, tab) == NULL)break;
		int cd = (int)gene_id[0];
		if (cd >= 65 && cd <= 90) // ASCD..XYZ
		{
			n_genes_tot++;
			int inx = -1;
			for (i = 0; i < ngen_rna_good; i++)
			{
				if (strcmp(gene_id, genes[i]) == 0)
				{
					inx = i;
					break;
				}
			}
			if (inx < 0)rnaseq[i].gen_type = bad_gen_type;
			rnaseq[i].num = inx;
		}
		else break;
	}	
	rewind(in_genelist);
	
	if (random == 1)
	{
		srand((unsigned)time(NULL));
		int ngen_rna_good1 = ngen_rna_good - 1;
		for (i = 0; i < ngen_rna_good1; i++)
		{
			if (rnaseq[i].gen_type >= 1 && rnaseq[i].num >= 0)
			{
				int trys = 20;
				int i2 = -1;
				for (j = 0; j < trys; j++)
				{
					int test = 1 + i + (rand() % (ngen_rna_good1 - i));
					if (rnaseq[test].gen_type >= 1 && rnaseq[test].num >= 0)
					{
						i2 = test;
						break;
					}
				}
				if (i2 != -1)
				{
					Mix(&rnaseq[i].gen_type, &rnaseq[i2].gen_type);
					Mix(&rnaseq[i].padj, &rnaseq[i2].padj);
					Mix(&rnaseq[i].log2fc, &rnaseq[i2].log2fc);
					//char gene_id[100];
					//memset(gene_id, '\0', sizeof(gene_id));
					//int g1 = (int)strlen(genes[i]);
					//int g2 = (int)strlen(genes[i2]);
					//strcpy(gene_id, genes[i]);
					//gene_id[g1] = '\0';
					//strncpy(genes[i], genes[i2], g2);
					//genes[i][g2] = '\0';
					//strncpy(genes[i2], gene_id, g1);
					//genes[i2][g1] = '\0';
				}
			}
		}
	}

	int unknown_gen_type = -100000;
	table_gene* gene;
	gene = new table_gene[n_genes_tot];
	if (gene == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < n_genes_tot; i++)
	{
		gene[i].gen_type = bad_gen_type;
		gene[i].log2fc = gene[i].padj = -10000;
		gene[i].num = -1;
		gene[i].type = gene[i].sim = 0;
		gene[i].moc = 0.5;
		for (j = 0; j < 10; j++)gene[i].dic[j] = 0.1;
		gene[i].rm2 = gene[i].rm_do = gene[i].rm_up = 1;
	}
	int ndegup = 0, ndegdo = 0, nnedeg_good = 0, nnedeg_bad = 0, other_bad = 0;
	int twin[6] = { 0,1,2,4,5,7 };
	int numg = 0;
	j = 0;
	while (fgets(d, sizeof(d), in_genelist) != NULL)
	{
		if (fgets(dmi, sizeof(dmi), in_moc) == NULL)
		{
			printf("Input file %s Error!", filei_moc);
			return -1;
		}
		if (fgets(ddi, sizeof(ddi), in_dic) == NULL)
		{
			printf("Input file %s Error!", filei_dic);
			return -1;
		}
		DelChar(d, '\n');
		gene[j].numg = j;
		char gene_id[100];
		memset(gene_id, 0, sizeof(gene_id));
		if (UnderStolStr(d, gene_id, sizeof(gene_id), col_genome, tab) == NULL)break;
		int cd = (int)gene_id[0];
		if (cd >= 65 && cd <= 90) // ASCD..XYZ
		{
			int inx = -1;
			for (i = 0; i < ngene_rnaseq; i++)
			{
				if (strcmp(gene_id, genes[i]) == 0)
				{
					inx = i;
					break;
				}
			}
			gene[j].num = inx;
			if (inx >= 0 && rnaseq[inx].gen_type >= 0)
			{
				gene[j].gen_type = rnaseq[inx].gen_type;
				if (gene[j].gen_type >= 0)
				{
					gene[j].log2fc = rnaseq[inx].log2fc;
					gene[j].padj = rnaseq[inx].padj;
				}
				switch (gene[j].gen_type)
				{
				case 0: { nnedeg_bad++; break; }
				case 1: { nnedeg_good++; break; }
				case 2: { ndegup++; break; }
				case 3: { ndegdo++; break; }
				default: {
					other_bad++;
				}
				}
			}
			else
			{
				other_bad++;
			}
			char buf[30];
			int mo[2], di[10];
			memset(buf, '\0', sizeof(buf));
			for (k = 0; k < 2; k++)
			{
				mo[k] = UnderStolInt(dmi, k, buf, sizeof(buf), tab);
			}
			memset(buf, '\0', sizeof(buf));
			for (k = 0; k < 10; k++)
			{
				di[k] = UnderStolInt(ddi, k, buf, sizeof(buf), tab);
			}
			int len_prom = mo[0] + mo[1];
			gene[j].moc = (double)mo[1] / len_prom;
			double mocf[2];
			for (k = 0; k < 2; k++)mocf[k] = (double)mo[k] / len_prom;
			double dexp[10];// di  0aa/tt 1ac/gt 2ag/ct 3at 4ca/tg 5cc/gg 6cg 7ga/tc 8gc 9ta .., popolam 012 457
			dexp[0] = dexp[3] = dexp[9] = mocf[0] * mocf[0];
			dexp[5] = dexp[6] = dexp[8] = mocf[1] * mocf[1];
			dexp[1] = dexp[2] = dexp[4] = dexp[7] = mocf[0] * mocf[1];
			int len_prom1 = len_prom - 1;
			if (ndity == 10)
			{
				for (k = 0; k < 10; k++)gene[j].dic[k] = (double)di[k] / len_prom1 / dexp[k];
				for (k = 0; k < 6; k++)gene[j].dic[twin[k]] /= 2;
			}
			else
			{
				if (ndity == 1)
				{
					gene[j].dic[0] = (double)di[6] / len_prom1 / dexp[6];
				}
				else
				{                   
					gene[j].dic[0] = (double)di[3] / len_prom1 / dexp[3];
					gene[j].dic[1] = (double)di[6] / len_prom1 / dexp[6];
					gene[j].dic[2] = (double)di[8] / len_prom1 / dexp[8];
					gene[j].dic[3] = (double)di[9] / len_prom1 / dexp[9];					
				}
			}
		}
		else break;
		j++;
	}
	fclose(in_genelist);
	//	printf("Up %d\tDown %d\tGood non-DEG %d\tBad non-DEGs %d Other %d\tTotal %d\n", ndegup, ndegdo, nnedeg_good, nnedeg_bad, other_bad, ndegup + ndegdo + nnedeg_good + nnedeg_bad + other_bad);
	int size_nedeg_min = (int)(nnedeg_good * size_nedeg_fraction);
	int ndeg = ndegup + ndegdo;
	int nseq_deg_nedeg = ndeg + nnedeg_good;
	// gene_type -1 0 1 2 3  = ngen_absent nnedeg_bad nnedeg_good ndegup ndegdo 	
	int nnedeg_good0 = ndeg;
	int ndegup0 = ndegdo;
	int ndegdo0 = 0;
	qsort(gene, n_genes_tot, sizeof(gene[0]), compare_genetype);
	/*	printf("All\n");
		for (i = 0; i < 10; i++)
		{
			int x = i + ndegdo0;
			printf("GType %d\tGinx %d\tpadj %g\tLog2FC %f\n", gene[x].gen_type, gene[x].num, gene[x].padj, gene[x].log2fc);
		}
	*/
	if (border == 1)
	{
		qsort(&gene[ndegup0], ndegup, sizeof(gene[0]), compare_fc_top);
		qsort(&gene[ndegdo0], ndegdo, sizeof(gene[0]), compare_fc_top);
	}
	else
	{
		qsort(&gene[ndegup0], ndegup, sizeof(gene[0]), compare_fc_bottom);
		qsort(&gene[ndegdo0], ndegdo, sizeof(gene[0]), compare_fc_bottom);
	}	
	/*printf("Down\n");
	for (i = 0; i < 10; i++)
	{
		int x = i + ndegdo0;
		printf("GType %d\tGinx %d\tpadj %g\tLog2FC %f\n", gene[x].gen_type, gene[x].num, gene[x].padj, gene[x].log2fc);
	}
	
	printf("Up\n");	
	for (i = -1; i <= ndegup; i++)
	{
		int x = i + ndegup0;
		printf("GType %d\tGinx %d\tpadj %g\tLog2FC %f\n", gene[x].gen_type, gene[x].num, gene[x].padj, gene[x].log2fc);
	}*/
	double* thr_log2fc_up1, * thr_log2fc_do1;
	thr_log2fc_up1 = new double[nseq_deg];
	if (thr_log2fc_up1 == NULL) { printf("Out of memory..."); return -1; };
	thr_log2fc_do1 = new double[nseq_deg];
	if (thr_log2fc_do1 == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < nseq_deg; i++)thr_log2fc_up1[i] = 1000;
	for (i = 0; i < nseq_deg; i++)thr_log2fc_do1[i] = 1000;

	int nseq_deg1 = nseq_deg - 1;
	int take_up = 0, take_do = 0, take_no = 0;
	int print_up = 1, print_do = 1;
	int iup_last = -1, ido_last = -1;
	int count_up = 0, count_do = 0;	
	for (i = ndegup0; i < ndegup0 + ndegup; i++)
	{
		if (gene[i].log2fc > 0)
		{
			if (print_up == 1)
			{
				take_up++;
				int rest = take_up % size_deg_seg;
				gene[i].type = 1 + count_up;
				if (rest == 0)
				{
					if (border != 1)count_up++;
					if (count_up <= nseq_deg1)thr_log2fc_up1[count_up] = gene[i].log2fc;
					if (border == 1)count_up++;
				}
				iup_last = i;
				if (take_up == size_deg_min)print_up = 0;
			}
		}
		if (print_up == 0)break;
	}
	for (i = ndegdo0; i < ndegdo0 + ndegdo; i++)
	{
		if (gene[i].log2fc < 0)
		{
			if (print_do == 1)
			{
				take_do++;
				int rest = take_do % size_deg_seg;
				gene[i].type = 1 + count_do;
				ido_last = i;
				if (rest == 0)
				{
					if (border != 1)count_do++;
					if (count_do <= nseq_deg1)thr_log2fc_do1[count_do] = gene[i].log2fc;
					if (border == 1)count_do++;
				}
				if (take_do == size_deg_min)print_do = 0;
			}
		}
		if (print_do == 0)break;
	}
	double thr_log2fc_up = 10, thr_log2fc_do = -10;
	if (take_up > 0)
	{
		thr_log2fc_up = gene[iup_last].log2fc;
	}
	if (take_do > 0)
	{
		thr_log2fc_do = gene[ido_last].log2fc;
	}
	take_no = nnedeg_good;
	/*
	double thr_log2fc_nemin = 0, thr_log2fc_nemax = 0;
	for (i = ndeg_nndeg - 1; i >= 0; i--)
	{
		double log2fc = rnaseq[i].log2fc;
		if (log2fc >= thr_log2fc_do && log2fc <= thr_log2fc_up)
		{
			if (log2fc > thr_log2fc_nemax)thr_log2fc_nemax = log2fc;
			if (log2fc < thr_log2fc_nemin)thr_log2fc_nemin = log2fc;
			take_no++;
			gene[rnaseq[i].num].deg_type = rnaseq[i].type = 0;
			if (take_no == size_nedeg_min)break;
		}
	}
	*/
	double** vec_up, ** vec_do, vec_no[SPACESIZE], vec_dif[SPACESIZE];
	vec_up = new double* [ndegup];
	if (vec_up == NULL) return -1;
	for (i = 0; i < ndegup; i++)
	{
		vec_up[i] = new double[space_size];
		if (vec_up[i] == NULL) return -1;
	}
	vec_do = new double* [ndegdo];
	if (vec_do == NULL) return -1;
	for (i = 0; i < ndegdo; i++)
	{
		vec_do[i] = new double[space_size];
		if (vec_do[i] == NULL) return -1;
	}
	int mah_up = 0, mah_do = 0;
	for (i = 0; i < ndeg + nnedeg_good; i++)
	{
		int ii = gene[i].gen_type;
		if (ii <= 1)continue;
		else
		{
			if (ii == 2 && mah_up < ndegup)
			{
				for (k = 0; k < space_size1; k++)vec_up[mah_up][k] = gene[i].dic[k];
				vec_up[mah_up][space_size1] = gene[i].moc;
				mah_up++;
			}
			if (ii == 3 && mah_do < ndegdo)
			{
				for (k = 0; k < space_size1; k++)vec_do[mah_do][k] = gene[i].dic[k];
				vec_do[mah_do][space_size1] = gene[i].moc;
				mah_do++;
			}
		}
	}
	if ((outlog = fopen(file_log, "wt")) == NULL)
	{
		fprintf(outlog, "Input file %s can't be opened!\n", file_log);
		exit(1);
	}
	double av_up[SPACESIZE], av_do[SPACESIZE];
	double uw_up[SPACESIZE][SPACESIZE], uw_do[SPACESIZE][SPACESIZE];
	for (k = 0; k < space_size; k++)av_up[k] = av_do[k] = 0;	
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)uw_up[k][i]  = uw_do[k][i] = 0;
	}
	//up
	for (j = 0; j < mah_up; j++)
	{
		for (k = 0; k < space_size; k++)
		{
			for (i = 0; i < space_size; i++)uw_up[k][i] += vec_up[j][k] * vec_up[j][i];
		}
		for (k = 0; k < space_size; k++)av_up[k] += vec_up[j][k];
	}
	for (k = 0; k < space_size; k++)av_up[k] /= mah_up;
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)
		{
			uw_up[k][i] -= space_size * av_up[k] * av_up[i];
		}
	}
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)
		{
			uw_up[k][i] /= space_size1;
		}
	}
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)uw[k][i] = uw_up[k][i];
	}
	if (BackMat(space_size, outlog) == -1)
	{
		printf("BackMat Error\n");
		back_err++;
		return -1;
	}
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)uw_up[k][i] = uw[k][i];
	}
	//do
	for (j = 0; j < mah_do; j++)
	{
		for (k = 0; k < space_size; k++)
		{
			for (i = 0; i < space_size; i++)uw_do[k][i] += vec_do[j][k] * vec_do[j][i];
		}
		for (k = 0; k < space_size; k++)av_do[k] += vec_do[j][k];
	}
	for (k = 0; k < space_size; k++)av_do[k] /= mah_do;
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)
		{
			uw_do[k][i] -= space_size * av_do[k] * av_do[i];
		}
	}
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)
		{
			uw_do[k][i] /= space_size1;
		}
	}
	//up
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)uw[k][i] = uw_do[k][i];
	}
	if (BackMat(space_size, outlog) == -1)
	{
		printf("BackMat Error\n");
		back_err++;
		return -1;
	}
	for (k = 0; k < space_size; k++)
	{
		for (i = 0; i < space_size; i++)uw_do[k][i] = uw[k][i];
	}

	//printf("nonDEGs\n");
	double thr_log2fc_nemin = 0, thr_log2fc_nemax = 0;
	for (i = nnedeg_good0; i < nnedeg_good0 + nnedeg_good; i++)
	{
		double log2fc = gene[i].log2fc;
		if (log2fc >= thr_log2fc_do && log2fc <= thr_log2fc_up)
		{
			if (log2fc > thr_log2fc_nemax)thr_log2fc_nemax = log2fc;
			if (log2fc < thr_log2fc_nemin)thr_log2fc_nemin = log2fc;
		}
	}
	double avr_up = 0, std_up = 0, avr_do = 0, std_do = 0, thr_md_up, thr_md_do;
	//for (i = nnedeg_good0; i < nnedeg_good0 + nnedeg_good; i++)
	for (i = 0; i < n_genes_tot; i++)
	{
		gene[i].rm_up = 0;
		for (k = 0; k < space_size1; k++)vec_no[k] = gene[i].dic[k];
		vec_no[space_size1] = gene[i].moc;
		for (k = 0; k < space_size; k++)vec_dif[k] = vec_no[k] - av_up[k];		
		for (k = 0; k < space_size; k++)
		{
			double buf = 0;
			for (j = 0; j < space_size; j++)buf += uw_up[k][j] * vec_dif[j];
			gene[i].rm_up += buf * vec_dif[k];
		}
		gene[i].rm_up = sqrt(gene[i].rm_up);
		gene[i].rm_do = 0;
		for (k = 0; k < space_size; k++)vec_dif[k] = vec_no[k] - av_do[k];
		for (k = 0; k < space_size; k++)
		{
			double buf = 0;
			for (j = 0; j < space_size; j++)buf += uw_do[k][j] * vec_dif[j];
			gene[i].rm_do += buf * vec_dif[k];
		}
		gene[i].rm_do  = sqrt(gene[i].rm_do);
		gene[i].rm2 = gene[i].rm_up + gene[i].rm_do;
	//	printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%g\t%f\n", i + 1, gene[i].gen_type,gene[i].num, gene[i].type, gene[i].rm_up, gene[i].rm_do, gene[i].rm2, gene[i].padj, gene[i].log2fc);
	}
	for (i = ndegup0; i < ndegup0 + ndegup; i++)
	{
		avr_up += gene[i].rm_up;
		std_up += gene[i].rm_up * gene[i].rm_up;
	}
	avr_up /= ndegup;
	std_up = (std_up - ndegup * avr_up * avr_up) / (ndegup - 1);
	std_up = sqrt(std_up);
	for (i = ndegdo0; i < ndegdo0 + ndegdo; i++)	
	{
		avr_do += gene[i].rm_do;
		std_do += gene[i].rm_do * gene[i].rm_do;
	}
	avr_do /= ndegdo;
	std_do = (std_do - ndegdo * avr_do * avr_do) / (ndegdo - 1);	
	std_do = sqrt(std_do);
	thr_md_up = avr_up + z_sco_thr * std_up;
	thr_md_do = avr_do + z_sco_thr * std_do;
	qsort(&gene[nnedeg_good0], nnedeg_good, sizeof(gene[0]), compare_mah);
	/*printf("downDEGs\n");
	for (i = 0; i < ndegdo; i++)
	{
		printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%g\t%f\n", i + 1, gene[i].gen_type, gene[i].num, gene[i].type, gene[i].rm_up, gene[i].rm_do, gene[i].rm2, gene[i].padj, gene[i].log2fc);
	}
	printf("upDEGs\n");
	for (i = ndegup0; i < ndeg; i++)
	{
		printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%g\t%f\n", i + 1, gene[i].gen_type, gene[i].num, gene[i].type, gene[i].rm_up, gene[i].rm_do, gene[i].rm2, gene[i].padj, gene[i].log2fc);
	}
	printf("non-DEGs\n");*/
	/*for (i = nnedeg_good0; i < nnedeg_good0 + nnedeg_good; i++)
	{
		printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%g\t%f\n", i + 1, gene[i].gen_type, gene[i].num, gene[i].type, gene[i].rm_up, gene[i].rm_do, gene[i].rm2, gene[i].padj, gene[i].log2fc);
	}*/
	/*printf("NonDEGs\n");
	for (i = 0; i < 20; i++)
	{
		int x = i + nnedeg_good0;
		printf("GType %d\tGinx %d\tpadj %g\tLog2FC %f\tR %f\n", gene[x].gen_type, gene[x].num, gene[x].padj, gene[x].log2fc, gene[x].rm2);
	}*/

	int nseq_nondeg_mdc = 0;//count of nonDEG with adopted mono- &di- content
	{
		int nb_max = Max(ndegup, ndegdo);
		nb_max = (int)(rat_nf_nb_max * nb_max);		
		for (i = nnedeg_good0; i < nnedeg_good0 + nnedeg_good; i++)
		{
			//int ii = gene[i].gen_type;
			//if (ii == 1)//non-DEG without mono- & di- content
			{
				if (gene[i].rm_up < thr_md_up && gene[i].rm_do < thr_md_do)
				{
					gene[i].sim = 1;
					nseq_nondeg_mdc++;
					if (nseq_nondeg_mdc >= nb_max)break;
				}
				//printf("%d\t%d\t%d\t%f\t%f\t%f\t%g\t%f\n", i + 1, gene[i].num, gene[i].type, gene[i].rm_up, gene[i].rm_do, gene[i].rm2, gene[i].padj, gene[i].log2fc);			
			}
		}
	}
	qsort(gene, n_genes_tot, sizeof(gene[0]), compare_numg);
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
			
	/*if (random == 1)
	{
		int nseq_mix = ndegup + ndegdo + nseq_nondeg_mdc;
		int* indices;
		indices = new int[nseq_mix];
		if (indices == NULL) { puts("Out of memory..."); exit(1); }
		srand((unsigned)time(NULL));			
		for (i = 0; i < nseq_mix; i++)
		{
			indices[i] = -1;
		}
		k = 0;
		for (i = 0; i < n_genes_tot; i++)
		{
			if (gene[i].gen_type >= 2 || (gene[i].gen_type == 1 && gene[i].sim == 1))indices[k++] = i;
		}
		int nseq_mix1 = nseq_mix - 1;
		for (i = 0; i < nseq_mix1; i++)
		{
			int i1 = indices[i];																		
			int i2 = indices[1 + i + (rand() % (nseq_mix1 - i))];													
			Mix(&gene[i1].gen_type, &gene[i2].gen_type);
			Mix(&gene[i1].sim, &gene[i2].sim);																		
		}
		delete[] indices;
	}*/
	int countcas = 0;
	int lost[5] = { 0,0,0,0,0 };
	for (i = 0; i < n_genes_tot; i++)
	{
	//	printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%g\t%f\t%d\n", gene[i].numg, gene[i].num, gene[i].type, gene[i].gen_type, gene[i].rm_up, gene[i].rm_do, gene[i].rm2, gene[i].padj, gene[i].log2fc, gene[i].sim);
		int ii = gene[i].gen_type;
		int cell = gene[i].type - 1;
		switch (ii)
		{
		case 0:
		{
			for (j = 0; j < nseq_deg; j++)
			{
				fprintf(out_up[j], "0\n");
				fprintf(out_do[j], "0\n");
				fprintf(out_no[j], "0\n");
			}
			lost[ii]++;
			countcas++;
			break;
		}
		case 1:
		{
			if (gene[i].sim == 1)
			{
				for (j = 0; j < nseq_deg; j++)
				{
					fprintf(out_up[j], "0\n");
					fprintf(out_do[j], "0\n");
					fprintf(out_no[j], "1\n");
				}
			}
			else
			{
				for (j = 0; j < nseq_deg; j++)
				{
					fprintf(out_up[j], "0\n");
					fprintf(out_do[j], "0\n");
					fprintf(out_no[j], "0\n");
				}
			}
			lost[ii]++;
			countcas++;
			break;
		}
		case 2:
		{			
			for (j = 0; j < nseq_deg; j++)
			{
				if (j == cell)fprintf(out_up[j], "1\n");
				else fprintf(out_up[j], "0\n");
				fprintf(out_do[j], "0\n");
				fprintf(out_no[j], "0\n");
			}
			lost[ii]++;
			countcas++;
			break;
		}
		case 3:
		{
			
			for (j = 0; j < nseq_deg; j++)
			{
				if (j == cell)fprintf(out_do[j], "1\n");
				else fprintf(out_do[j], "0\n");
				fprintf(out_up[j], "0\n");
				fprintf(out_no[j], "0\n");
			}
			lost[ii]++;
			countcas++;
			break;
		}
		default:
		{
			for (j = 0; j < nseq_deg; j++)
			{
				fprintf(out_up[j], "0\n");
				fprintf(out_do[j], "0\n");
				fprintf(out_no[j], "0\n");
			}
			lost[4]++;
			countcas++;
			break;
		}
		}				
	}
	for (j = 0; j < nseq_deg; j++)
	{
		fclose(out_up[j]);
		fclose(out_do[j]);
		fclose(out_no[j]);
	}
	for (j = 0; j < nseq_deg; j++)
	{
		if (thr_log2fc_up1[j] != 1000)thr_log2fc_up1[j] = pow((double)2, thr_log2fc_up1[j]);
		if (thr_log2fc_do1[j] != 1000)thr_log2fc_do1[j] = pow((double)2, thr_log2fc_do1[j]);
	}
	if ((out_sta = fopen(fileo_sta, "at")) == NULL)
	{
		printf("Input file %s can't be opened!\n", fileo_sta);
		return -1;
	}
	if (rev == 0)
	{
		fprintf(out_sta, "%s\t%s\t%s\t%s\t%s\t", filei_rnaseq, filei_genelist, argv[7], argv[8], argv[9]);
		fprintf(out_sta, "Top genes\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t", n_genes_tot, argv[10], ndegup, take_up, argv[11], ndegdo, take_do, argv[12], nnedeg_good, nnedeg_bad, nseq_nondeg_mdc);
		fprintf(out_sta, "UpFC\t%f\tDoFC\t%f\tNoFC\t%f\t%f\t", thr_log2fc_up, thr_log2fc_do, thr_log2fc_nemax, thr_log2fc_nemin);
		fprintf(out_sta, "UPs");
		for (j = 0; j < nseq_deg; j++)fprintf(out_sta, "\t%f", thr_log2fc_up1[j]);
		fprintf(out_sta, "\tDOs");
		for (j = 0; j < nseq_deg; j++)fprintf(out_sta, "\t%f", thr_log2fc_do1[j]);
		fprintf(out_sta, "\tNdi\t%d\tRand\t%d", ndity, random);
		fprintf(out_sta, "\n");		
	}
	else
	{
		fprintf(out_sta, "%s\t%s\t%s\t%s\t%s\t", filei_rnaseq, filei_genelist, argv[7], argv[8], argv[9]);
		fprintf(out_sta, "Bottom genes\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t", n_genes_tot, argv[10], ndegup, take_up, argv[11], ndegdo, take_do, argv[12], nnedeg_good, nnedeg_bad, nseq_nondeg_mdc);
		fprintf(out_sta, "UpFC\t%f\tDoFC\t%f\tNoFC\t%f\t%f\t", thr_log2fc_up, thr_log2fc_do, thr_log2fc_nemax, thr_log2fc_nemin);
		fprintf(out_sta, "UPs");
		for (j = 0; j < nseq_deg; j++)fprintf(out_sta, "\t%f", thr_log2fc_do1[j]);
		fprintf(out_sta, "\tDOs");
		for (j = 0; j < nseq_deg; j++)fprintf(out_sta, "\t%f", thr_log2fc_up1[j]);
		fprintf(out_sta, "\tNdi\t%d\tRand\t%d", ndity, random);
		fprintf(out_sta, "\n");
	}
	fclose(out_sta);
	fclose(outlog);
	if (back_err == 0)remove(file_log);
	for (i = 0; i < ngene_rnaseq; i++)
	{
		delete[] genes[i];
	}
	delete[] genes;
	delete[] out_up;
	delete[] out_do;
	delete[] thr_log2fc_up1;
	delete[] thr_log2fc_do1;
	for (i = 0; i < ndegup; i++)
	{
		delete[] vec_up[i];
	}
	delete[] vec_up;
	for (i = 0; i < ndegdo; i++)
	{
		delete[] vec_do[i];
	}
	delete[] vec_do;
	return 1;
}


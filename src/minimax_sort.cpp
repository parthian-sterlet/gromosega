#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
// ranging TF families from Gromosega prediction
 
struct rang //for selected 8 points
{
	double sum;
	int cou;
	int inx;	
	char name[100];
};
int compare_rang(const void* X1, const void* X2)
{
	struct rang* S1 = (struct rang*)X1;
	struct rang* S2 = (struct rang*)X2;
	if (S1->cou - S2->cou > 0)return -1;
	if (S1->cou - S2->cou < 0)return 1;
	if (S1->sum - S2->sum > 0)return -1;
	if (S1->sum - S2->sum < 0)return 1;
	return 0;
}
int compare_sum(const void* X1, const void* X2)
{
	struct rang* S1 = (struct rang*)X1;
	struct rang* S2 = (struct rang*)X2;
	if (S1->sum - S2->sum > 0)return -1;
	if (S1->sum - S2->sum < 0)return 1;
	if (S1->cou - S2->cou > 0)return -1;
	if (S1->cou - S2->cou < 0)return 1;
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
	return -1;
}
char* UnderStolStr(char* str, char* ret, size_t size, int nstol, char razd)
{
	memset(ret, 0, size);
	int p1, p2;
	if (nstol == 0)p1 = -1;
	else p1 = StrNStr(str, razd, nstol);
	p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -1)
	{
		p2 = (int)strlen(str);
	}
	if (p1 == -1 || p2 == -1) return NULL;
	int len = p2 - p1 - 1;
	strncpy(ret, &str[p1 + 1], len);
	ret[len] = '\0';
	return ret;
}
double UnderStolDouble(char* str, int nstol, char* ret, size_t size, const char razd)
{
	memset(ret, 0, size);
	if (nstol == 0)return atof(str);
	int p1 = StrNStr(str, razd, nstol);
	int p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -10000)
	{
		p2 = (int)strlen(str);
	}
	if (p1 == -1 || p2 == -1) return -1;
	int len = p2 - p1 - 1;
	strncpy(ret, &str[p1 + 1], len);
	ret[len] = '\0';
	int cd = (int)ret[0];
	if (isdigit(cd) || cd == 45)return atof(ret);//0123456789 or -
	else return -1;
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
int GetTerm(char* d, char* id, int const pos_sta1, int const pos_end1, int const pos_sta2, int const pos_end2, const char end1, const char end2, const char end3)
{
	int i;
	int ista = 0;
	int j = 0, iend = (int)strlen(d);	
	for (i = ista;i < iend; i++)
	{		
		int c = int(d[i]);
		if ((c >= pos_sta1 && c <= pos_end1) || (c >= pos_sta2 && c <= pos_end2))
		{
			id[j++] = d[i];
			continue;
		}
		if (d[i] == end1 || (d[i] == end2 || d[i] == end3))
		{
			if (j > 0)break;
		}
	}
	id[j] = '\0';
	return j;
}
int main(int argc, char* argv[])
{
	int i, j, k;
	char d[5000], s[100], track_name[2][300], track_name0[300], path[300], filei_name[300], filei[300], fileo[300], fileo_base[300];
	FILE* out, * in;
	
	if (argc != 6)
	{
		printf("%s 1path_input 2file_input 3file_base output 4int n_parts (1/2 = one/two part(s) 5int ntop families\n", argv[0]);// f f 9 2 10 5 5
		return -1;
	}
	strcpy(path, argv[1]);//in_file
	strcpy(filei_name, argv[2]);//in_file
	strcpy(fileo_base, argv[3]);//out_file	
	int npar = atoi(argv[4]);//1 2
	int ntop = atoi(argv[5]);//10 default
	
	memset(filei, '\0', sizeof(filei));
	strcpy(filei, path);
	strcat(filei, filei_name);
	if ((in = fopen(filei, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", filei);
		return -1;
	}
	const char tab = '\t', end1 = '_', end2 = ' ', end3 = '.';
	int const pos_sta1 = 65, pos_end1 = 90, pos_sta2 = 97, pos_end2 = 122;// A, Z, a, z
	char buf[50];
	char types[2][20];
	strcpy(types[0], "sole");
	strcpy(types[1], "slice");
	fgets(d, sizeof(d), in);//header with counts of TF family TFBS motifs
	fgets(d, sizeof(d), in);//header with TF family names
	DelChar(d, '\n');
	int nfamilies = 0; // total number of TF families
	int dlen = (int)strlen(d);
	for (i = 0; i < dlen; i++)
	{
		if (d[i] == tab)nfamilies++;
	}
	nfamilies--;
	if (ntop > nfamilies)ntop = nfamilies;

	char** family_names;
	const size_t lens = 100;
	size_t sizemot = lens * sizeof(track_name[0][0]);
	family_names = new char* [nfamilies];
	if (family_names == NULL) { printf("Out of memory..."); return -1; };
	for (i = 0; i < nfamilies; i++)
	{
		family_names[i] = new char[lens];
		if (family_names[i] == NULL) { puts("Out of memory..."); exit(1); }
	}	
	for (i = 0; i < nfamilies; i++)
	{
		memset(family_names[i], '\0', sizeof(sizemot));
		if (UnderStolStr(d, family_names[i], sizemot, i+2, tab) == NULL) 
		{ 
			printf("Wrong format %s\n", filei); 
			return(-1); 
		}
	}
	rang* family_ratio;
	family_ratio = new rang[nfamilies];
	if (family_ratio == NULL) { puts("Out of memory..."); exit(1); }

	int nstr[2] = { 0,0 };
	k = 0;
	i = 0;
	while (fgets(d, sizeof(d), in) != NULL)
	{		
		memset(track_name0, '\0', sizeof(sizemot));
		int test = GetTerm(d, track_name0, pos_sta1, pos_end1, pos_sta2, pos_end2, end1, end2, end3);
		if (test == 0)
		{
			printf("File %s reading error!\n", filei);
			exit(1);
		}
		if (k == 0 && i == 0)strcpy(track_name[k], track_name0);
		else
		{
			if (strcmp(track_name[k], track_name0) != 0)
			{
				k++;
				strcpy(track_name[k], track_name0);
				i = 0;
			}
		}
		i++;
		double val = UnderStolDouble(d, 1, buf, sizeof(buf), tab);
		if (val >= 1)nstr[k]++;				
	}	
	rewind(in);
	fgets(d, sizeof(d), in);//header with counts of TF family TFBS motifs
	fgets(d, sizeof(d), in);//header with TF family names
	DelChar(d, '\n');
	for (k = 0; k < npar; k++)
	{
		memset(fileo, '\0', sizeof(fileo));
		strcpy(fileo, fileo_base); 
		strcat(fileo, "_");
		if(npar == 1)strcat(fileo, types[0]);
		else strcat(fileo, types[1]);
		strcat(fileo, "_");
		strcat(fileo, track_name[k]);
		if ((out = fopen(fileo, "at")) == NULL)
		{
			printf("Output file %s can't be opened!\n", fileo);
			return -1;
		}
		for (i = 0; i < nfamilies; i++)
		{
			family_ratio[i].cou = 0;
			family_ratio[i].sum = 0;
			family_ratio[i].inx = i;
			memset(family_ratio[i].name, '\0', sizeof(sizemot));
			strcpy(family_ratio[i].name, family_names[i]);
		}
		for (j = 0; j < nstr[k]; j++)
		{
			fgets(d, sizeof(d), in);
			DelChar(d, '\n');			
			for (i = 0; i < nfamilies; i++)
			{
				memset(s, '\0', sizeof(s));
				double val = UnderStolDouble(d, i + 2, buf, sizeof(buf), tab);
				if (val > 0)
				{
					family_ratio[i].cou++;
					family_ratio[i].sum += val;
				}
			}
		}
		qsort(family_ratio, nfamilies, sizeof(family_ratio[0]), compare_sum);
		fprintf(out, "%s", path);
		fprintf(out, "\t%s", filei_name);
		fprintf(out, "\t%s", track_name[k]);
		for (i = 0; i < ntop; i++)fprintf(out, "\t%s", family_ratio[i].name);
		fprintf(out, "\t\tFraction sets, %%");
		for (i = 0; i < ntop; i++)fprintf(out, "\t%.3f", 100*(double)family_ratio[i].cou / nstr[k]);
		fprintf(out, "\t\tSum ratios, *1000");
		for (i = 0; i < ntop; i++)fprintf(out, "\t%.3f", 1000*family_ratio[i].sum / nstr[k]);
		fprintf(out, "\n");
		fclose(out);
	}	
	fclose(in);	
	delete[] family_ratio;
	for (k = 0; k < nfamilies; k++)
	{
		delete[] family_names[k];		
	}	
	delete[] family_names;
	return 1;
}

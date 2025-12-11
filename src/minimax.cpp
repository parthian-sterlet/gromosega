#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  <time.h>
#include <ctype.h>

#define Min(a,b) ((a)>(b))? (b):(a);
#define Max(a,b) ((a)>(b))? (a):(b);
#define POOL 200 // size of population = 4islo osobey populyacii
#define ELIT 100 // size of population = 4islo osobey elity populyacii
#define MMAX 1700  //max total number of motifs

void MixI(int* a, int* b)
{
	int buf = *a;
	*a = *b;
	*b = buf;
}
void BigMixI(int* d1, int len) // pereme6ivanie stroki
{
	int r;
	for (r = 0; r < len - 1; r++)
	{
		MixI(&d1[r], &d1[1 + r + (rand() % (len - 1 - r))]);
	}
}
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
int StrNStr(char* str, char c, int n)
{
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
double UnderStol(char* str, int nstol, char razd)
{
	if (nstol == 0)return atof(str);
	char ret[100];
	memset(ret, 0, sizeof(ret));
	int p1 = StrNStr(str, razd, nstol);
	int p2 = StrNStr(str, razd, nstol + 1);
	if (p2 == -1)
	{
		p2 = (int)strlen(str);
	}
	if (p1 == -1 || p2 == -1) return -1;
	int len = p2 - p1 - 1;
	strncpy(ret, &str[p1 + 1], len);
	ret[len] = '\0';
	return atof(ret);
}

struct combi {
	int* mot;// presence of motifs in individual
	double fit;
	double xp;
	double xn;
	double x2p;
	double x2n;
	double prec[20]; // Precision at TP rates 0.05 0.1.. 0.9 0.95 1.0
	void get_copy(combi* a, int mtot);
	void init0(int mot);
	void init1(int msel, int mot);
	void term(void);
	void fprint_all(int msel, FILE *out);
	void fprint_nam(int mtot, char **nam_mot, char** nam_class, FILE *out);
} *pop, det1, det2[2];//population
int compare_pop(const void* X1, const void* X2)
{
	struct combi* S1 = (struct combi*)X1;
	struct combi* S2 = (struct combi*)X2;
	if (S1->fit - S2->fit > 0)return -1;
	if (S1->fit - S2->fit < 0)return 1;
	return 0;
}
void combi::init0(int mtot)
{	
	mot = new int[mtot];
	if (mot == NULL)
	{
		puts("Out of memory...");
		exit(1);
	}
}
void combi::init1(int msel, int mtot)
{
	fit = -1;
	int j, ord[MMAX];
	for (j = 0; j < msel; j++)ord[j] = 1;
	for (j = msel; j < mtot; j++)ord[j] = 0;
	BigMixI(ord, mtot);
	for (j = 0; j < mtot; j++)mot[j] = ord[j];
}
void combi::term(void)
{
	delete[] mot;
}
void combi::get_copy(combi* a, int mtot)
{
	int i;
	a->fit = fit;
	for (i = 0; i < mtot; i++)a->mot[i] = mot[i];	
}
void combi::fprint_all(int mtot, FILE* out)
{
	int i;
	fprintf(out,"Fit %f Motifs#: ", fit);
	for (i = 0; i < mtot; i++)fprintf(out,"%2d", mot[i]);	
	fprintf(out,"\n");
}
void combi::fprint_nam(int mtot, char** nam_mot, char** nam_class, FILE* out)
{
	int i;
	fprintf(out,"%f\tMotifs\t", fit);
	for (i = 0; i < mtot; i++)if (mot[i] == 1)fprintf(out,"%s\t", nam_mot[i]);
	fprintf(out, "ClassesFamilies\t");
	for (i = 0; i < mtot; i++)if (mot[i] == 1)fprintf(out, "%s\t", nam_class[i]);
	for (i = 0; i < 20; i++)fprintf(out,"\t%.3f", prec[i]);
	fprintf(out,"\n");
}
struct qbs {
	double err;//ERR score
	int pn;// 1 pos 0 neg
};
int compare_qq(const void* X1, const void* X2)//increase
{
	double X = (*(double*)X1 - *(double*)X2);
	if (X > 0)return -1;
	if (X < 0)return 1;
	return 0;
}
int compare_qbs(const void* X1, const void* X2)//decrease
{
	struct qbs* S1 = (struct qbs*)X1;
	struct qbs* S2 = (struct qbs*)X2;
	if (S1->err - S2->err > 0)return -1;
	if (S1->err - S2->err < 0)return 1;
	return 0;
}
int MutCo(combi* a, int msel, int mtot)
{
	int i, r1, r2, mtot1 = mtot - 1, trys = mtot * 5, gom = 0;
	for(i=0;i<trys;i++)
	{
		r1 = rand() % mtot;
		r2 = rand() % mtot1;
		if (r2 >= r1)r2++;
		if (a->mot[r1] != a->mot[r2])
		{
			gom = 1;
			break;
		}
	}
	if(gom == 1)
	{
		int buf = a->mot[r1];
		a->mot[r1] = a->mot[r2];
		a->mot[r2] = buf;
		return 1;
	}
	puts("Mutation error...");
	return -1;
}
int RecCo(combi* a, combi* b, int msel, int mnsel, int mtot)
{
	int ac[MMAX], bc[MMAX], ii[MMAX], i, ra = -1, rb = -1, gom = 0;
	for (i = 0; i < mtot; i++)
	{
		ac[i] = a->mot[i];
		bc[i] = b->mot[i];
		ii[i] = i;
	}
	BigMixI(ii, mtot);
	for (i = 0; i < mtot; i++)
	{
		int ix = ii[i];
		if (a->mot[ix] != b->mot[ix])
		{
			if (a->mot[ix] > b->mot[ix] && ra == -1)
			{
				gom++;
				ra = ix;
				if (gom == 2)break;
			}
			if (a->mot[ix] < b->mot[ix] && rb == -1)
			{
				gom++;
				rb = ix;
				if (gom == 2)break;
			}
		}
	}
	if (gom == 2)
	{
		int buf;
		buf = a->mot[ra];
		a->mot[ra] = b->mot[ra];
		b->mot[ra] = buf;
		buf = a->mot[rb];
		a->mot[rb] = b->mot[rb];
		b->mot[rb] = buf;
		return 1;
	}
	puts("Recombination error...");
	return -1;
}
int GomTown(combi a, combi b, int mtot)
{
	int i;
	for (i = 0; i < mtot; i++)
	{
		if (b.mot[i] != a.mot[i])return 0;
	}
	return 1;
}
int EvalFit(combi* a, int msel, int mtot, double** errp, double** errn, int nseqp, int nseqn, int nseqtot, double nseqrat, double prec_exp, double fp2, qbs *seqtot)
{
	int i, j, k;
	int mh[MMAX];
	for (i = 0; i < MMAX; i++)mh[i] = 0;
	k = 0;
	for (i = 0; i < mtot; i++)
	{
		if (a->mot[i] == 1)mh[k++] = i;
	}	
	k = 0;
	for (j = 0; j < nseqp; j++)
	{
		double max = 0;
		for (i = 0; i < msel; i++)
		{
			int ii = mh[i];
			if (errp[ii][j] > max)max = errp[ii][j];						
		}
		seqtot[k].pn = 1;
		seqtot[k++].err = max;
	}
	for (j = 0; j < nseqn; j++)
	{
		double max = 0;
		for (i = 0; i < msel; i++)
		{
			int ii = mh[i];
			if (errn[ii][j] > max)max = errn[ii][j];			
		}
		seqtot[k].pn = 0;
		seqtot[k++].err = max;
	}
	qsort(seqtot, nseqtot, sizeof(seqtot[0]), compare_qbs);
	int tp = 0, fp = 0, dtp = 0, dfp = 0;
	int nseqtot1 = nseqtot - 1;
	double prec_pred = 1, auc_pr = 0, tpc_pred = 0, fpc_pred = 0, err_pred =0, mnoj = (1 - prec_exp) / nseqp;
	double fp2pow = pow((double)10, -fp2);
	for (i = 0; i < nseqtot; i++)
	{
		if (seqtot[i].pn == 1)dtp++;
		else dfp++;
		int i1 = i + 1;
		if ((i == nseqtot1 || (seqtot[i].err != seqtot[i1].err)) && dtp > 0) //seqtot[i1].pn == 0 && 
		{
			double tpc, fpc;
			if (seqtot[i].err < fp2)
			{
				double err_cur = pow((double)10, -seqtot[i].err);
				double wei = (fp2pow - err_pred) / (err_cur - err_pred);
				tpc = tpc_pred + wei * dtp;
				fpc = fpc_pred + wei * dfp;
			}
			else
			{
				tpc = tpc_pred + dtp;
				fpc = fpc_pred + dfp;
			}
			double prec_cur = tpc / (tpc + nseqrat * fpc);
			double prec_av = (prec_pred + prec_cur) / 2;
			double dauc = (prec_av - prec_exp) * (tpc - tpc_pred) * mnoj;
			auc_pr += dauc;
			tpc_pred = tpc;
			fpc_pred = fpc;
			prec_pred = prec_cur;
			err_pred = pow((double)10, -seqtot[i].err);
			if (seqtot[i].err < fp2)break;
			dtp = dfp = 0;
		}
	}
	a->fit = auc_pr;
	return 1;
}
int Precision(combi* a, int msel, int mtot, double** errp, double** errn, int nseqp, int nseqn, int nseqtot, double nseqrat, double prec_exp, double fp2, qbs* seqtot, 
	FILE *outh, char **motnames, int rank)
{
	int i, j, k;
	int mh[MMAX];
	for (i = 0; i < MMAX; i++)mh[i] = 0;
	k = 0;
	for (i = 0; i < mtot; i++)
	{
		if (a->mot[i] == 1)mh[k++] = i;
	}
	k = 0;
	for (j = 0; j < nseqp; j++)
	{
		double max = 0;
		for (i = 0; i < msel; i++)
		{
			int ii = mh[i];
			if (errp[ii][j] > max)max = errp[ii][j];
		}
		seqtot[k].pn = 1;
		seqtot[k++].err = max;
	}
	for (j = 0; j < nseqn; j++)
	{
		double max = 0;
		for (i = 0; i < msel; i++)
		{
			int ii = mh[i];
			if (errn[ii][j] > max)max = errn[ii][j];
		}
		seqtot[k].pn = 0;
		seqtot[k++].err = max;
	}
	qsort(seqtot, nseqtot, sizeof(seqtot[0]), compare_qbs);
	for (j = 0; j < 10; j++)a->prec[j] = prec_exp;
	fprintf(outh, "%d", rank);
	for (i = 0; i < msel; i++)
	{
		fprintf(outh, "\t%s", motnames[mh[i]]);		
	}
	fprintf(outh, "\nRecall\n\tPrecision\tAUC\tmLog10(ERR)\n");//\tTP\tFP\tTPC\tFPC
	double tpr20[20];
	tpr20[0] = 0.05;
	for (i = 1; i < 20; i++)tpr20[i] = tpr20[i - 1] + 0.05;
	int tp = 0, fp = 0, dtp = 0, dfp = 0;
	int nseqtot1 = nseqtot - 1;
	double prec_pred = 1, tpc_pred = 0, fpc_pred = 0, err_pred = 0, auc_pr = 0, mnoj = (1 - prec_exp)/nseqp;
	double fp2pow = pow((double)10, -fp2);
	//int count = 0;
	for (i = 0; i < nseqtot; i++)
	{
		if (seqtot[i].pn == 1)dtp++;
		else dfp++;
		int i1 = i + 1;
		if ((i == nseqtot1 || (seqtot[i].err != seqtot[i1].err)) && dtp > 0)//seqtot[i1].pn == 0 && 
		{
			double tpc, fpc;
			if (seqtot[i].err < fp2)
			{
				double err_cur = pow((double)10, -seqtot[i].err);
				double wei = (fp2pow - err_pred) / (err_cur - err_pred);
				tpc = tpc_pred + wei * dtp;
				fpc = fpc_pred + wei * dfp;
			}			
			else
			{
				tpc = tpc_pred + dtp;
				fpc = fpc_pred + dfp;
			}
			double prec_cur = tpc / (tpc + nseqrat * fpc);
			double prec_av = (prec_pred + prec_cur) / 2;
			double tpr = (tpc_pred + (tpc - tpc_pred)/2) / nseqp;	
			double dauc = (prec_av - prec_exp) * (tpc - tpc_pred) * mnoj;
			auc_pr += dauc;
			tp += dtp;
			fp += dfp;
			fprintf(outh, "%f\t%f\t%f\t%.18f\n", tpr, prec_av, auc_pr, seqtot[i].err);// ,tp,fp,tpc,fpc
			/*if (count % 50 == 0)
			{
				int yy = 0;
			}
			count++;*/
			for (j = 0; j < 20; j++)
			{
				if (tpr <= tpr20[j])a->prec[j] = prec_av;				
			}
			tpc_pred = tpc;
			fpc_pred = fpc;
			prec_pred = prec_cur;
			err_pred = pow((double)10, -seqtot[i].err);
			if (seqtot[i].err < fp2)break;
			dtp = dfp = 0;
		}
	}
	return 1;
}
int Modules(combi* a, int msel, int mtot, double** errp, double** errn, int nseqp, int nseqn, double fp2)
{
	int i, j, k;
	int mh[MMAX];
	for (i = 0; i < MMAX; i++)mh[i] = 0;
	k = 0;
	for (i = 0; i < mtot; i++)
	{
		if (a->mot[i] == 1)mh[k++] = i;
	}
	double xp = 0, xn = 0, x2p = 0, x2n = 0;
	k = 0;
	for (j = 0; j < nseqp; j++)
	{
		double max = fp2;
		for (i = 0; i < msel; i++)
		{
			int ii = mh[i];
			if (errp[ii][j] > max)max = errp[ii][j];
		}
		double dmax = max - fp2;
		xp += dmax;
		x2p += dmax * dmax;
	}
	for (j = 0; j < nseqn; j++)
	{
		double max = fp2;
		for (i = 0; i < msel; i++)
		{
			int ii = mh[i];
			if (errn[ii][j] > max)max = errn[ii][j];
		}
		double dmax = max - fp2;
		xn += dmax;
		x2n += dmax * dmax;
	}
	a->xp = xp;
	a->xn = xn;
	a->x2p = x2p;
	a->x2n = x2n;
	return 1;
}
int PearsonInternal(combi* a, int msel, int mtot, double** errp, double** errn, int nseqp, int nseqn, char** nam_mot, FILE *out)
{
	int i, j, k,m,n,t;
	int msel1 = msel - 1, msel2 = msel * (msel1) / 2;
	double *x, *x2, **xyp, **xyn;
	int* mha;
	x = new double[msel];
	if (x == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	x2 = new double[msel];
	if (x2 == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	xyp = new double*[msel1];
	if (xyp == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < msel1; i++)
	{
		xyp[i] = new double[msel1 - i];
		if (xyp[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	xyn = new double* [msel1];
	if (xyn == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	for (i = 0; i < msel1; i++)
	{
		xyn[i] = new double[msel1 - i];
		if (xyn[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	mha = new int[msel];
	if (mha == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }

	for (i = 0; i < msel; i++)mha[i] = 0;	
	k = 0;
	for (i = 0; i < mtot; i++)
	{
		if (a->mot[i] == 1)mha[k++] = i;
	}
	for (i = 0; i < msel1; i++)for (j = 0; j < msel1-i; j++)xyp[i][j] = xyn[i][j] = 0;
	//Positive
	for (i = 0; i < msel; i++)x[i] = x2[i] = 0;
	k = 0;
	for (m = 0; m < msel; m++)
	{
		int mk = mha[m];
		for (j = 0; j < nseqp; j++)
		{
			x[k] += errp[mk][j];
			x2[k] += errp[mk][j] * errp[mk][j];
		}
		x[k] /= nseqp;		
		k++;
	}
	k = 0;
	for (m = 0; m < msel; m++)
	{
		double mnoj1 = sqrt(x2[k] - nseqp * x[k] * x[k]);
		int mk = mha[m];
		t = 0;
		int t1 = k + 1;
		for (n = m + 1; n < msel; n++)
		{
			int nk = mha[n];
			for (j = 0; j < nseqp; j++)
			{
				xyp[k][t] += errp[mk][j] * errp[nk][j];
			}			
			xyp[k][t] -= nseqp * x[k] * x[t1];			
			double mnoj2 = sqrt(x2[t1] - nseqp * x[t1] * x[t1]);
			xyp[k][t] /= (mnoj1 * mnoj2);
			t++;
			t1++;
		}
		k++;
	}
	//Negative
	for (i = 0; i < msel; i++)x[i] = x2[i] = 0;
	k = 0;
	for (m = 0; m < msel; m++)
	{
		int mk = mha[m];
		for (j = 0; j < nseqn; j++)
		{
			x[k] += errn[mk][j];
			x2[k] += errn[mk][j] * errn[mk][j];
		}
		x[k] /= nseqn;		
		k++;
	}
	k = 0;
	for (m = 0; m < msel; m++)
	{
		double mnoj1 = sqrt(x2[k] - nseqn * x[k] * x[k]);
		int mk = mha[m];
		t = 0;
		int t1 = k + 1;
		for (n = m + 1; n < msel; n++)
		{
			int nk = mha[n];
			for (j = 0; j < nseqn; j++)
			{
				xyn[k][t] += errn[mk][j] * errn[nk][j];
			}			
			xyn[k][t] -= nseqn * x[k] * x[t1];
			double mnoj2 = sqrt(x2[t1] - nseqn * x[t1] * x[t1]);
			xyn[k][t] /= (mnoj1 * mnoj2);
			t++;
			t1++;
		}
		k++;
	}
	for (i = 0; i < msel; i++)fprintf(out, "\t%s", nam_mot[mha[i]]);
	fprintf(out, "\t\t");
	for (i = 0; i < msel; i++)fprintf(out, "\t%s", nam_mot[mha[i]]);
	fprintf(out, "\n");
	for (i = 0; i < msel1; i++)
	{
		fprintf(out, "%s",nam_mot[mha[i]]);
		for (j = 0; j <=  i; j++)fprintf(out, "\t");
		for (j = 0; j < msel1 - i; j++)
		{
			fprintf(out, "\t%f", xyp[i][j]);
		}
		fprintf(out, "\t\t");
		fprintf(out, "%s", nam_mot[mha[i]]);
		for (j = 0; j <= i; j++)fprintf(out, "\t");
		for (j = 0; j < msel1 - i; j++)
		{
			fprintf(out, "\t%f", xyn[i][j]);
		}
		fprintf(out, "\n");
	}
	fprintf(out, "%s", nam_mot[mha[msel1]]);
	for (j = 0; j < msel; j++)fprintf(out, "\t");
	fprintf(out, "\t\t");
	fprintf(out, "%s\n", nam_mot[mha[msel1]]);
	delete[] x;
	delete[] x2;
	delete[] mha;
	for (i = 0; i < msel1; i++)delete[] xyp[i];
	delete[] xyp;
	for (i = 0; i < msel1; i++)delete[] xyn[i];
	delete[] xyn;
	return 1;
}
int PearsonExternal(combi* a, combi* b, int msel, int mtot, double** errp, double** errn, int nseqp, int nseqn, double fp2, double &abp, double& abn)
{
	int i, j, k;
	int mha[MMAX], mhb[MMAX];
	for (i = 0; i < MMAX; i++)mha[i] = mhb[i] = 0;
	k = 0;
	for (i = 0; i < mtot; i++)
	{
		if (a->mot[i] == 1)mha[k++] = i;
	}
	k = 0;
	for (i = 0; i < mtot; i++)
	{
		if (b->mot[i] == 1)mhb[k++] = i;
	}
	double xyp =0, xyn =0;
	k = 0;
	for (j = 0; j < nseqp; j++)
	{
		double maxa = fp2;
		for (i = 0; i < msel; i++)
		{
			int ii = mha[i];
			if (errp[ii][j] > maxa)maxa = errp[ii][j];
		}
		double maxb = fp2;
		for (i = 0; i < msel; i++)
		{
			int ii = mhb[i];
			if (errp[ii][j] > maxb)maxb = errp[ii][j];
		}				
		xyp += (maxa - fp2) * (maxb - fp2);
	}
	for (j = 0; j < nseqn; j++)
	{
		double maxa = fp2;
		for (i = 0; i < msel; i++)
		{
			int ii = mha[i];
			if (errn[ii][j] > maxa)maxa = errn[ii][j];
		}
		double maxb = fp2;
		for (i = 0; i < msel; i++)
		{
			int ii = mhb[i];
			if (errn[ii][j] > maxb)maxb = errn[ii][j];
		}
		xyn += (maxa - fp2) * (maxb - fp2);
	}
	abp = (nseqp * xyp - a->xp * b->xp) / sqrt((nseqp * a->x2p - a->xp * a->xp) * (nseqp * b->x2p - b->xp * b->xp));
	abn = (nseqn * xyn - a->xn * b->xn) / sqrt((nseqn * a->x2n - a->xn * a->xn) * (nseqn * b->x2n - b->xn * b->xn));
	return 1;
}
int main(int argc, char* argv[])
{
	char filei_tabp[300], filei_tabn[300], filei_mot_names[300], filei_class_names[300], fileo_prc[300], d[50000], fileo[300], fileo_corr_ext[300], fileo_corr_int[300], file_log[300];
	int i, j, k, m, msel, mtot = 0, nseqp = 0, nseqn = 0;
	FILE* inp, * inn, * inmotnam, * inclassnam, * outh, *outlog, *out, *out_corr_ext, * out_corr_int;

	if (argc != 12)
	{
		printf("Syntax: 1filei_tabp 2filei_tabn 3filei_motif_names 4filei_motif_class_names 5int motif_count 6double ERRthresh ");
		printf("7fileo_prc 8fileo_corr_external 9fileo_corr_internal 10fileo_results 11file_log");
		exit(1);
	}
	strcpy(filei_tabp, argv[1]);
	strcpy(filei_tabn, argv[2]);
	strcpy(filei_mot_names, argv[3]);
	strcpy(filei_class_names, argv[4]);
	msel = atoi(argv[5]);//4islo motivov
	double fp2 = atof(argv[6]); //ERR threshold for pAUPRC	
	strcpy(fileo_prc, argv[7]);	
	strcpy(fileo_corr_ext, argv[8]);
	strcpy(fileo_corr_int, argv[9]);
	strcpy(fileo, argv[10]);
	strcpy(file_log, argv[11]);
	strcat(fileo, "_");
	strcat(file_log, "_");
	strcat(fileo_prc, "_");
	strcat(fileo_corr_ext, "_");
	strcat(fileo_corr_int, "_");
	strcat(fileo_prc, argv[5]);
	strcat(fileo, argv[5]);
	strcat(file_log, argv[5]);
	strcat(fileo_corr_ext, argv[5]);
	strcat(fileo_corr_int, argv[5]);

	srand((unsigned)time(NULL));
	if ((inp = fopen(filei_tabp, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", filei_tabp);
		exit(1);
	}
	while (fgets(d, sizeof(d), inp) != NULL)
	{
		if (isdigit(d[0]) == 0) break;
		nseqp++;
	}
	rewind(inp);
	if ((inn = fopen(filei_tabn, "rt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", filei_tabn);
		exit(1);
	}
	while (fgets(d, sizeof(d), inn) != NULL)
	{
		if (isdigit(d[0]) == 0) break;
		nseqn++;
	}
	rewind(inn);

	int nseqtot = nseqp + nseqn;
	double prec_exp = (double)0.5, nseqrat = (double)nseqp / nseqn;
	qbs* seqtot;
	seqtot = new qbs[nseqtot];
	if (seqtot == NULL) { printf("Out of memory..."); return -1; };

	if ((inmotnam = fopen(filei_mot_names, "rt")) == NULL)
	{
		printf("Input file %s can't be opened\n", filei_mot_names);
		return -1;
	}
	if ((inclassnam = fopen(filei_class_names, "rt")) == NULL)
	{
		printf("Input file %s can't be opened\n", filei_class_names);
		return -1;
	}
	while (fgets(d, sizeof(d), inmotnam) != NULL)
	{
		if (*d != '\n' && *d != '\t')mtot++;
	}
	rewind(inmotnam);
	if ((inclassnam = fopen(filei_class_names, "rt")) == NULL)
	{
		puts("Input file can't be opened");
		return -1;
	}
	{
		int nclass = 0;
		while (fgets(d, sizeof(d), inclassnam) != NULL)
		{
			if (*d != '\n' && *d != '\t')nclass++;
		}
		rewind(inclassnam);
		if (nclass != mtot)
		{
			printf("Files %s & %s are not accordant %d & %d\n", filei_mot_names, filei_class_names, mtot, nclass);
			fclose(inclassnam);
			exit(1);
		}
	}	
	if (mtot > MMAX)
	{
		printf("Number of motifs %d above the upper limit %d", mtot, MMAX);
		exit(1);
	}
	char** motnames;
	char** class_names;
	{
		const size_t lens = 50;
		motnames = new char* [mtot];
		if (motnames == NULL) { printf("Out of memory..."); return -1; };
		for (i = 0; i < mtot; i++)
		{
			motnames[i] = new char[lens];
			if (motnames[i] == NULL) { puts("Out of memory..."); exit(1); }
		}
		size_t sizemot = lens * sizeof(motnames[i][0]);
		i = 0;
		for (i = 0; i < mtot; i++)
		{
			memset(motnames[i], '\0', sizeof(sizemot));
			fgets(d, sizeof(d), inmotnam);
			if (*d != '\n' && *d != '\t')
			{
				DelChar(d, '\n');
				int dlen = (int)strlen(d);
				strncpy(motnames[i], d, dlen);
				motnames[i][dlen] = '\0';
			}
			else
			{
				printf("Wrong format %s\n", filei_mot_names); return(-1);
			}
		}		
		fclose(inmotnam);
		class_names = new char* [mtot];
		if (class_names == NULL) { printf("Out of memory..."); return -1; };
		for (i = 0; i < mtot; i++)
		{
			class_names[i] = new char[lens];
			if (class_names[i] == NULL) { puts("Out of memory..."); exit(1); }
		}
		sizemot = lens * sizeof(class_names[i][0]);
		i = 0;
		for (i = 0; i < mtot; i++)
		{
			memset(class_names[i], '\0', sizeof(sizemot));
			fgets(d, sizeof(d), inclassnam);
			if (*d != '\n' && *d != '\t')
			{
				DelChar(d, '\n');
				int dlen = (int)strlen(d);
				strncpy(class_names[i], d, dlen);
				class_names[i][dlen] = '\0';
			}
			else
			{
				printf("Wrong format %s\n", filei_class_names); return(-1);
			}
		}

	}	
	fclose(inclassnam);
	double** errp;
	errp = new double* [mtot];
	if (errp == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < mtot; i++)
	{
		errp[i] = new double[nseqp];
		if (errp[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	double** errn;
	errn = new double* [mtot];
	if (errn == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < mtot; i++)
	{
		errn[i] = new double[nseqn];
		if (errn[i] == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	}
	{
		char sep = '\t';
		for (j = 0; j < nseqp; j++)
		{
			fgets(d, sizeof(d), inp);
			for (i = 0; i < mtot; i++)
			{
				double test = UnderStol(d, i, sep);
				if (test == -1) { printf("Wrong format %s\n", filei_tabp); return(-1); }
				errp[i][j] = test;
			}
		}
		fclose(inp);
		for (j = 0; j < nseqn; j++)
		{
			fgets(d, sizeof(d), inn);
			for (i = 0; i < mtot; i++)
			{
				double test = UnderStol(d, i, sep);
				if (test == -1) { printf("Wrong format %s\n", filei_tabn); return(-1); }
				errn[i][j] = test;
			}
		}
		fclose(inn);
	}
	int pool_act, elit_act;
	if (msel == 1)
	{
		pool_act = mtot;
		elit_act = mtot;
	}
	else
	{
		pool_act = POOL;
		elit_act = ELIT;
		if (msel == 2)
		{
			int mpa = mtot * (mtot - 1) / 2;
			if (mpa < pool_act)
			{
				pool_act = elit_act = 3*mpa / 4;
			}
		}
	}
	pop = new combi[pool_act];
	if (pop == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	int point_test;
	if (msel <= 5) point_test = elit_act - 1;
	else
	{
		int limit = pool_act / msel;
		point_test = (5* (elit_act - 1) + limit * (msel - 5)) / msel;
	}
	if ((out = fopen(fileo, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", fileo);
		exit(1);
	}
	if ((outlog = fopen(file_log, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", file_log);
		exit(1);
	}
	double* fit_prev_iter;
	fit_prev_iter = new double[elit_act];
	if (fit_prev_iter == NULL) { fprintf(stderr, "Error: Out of memory..."); return -1; }
	int empty_cycle = 0;
	det1.init0(mtot);
	for (i = 0; i < pool_act; i++)pop[i].init0(mtot);
	for (i = 0; i < 2; i++)det2[i].init0(mtot);
	if (msel > 1)
	{
		for (i = 0; i < pool_act; i++)
		{
			int gom, trys = 0;
			do
			{
				gom = 0;
				//	printf("A \n");
				det1.init1(msel, mtot);
				//printf("B \n");
				for (j = 0; j < i; j++)
				{
					if (GomTown(det1, pop[j], mtot) == 1)
					{
						gom = 1;
						break;
					}
				}
				//printf("C \n");
				if (gom == 1)trys++;
			} while (gom == 1);
			//printf("D \n");
			//det1.fprint_all(mtot);
			//printf("E \n");
			det1.get_copy(&pop[i], mtot);
			//printf("F \n");
			int test = EvalFit(&pop[i], msel, mtot, errp, errn, nseqp, nseqn, nseqtot, nseqrat, prec_exp, fp2, seqtot);
			//printf("Rank %d\tTry %d\t", i + 1, trys + 1);
			//pop[i].fprint_all(mtot);
		}
	}
	else
	{
		{			
			for (i = 0; i < mtot; i++)
			{
				pop[i].fit = 0;
				for (j = 0; j < mtot; j++)
				{
					if (i == j)pop[i].mot[j] = 1;
					else pop[i].mot[j] = 0;
				}
				int test = EvalFit(&pop[i], msel, mtot, errp, errn, nseqp, nseqn, nseqtot, nseqrat, prec_exp, fp2, seqtot);				
			}			
		}
	}
	qsort((void*)pop, pool_act, sizeof(pop[0]), compare_pop);
	for (i = 0; i < elit_act; i++)fit_prev_iter[i] = pop[i].fit;
	printf("Initiated & sorted: %d motifs, population of %d motif groups (elit of %d groups), %d motif(s) in each group\tPoint test %d\n", mtot, pool_act, elit_act, msel, point_test);
	fprintf(outlog,"Initiated & sorted: %d motifs, population of %d motif groups (elit of %d groups), %d motif(s) in each group\tPoint test %d\n",mtot,pool_act,elit_act,msel, point_test);
//	for (i = 0; i < elit_act; i++)pop[i].fprint_all(mtot,outlog);
	time_t tnow = time(NULL);
	fprintf(outlog, "%s", ctime(&tnow));
	if(msel > 1)
	{
		fprintf(outlog,"GA started\n");
		printf("GA started\n");
		int pool_last = pool_act - 1;
		int pair_all_max = 0;
		for (j = 0; j < pool_act; j++)
		{
			for (k = j + 1; k < pool_act; k++)
			{
				pair_all_max++;
			}
		}
		int** pair_d; // pair 1 & 2 population ranks of individuals
		pair_d = new int* [pair_all_max];
		if (pair_d == NULL) return -1;
		for (m = 0; m < pair_all_max; m++)
		{
			pair_d[m] = new int[2];
			if (pair_d[m] == NULL) return -1;
		}
		m = 0;
		for (j = 0; j < pool_act; j++)
		{
			for (k = j + 1; k < pool_act; k++)
			{
				pair_d[m][0] = j;
				pair_d[m][1] = k;
				m++;
			}
		}
		int* pair_take;
		pair_take = new int[pair_all_max];
		if (pair_take == NULL) { puts("Out of memory..."); exit(1); }
		for (k = 0; k < pair_all_max; k++)pair_take[k] = k;
		//GA
		int m_success_no_max = pool_last / 2, r_success_no_max = 1;
		int iter = 1;
		int total_m_success, total_r_success;
		double fit_prev = pop[0].fit, fit_med = pop[point_test].fit, ga_exit = 1;
		do
		{
			total_m_success = 0;
			total_r_success = 0;
			fprintf(outlog, "Mutation %d\n", iter);
			for (i = 0; i < pool_act; i++)
			{
				int m_success_no = 0;
				if (pop[i].fit == 0)continue;
				do
				{
					pop[i].get_copy(&det1, mtot);
					int mut = MutCo(&det1, msel, mtot);
					if (mut == -1)
					{
						m_success_no++;
						continue;
					}
					int gom = 0;
					for (j = 0; j < pool_act; j++)
					{
						if (GomTown(det1, pop[j], mtot) == 1) { gom = 1; break; }
					}
					if (gom == 1)
					{
						m_success_no++;
						continue;
					}
					//pop[i].fprint_all(mtot);				
					double fit = EvalFit(&det1, msel, mtot, errp, errn, nseqp, nseqn, nseqtot, nseqrat, prec_exp, fp2, seqtot);
					//det1.fprint_all(mtot);
					if (fit == -1)m_success_no++;
					else
					{
						if (det1.fit > pop[i].fit)
						{
							det1.get_copy(&pop[i], mtot);
							m_success_no = 0;
							total_m_success++;
							/*int check = pop[i].check(n_olig, max_olig);
							if (check == -1)
							{
								fprintf(outlog, "After Mut!\n");
								exit(1);
							}*/
						}
						else m_success_no++;
					}
				} while (m_success_no < m_success_no_max);
				//printf("Mutation %d\tSuccess %d\n", i + 1, total_m_success);
			}
			if(total_m_success>0)qsort((void*)pop, pool_act, sizeof(pop[0]), compare_pop);
			fprintf(outlog,"Mut Iteration %d\t%d\tRat1st %f\tRatPointTest %f\tAt rank %d\n", iter, total_m_success, pop[0].fit / fit_prev, pop[point_test].fit / fit_med, point_test);
			printf("Mut Iteration %d\t%d\tRat1st %f\tRatPointTest %f\tAt rank %d\n", iter, total_m_success, pop[0].fit/fit_prev,pop[point_test].fit/fit_med,point_test);			
			//for (i = 0; i < elit_act; i++)pop[i].fprint_all(mtot, outlog);
			//fprintf(outlog, "First\t"); pop[0].fprint_all(mtot, outlog);
			//fprintf(outlog, "PointTest %d\t", point_test); pop[point_test].fprint_all(mtot, outlog);
			//Rec
		//	fprintf(outlog, "Recombination %d\n", iter);
		//	printf("Recombination %d\n", iter);
			BigMixI(pair_take, pair_all_max);
			for (i = 0; i < pair_all_max; i++)
			{
				int kk[2];
				for (m = 0; m < 2; m++)kk[m] = pair_d[pair_take[i]][m];
				for (m = 0; m < 2; m++)
				{
					pop[kk[m]].get_copy(&det2[m], mtot);
				}
				double fit_parent_max;
				fit_parent_max = Max(det2[0].fit, det2[1].fit);
				int r_success_no = 0;
				do
				{
					int rec = RecCo(&det2[0], &det2[1], msel, msel, mtot);
					//for (k = 0; k < 2; k++)det2[k].fprint_all(mtot,outlog);
					if (rec == -1)
					{
						r_success_no++;
						continue;
					}
					int gom = 0;
					for (k = 0; k < 2; k++)
					{
						for (j = 0; j < pool_act; j++)
						{
							if (GomTown(det2[k], pop[j], mtot) == 1) { gom = 1; break; }
						}
						if (gom == 1)break;
					}
					if (gom == 1)
					{
						r_success_no++;
						break;
					}
					for (k = 0; k < 2; k++)EvalFit(&det2[k], msel, mtot, errp, errn, nseqp, nseqn, nseqtot, nseqrat, prec_exp, fp2, seqtot);
					double fitmax = Max(det2[0].fit, det2[1].fit);					
					if (fitmax > fit_parent_max)
					{						
						for (k = 0; k < 2; k++)det2[k].get_copy(&pop[kk[k]], mtot);
						total_r_success++;
					}
					else r_success_no++;
				} while (r_success_no < r_success_no_max);
				//printf("Recombination %d\tSuccess %d\n", i + 1, total_r_success);
			}
			if (total_r_success > 0)qsort((void*)pop, pool_act, sizeof(pop[0]), compare_pop);			
			double rat_fit = pop[0].fit / fit_prev, rat_med = pop[point_test].fit / fit_med;
			fit_prev = pop[0].fit;
			fit_med = pop[point_test].fit;
			int change_rank = point_test + 1;
			for (i = 0; i < elit_act; i++)
			{
				if(pop[i].fit!=fit_prev_iter[i])
				{
					change_rank = i;
					break;
				}
			}
			for (i = 0; i < elit_act; i++)fit_prev_iter[i] = pop[i].fit;
			if (rat_fit <= ga_exit && change_rank >= point_test)empty_cycle++;
			//if (rat_fit <= ga_exit && rat_med <= ga_exit)empty_cycle++;
			else empty_cycle = 0;
			fprintf(outlog, "Rec Iteration %d\t%d\tRat1st %f\tRatPointTest %f\tat rank %d\tChangeRank %d\tEmptyCycle %d\n", iter, total_r_success,rat_fit,rat_med, point_test,change_rank,empty_cycle);
			printf("Rec Iteration %d\t%d\tRat1st %f\tRatPointTest %f\tAt rank %d\tChangeRank %d\tEmptyCycle %d\n", iter, total_r_success, rat_fit, rat_med, point_test, change_rank, empty_cycle);
			//for (i = 0; i < elit_act; i++)pop[i].fprint_all(mtot, outlog);
		//	fprintf(outlog, "First\t"); pop[0].fprint_all(mtot, outlog);
			//fprintf(outlog, "Point Test%d\t",point_test); pop[point_test].fprint_all(mtot, outlog);
	//		fprintf(outlog, "Iteration %d Mut %d Rec %d FitRatio %f\tRatPointTest %f\n", iter, total_m_success, total_r_success, rat_fit,rat_med);
			printf("Iteration %d Mut %d Rec %d\n", iter, total_m_success, total_r_success);
			tnow = time(NULL);
			fprintf(outlog, "%s", ctime(&tnow));
			iter++;
		} while (empty_cycle < 2);
		for (i = 0; i < pair_all_max; i++)
		{
			delete[] pair_d[i];
		}
		delete[] pair_d;
		delete[] pair_take;
		/*fprintf(outlog, "Top %d groups\n", elit_act);
		for (i = 0; i < elit_act; i++)
		{
			fprintf(outlog, "Rank %d\t", i + 1);
			pop[i].fprint_all(mtot,outlog);
		}*/
		/*fprintf(outlog, "Bottom %d groups\n", elit_act);
		for (i = pool_act - elit_act; i < pool_act; i++)
		{
			fprintf(outlog, "Rank %d\t", i + 1);
			pop[i].fprint_all(mtot, outlog);
		}*/
	}
	if ((outh = fopen(fileo_prc, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", fileo_prc);
		exit(1);
	}
	for (i = 0; i < elit_act; i++)
	{
		Precision(&pop[i], msel, mtot, errp, errn, nseqp, nseqn, nseqtot, nseqrat, prec_exp, fp2, seqtot,outh, motnames,i+1);
	}
	fclose(outh);
	for (i = 0; i < elit_act; i++)
	{
		fprintf(out,"Rank %d\t", i + 1);
		pop[i].fprint_nam(mtot, motnames, class_names,out);
	}
	int pair_elit = elit_act * (elit_act - 1) / 2;
	double* abp;
	abp = new double[pair_elit];
	if (abp == NULL) { puts("Out of memory..."); exit(1); }
	double* abn;
	abn = new double[pair_elit];
	if (abn == NULL) { puts("Out of memory..."); exit(1); }
	for (i = 0; i < elit_act; i++)Modules(&pop[i], msel, mtot, errp, errn, nseqp, nseqn, fp2);
	k = 0;
	for (i = 0; i < elit_act; i++)
	{
		for (j = i+1; j < elit_act; j++)
		{
			PearsonExternal(&pop[i], &pop[j], msel, mtot, errp, errn, nseqp, nseqn, fp2, abp[k], abn[k]);
			k++;
		}
	}
	fclose(out);
	if (msel > 1)
	{
		if ((out_corr_int = fopen(fileo_corr_int, "wt")) == NULL)
		{
			printf("Input file %s can't be opened!\n", fileo_corr_int);
			exit(1);
		}
		for (i = 0; i < elit_act; i++)
		{
			fprintf(out_corr_int, "Rank %d\tPositive\t", i + 1);
			for (j = 0; j < msel; j++)fprintf(out_corr_int, "\t");
			fprintf(out_corr_int, "Rank %d\tNegative\n", i + 1);
			PearsonInternal(&pop[i], msel, mtot, errp, errn, nseqp, nseqn, motnames, out_corr_int);
		}
		fclose(out_corr_int);
	}
	if ((out_corr_ext = fopen(fileo_corr_ext, "wt")) == NULL)
	{
		printf("Input file %s can't be opened!\n", fileo_corr_ext);
		exit(1);
	}
	fprintf(out_corr_ext,"Positive\n");
	for (i = 0; i < elit_act; i++)fprintf(out_corr_ext,"\t%d", i + 1);
	fprintf(out_corr_ext, "\n");
	k = 0;
	for (i = 0; i < elit_act; i++)
	{
		fprintf(out_corr_ext, "%d", i + 1);
		for (j = 0; j <= i; j++)fprintf(out_corr_ext, "\t");
		for (j = i + 1; j < elit_act; j++)
		{
			fprintf(out_corr_ext, "\t%.4f", abp[k++]);
		}
		fprintf(out_corr_ext, "\n");
	}
	fprintf(out_corr_ext,"Negative\n");
	for (i = 0; i < elit_act; i++)fprintf(out_corr_ext,"\t%d", i + 1);
	fprintf(out_corr_ext,"\n");
	k = 0;
	for (i = 0; i < elit_act; i++)
	{
		fprintf(out_corr_ext,"%d", i + 1);
		for (j = 0; j <= i; j++)fprintf(out_corr_ext,"\t");
		for (j = i + 1; j < elit_act; j++)
		{
			fprintf(out_corr_ext,"\t%.4f", abn[k++]);
		}
		fprintf(out_corr_ext,"\n");
	}
	fclose(out_corr_ext);
	fclose(outlog);
	for (i = 0; i < pool_act; i++)pop[i].term();
	for (k = 0; k < mtot; k++)
	{
		delete[] errp[k];
		delete[] errn[k];
		delete[] motnames[k];
	}
	delete[] pop;
	delete[] abp;
	delete[] abn;
	delete[] errp;
	delete[] errn;
	delete[] motnames;
	delete[] seqtot;
	delete[] fit_prev_iter;
}
#define _CRT_SECURE_NO_WARNINGS

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include  <time.h>
#include  <ctype.h>
#define SEQLEN 50000
//take input file and print in output only lines that correspond to that listed in input list_file 0 = no print, 1 = print
int main(int argc, char *argv[])
{
	FILE *in, *innum, *out;
	char d[SEQLEN], s[40];

	if(argc!=4)
	{
		printf ("%s 1file_input 2file_line_list 3file_output",argv[0]);//-shift_min -shift_max
        exit (1);
	}
	if((in=fopen(argv[1],"rt"))==NULL)
	{
		printf("Input file %s can't be opened\n",argv[1]);
		exit(1);
	}
	if((innum=fopen(argv[2],"rt"))==NULL)
	{
		printf("Input file %s can't be opened\n",argv[2]);
		exit(1);
	}
	if((out=fopen(argv[3],"wt"))==NULL)
	{
		printf("Output file %s can't be opened\n",argv[3]);
		exit(1);
	}
	while(fgets(s,sizeof(s),innum)!=NULL)
	{
		if(isdigit(s[0]))
		{
			int val=atoi(s);
			if(fgets(d,sizeof(d),in)==NULL)
			{
				printf("Unexpected end of file %s\n", argv[1]);
				exit(1);
			}
			if(val==1)fprintf(out,"%s",d);
		}
	}
	fclose(in);
	fclose(out);
	return 1;
}
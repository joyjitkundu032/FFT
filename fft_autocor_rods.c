#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <fftw3.h>
#include</home/joyjit/mt19937ar.h>
#include</home/joyjit/mt19937ar.c>

#define K 7  /* length of K-mer*/
#define L 336	/*Number of sites in one direction, should be even*/
#define N (L*L) /* The total number of sites */
#define DIR 2 /*Number of positive directions */
#define NN 2 /*Number of nearest neighbours in forward direction */
#define T 8000000	/*# of MC steps for equilibration*/
#define BLOCKS 20 /*Number of blocks in which average is taken*/
#define AVE 10000 /* # of readings in each block */
#define PSTART 0.7 /*starting probability*/
#define PDIFF 0.1 /*increment in probability*/
#define PEND (PSTART+0.5*PDIFF) /*ending probability */
#define GAP2 50 /*GAP after which output is written during equilibration */
#define INITIAL_FLAG 0 /*0: all empty 1: all filled */
#define BINS 101 /*Keep it odd*/
#define BINSIZE (2.0/(1.0*BINS))
#define STEP 10
#define R 2097152
#define NUM 1000

FILE *fp2;
int lat[N],horizontal[N],vertical[N],lasth,lastv,ln[N],rn[N],bn[N],tn[N];
int starth,finalh,endh,startv,finalv,endv;
double prob,count[BINS],periodic[K],accept[L+1],avgcorr_rho[R/2],avgcorr_Q[R/2];
double auto_rho[R/2],auto_Q[R/2],rho[R/2],Q[R/2];
double meanrho[BLOCKS],meanabs[BLOCKS],fft_rho[R],fft_Q[R];
double meanm1[BLOCKS],meanm2[BLOCKS],meanm4[BLOCKS];
double mass[BINS],factor[BINS],bincount[BINS];
char outfile1[100],outfile2[100],outfile3[100],fft_corr[100],avgvalue[100],Kspace_corr[100];
char readfile1[100],readfile2[100];
int k,l,m,n,x,y;
double avg_rho,avg_Q;

void complex_1d_rhodft()
{
	int i,j;
	double sum_rho;
	fftw_complex* f, * g, * h;
	fftw_plan     plan2_f, plan2_b;
	f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * R);
	g = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * R);
	h = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * R);
	sum_rho=0;
	for(i=0;i<R/2;i++)
	{
		f[i][0]=rho[i];
		f[i][1]=0.0;
		sum_rho=sum_rho+rho[i];
	}
	for(i=R/2;i<R;i++)
	{
		f[i][0]=0.0;
		f[i][1]=0.0;
	}

	plan2_f = fftw_plan_dft_1d(R, f, g, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan2_f);
	
	for(i=0;i<R;i++)
	{
		g[i][0]=(g[i][0]*g[i][0]+g[i][1]*g[i][1]);
		g[i][1]=0.0;
		fft_rho[i]=fft_rho[i]+g[i][0];
	}
 
	plan2_b = fftw_plan_dft_1d(R, g, h, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan2_b);
	
	for(i=0;i<R;i++)
	{
		h[i][0]=h[i][0]/pow(R,2.0);
		h[i][1]=h[i][1]/pow(R,2.0);
		//avgcorr_rho[i]=avgcorr_rho[i]+h[i][0];
	}
	for(i=0;i<R/2;i++)
	{
		h[i][0]=1.0*h[i][0]*R/(R/2-i);
		avgcorr_rho[i]=avgcorr_rho[i]+h[i][0];
	}

	fftw_destroy_plan(plan2_f);
	fftw_destroy_plan(plan2_b);
	fftw_free(f);
	fftw_free(g);
	fftw_free(h);

	/*printf("%e\n",sumavg2);*/
	sum_rho=2.0*sum_rho/(1.0*R);
	avg_rho=avg_rho+sum_rho;

}
void complex_1d_orderdft()
{
	int i,j;
	double sum_Q;
	fftw_complex* f, * g, * h;
	fftw_plan     plan2_f, plan2_b;
	f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * R);
	g = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * R);
	h = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * R);
	sum_Q=0;
	for(i=0;i<R/2;i++)
	{
		f[i][0]=Q[i];
		f[i][1]=0.0;
		sum_Q=sum_Q+Q[i];
	}
	for(i=R/2;i<R;i++)
	{
		f[i][0]=0.0;
		f[i][1]=0.0;
	}

	plan2_f = fftw_plan_dft_1d(R, f, g, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan2_f);
	
	for(i=0;i<R;i++)
	{
		g[i][0]=(g[i][0]*g[i][0]+g[i][1]*g[i][1]);
		g[i][1]=0.0; 
		fft_Q[i]=fft_Q[i]+g[i][0];
	}
		
	plan2_b = fftw_plan_dft_1d(R, g, h, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan2_b);
	
	for(i=0;i<R;i++)
	{
		h[i][0]=h[i][0]/pow(R,2.0);
		h[i][1]=h[i][1]/pow(R,2.0);
		//avgcorr_Q[i]=avgcorr_Q[i]+h[i][0];
	}
	for(i=0;i<R/2;i++)
	{
		h[i][0]=1.0*R*h[i][0]/(R/2-i);
		avgcorr_Q[i]=avgcorr_Q[i]+h[i][0];
	}

	fftw_destroy_plan(plan2_f);
	fftw_destroy_plan(plan2_b);
	fftw_free(f);
	fftw_free(g);
	fftw_free(h);

	/*printf("%e\n",sumavg2);*/
	sum_Q=2.0*sum_Q/(1.0*R);
	avg_Q=avg_Q+sum_Q;

}
void print_config()
{
	/* Prints out the configuration */
	int x,y;

	for(y=L-1;y>=0;y--)
	{
		for(x=0;x<L;x++)
			printf("%2d",lat[find_site(x,y)]);
		printf("\n");
	}
	printf("***************************\n");
	sleep(2);
}

int find_site(int x, int y)
{
	/* Given x, y it gives site index */
	x=(x+L) % L;
	y=(y+L) % L;
	return x+y*L;
}

void take_input()
{
	/*reads in seed for random number generator */
	long seedval;

	/*printf("\nEnter value of seed : ");*/
	//scanf("%ld",&seedval);
	seedval=715362;
	init_genrand(seedval);
}

void initialize()
{
	/* initializes nbr list and output file names */
	int i,j,x,y;
	double tmp,tmp1;

	for(y=0;y<L;y++)
	{
		for(x=0;x<L;x++)
		{
			i=find_site(x,y);
			j=find_site(x-1,y); ln[i]=j;
			j=find_site(x+1,y); rn[i]=j;
			j=find_site(x,y-1); bn[i]=j;
			j=find_site(x,y+1); tn[i]=j;
		}
	}
	if(INITIAL_FLAG==0)
		sprintf(outfile1,"emptyK%dL%d",K,L);
	if(INITIAL_FLAG==1)
		sprintf(outfile1,"filledK%dL%d",K,L);
	/*initializing bin parameters*/

	for(j=0;j<BINS;j++)
	{
		bincount[j]=0;
		mass[j]=0;
	}
	for(j=-L*L;j<=L*L;j=j+K)
	{
		tmp=(1.0*j)/(1.0*N);
		i=floor((tmp+1.0)/BINSIZE);
		if(j==-N)
			i=0;
		if(j==N)
			i=BINS-1;
		bincount[i]++;
		mass[i]=mass[i]+tmp;
	}
	tmp1=(K*1.0)/(1.0*N);
	for(j=0;j<BINS;j++)
	{
		mass[j]=mass[j]/bincount[j];
		factor[j]=tmp1*bincount[j];
	}
}

void lat_init()
{
	/* Initializes quantities that have to initialized for every value
	 * of probability p */
	int i;
	double x;
	FILE *fp;

	lastv=0;lasth=0; 
	for(i=0;i<BLOCKS;i++)
	{
		meanrho[i]=0;meanabs[i]=0;meanm1[i]=0; meanm2[i]=0; meanm4[i]=0;
	}
	if((L % 2) !=0)
	{
		printf("ERROR IN DIVISIBILITY\n");
		exit(0);
	}
	for(i=0;i<BINS;i++)
		count[i]=0;
	if(INITIAL_FLAG==0)
	{
		for(i=0;i<N;i++)
			lat[i]=0;
		sprintf(outfile2,"emptyK%dL%dP%5.4lf_t",K,L,prob);
		sprintf(outfile3,"emptyK%dL%dP%5.4lf_p",K,L,prob);
	}
	if(INITIAL_FLAG==1)
	{
		for(i=0;i<N/2;i++)
		{
			lat[i]=1;
			horizontal[lasth]=i;
			lasth++;
		}
		for(i=N/2;i<N;i++)
		{
			lat[i]=2;
			vertical[lastv]=i;
			lastv++;
		}
		sprintf(outfile2,"filledK%dL%dP%5.4lf_t",K,L,prob);
		sprintf(outfile3,"filledK%dL%dP%5.4lf_p",K,L,prob);
	}
	fp=fopen(outfile2,"w");
	fprintf(fp,"#t rho abs(m) m\n");
	fclose(fp);

	sprintf(readfile1,"acceptprobK%dL%dP%5.4lf",K,L,prob);
	fp=fopen(readfile1,"r");
	if(fp==NULL)
	{
		printf("The FILE [%s] DOES NOT EXIST\n",readfile1);
		exit(0);
	}
	while(fscanf(fp,"%d%lf",&i,&x)!=EOF)
		accept[i]=x;
	fclose(fp);
	if(i!=L)
	{
		printf("Error in FILE [%s]\n",readfile1);
		exit(0);
	}

	sprintf(readfile2,"periodicprobK%dL%dP%5.4lf",K,L,prob);
	fp=fopen(readfile2,"r");
	if(fp==NULL)
	{
		printf("The FILE [%s] DOES NOT EXIST\n",readfile2);
		exit(0);
	}
	while(fscanf(fp,"%d%lf",&i,&x)!=EOF)
		periodic[i]=x;
	fclose(fp);
	if(i!=K-1)
	{
		printf("Error in FILE [%s]\n",readfile2);
		exit(0);
	}

}

void remove_hor()
{
	/* removes all horizontal kmers */
	int i;

	for(i=0;i<lasth;i++)
		lat[horizontal[i]]=0;
	lasth=0;
}

void remove_ver()
{
	/* removes all vertical kmers */
	int i;

	for(i=0;i<lastv;i++)
		lat[vertical[i]]=0;
	lastv=0;
}

void deposit_hor(int i)
{
	/*puts a horizontal kmer with head at i*/

	int j;
	
	for(j=1;j<=K;j++)
	{
		lat[i]=1;
		horizontal[lasth]=i;
		lasth++;
		i=rn[i];
	}
}

void deposit_ver(int i)
{
	/*puts a vertical kmer with head at i*/

	int j;
	
	for(j=1;j<=K;j++)
	{
		lat[i]=2;
		vertical[lastv]=i;
		lastv++;
		i=tn[i];
	}
}

void fill_periodic_hor(int j)
{
	/*Given an empty horizontal line, the function checks the
	 * occupation probability of first K-1 sites */

	int i,k;

	starth=j;finalh=ln[j];
	for(i=0;i<K-1;i++)
	{
		if(genrand_real3() < periodic[i])
		{
			for(k=1;k<=K;k++)
			{
				lat[starth]=1;
				horizontal[lasth]=starth;
				lasth++;
				starth=rn[starth];
			}
			return;
		}
		else
		{
			finalh=starth;
			starth=rn[starth];
		}
	}
	return;
}

void fill_periodic_ver(int j)
{
	/*Given an empty vertical line, the function checks the
	 * occupation probability of first K-1 sites */

	int i,k;

	startv=j;finalv=bn[j];
	for(i=0;i<K-1;i++)
	{
		if(genrand_real3() < periodic[i])
		{
			for(k=1;k<=K;k++)
			{
				lat[startv]=2;
				vertical[lastv]=startv;
				lastv++;
				startv=tn[startv];
			}
			return;
		}
		else
		{
			finalv=startv;
			startv=tn[startv];
		}
	}
	return;
}

int find_starthfinalh(int row)
{
	/*in row 'row' finds out starth and finalh*/

	int i;

	starth=row*L; finalh=row*L;
	if(lat[starth]!=0)
	{
		while(lat[starth]!=0)
		{
			starth=rn[starth];
			if(starth==row*L)
				return 0;
		}
		while(lat[finalh]!=0)
			finalh=ln[finalh];
		return 1;
	}

	while(lat[ln[starth]]==0)
	{
		starth=ln[starth];
		if(starth==row*L)
		{
			fill_periodic_hor(starth);
			return 1;
		}
	}
	finalh=ln[starth];
	while(lat[finalh]!=0)
		finalh=ln[finalh];
	return 1;
}

void fill_row(int row)
{
	/*fills row 'row' with horizontal kmers*/

	int i,len;

	if(find_starthfinalh(row)==0)
		return;
	do
	{
		endh=starth;len=1;
		while(lat[rn[endh]]==0)
		{
			endh=rn[endh];
			len++;
			if(endh==finalh)
				break;
		}
		while(len>=K)
		{
			if(genrand_real3() < accept[len])
			{
				deposit_hor(starth);
				for(i=0;i<K;i++)
					starth=rn[starth];
				len=len-K;
			}
			else
			{
				starth=rn[starth];
				len--;
			}
		}
		if(endh==finalh)
			return;
		starth=rn[endh];
		while(lat[starth]!=0)
			starth=rn[starth];
	}
	while(endh!=finalh);
}

int find_startvfinalv(int col)
{
	/*in row 'row' finds out startv and finalv*/

	int i;

	startv=col; finalv=col;
	if(lat[startv]!=0)
	{
		while(lat[startv]!=0)
		{
			startv=tn[startv];
			if(startv==col)
				return 0;
		}
		while(lat[finalv]!=0)
			finalv=bn[finalv];
		return 1;
	}

	while(lat[bn[startv]]==0)
	{
		startv=bn[startv];
		if(startv==col)
		{
			fill_periodic_ver(startv);
			return 1;
		}
	}
	finalv=bn[startv];
	while(lat[finalv]!=0)
		finalv=bn[finalv];
	return 1;
}

void fill_col(int col)
{
	/*fills col 'col' with horizontal kmers*/

	int i,len;

	if(find_startvfinalv(col)==0)
		return;
	do
	{
		endv=startv;len=1;
		while(lat[tn[endv]]==0)
		{
			endv=tn[endv];
			len++;
			if(endv==finalv)
				break;
		}
		while(len>=K)
		{
			if(genrand_real3() < accept[len])
			{
				deposit_ver(startv);
				for(i=0;i<K;i++)
					startv=tn[startv];
				len=len-K;
			}
			else
			{
				startv=tn[startv];
				len--;
			}
		}
		if(endv==finalv)
			return;
		startv=tn[endv];
		while(lat[startv]!=0)
			startv=tn[startv];
	}
	while(endv!=finalv);
}

void evolve()
{
	int row,col,i;

	remove_hor();
	for(row=0;row<L;row++)
		fill_row(row);


	remove_ver();
	for(col=0;col<L;col++)
		fill_col(col);
}

void chk_config()
{
	int i,j;
	j=0;
	for (i=0;i<N;i++)
	{
		if(lat[i]!=0)
			j++;
		if((lat[i]>= DIR)||(lat[i]<0))
		{
			printf("ERROR IN STATE\n");
			sleep(3);
		}
	}
	if(j!=(lastv+lasth))
	{
		printf("ERROR IN BOOK KEEPING\n");
		sleep(3);
	}
}

void print_corr(int dev)
{
	FILE *fp3,*fp5;	
	int id;	
	fp3=fopen(fft_corr,"w");
	fprintf(fp3,"# of avg = %d\n",dev);
	for(id=0;id<R/2;id++)
	{
		if(id % STEP == 0)
		{		
			auto_rho[id]=(double)(avgcorr_rho[id])/dev;
			auto_Q[id]=(double)(avgcorr_Q[id])/dev;
			fprintf(fp3,"%d\t%e\t%e\n",id,auto_rho[id],auto_Q[id]);
		}
	}
	fclose(fp3);
	fp5=fopen(avgvalue,"w");
	fprintf(fp5,"# of avg = %d\n",dev);
	fprintf(fp5,"%e\t%e",avg_rho/(1.0*dev),avg_Q/(1.0*dev));
	fclose(fp5);
	
}

main()
{
	int ms,j,p,l,i,t;
	avg_rho=0;
	avg_Q=0;

	for(i=0;i<R/2;i++)
	{
		avgcorr_rho[i]=0.0;auto_rho[i]=0.0;
		avgcorr_Q[i]=0.0;auto_Q[i]=0.0;
		fft_rho[i]=0.0;
		fft_Q[i]=0.0;
	}
	for(i=R/2;i<R;i++)
	{
		fft_rho[i]=0.0;
		fft_Q[i]=0.0;
	}
	take_input();
	initialize();
	//scanf("%lf",&prob);
	prob=0.9995;
	sprintf(fft_corr,"2DcorrL%dP%5.4lf.dat",L,prob);
	sprintf(Kspace_corr,"K_space_corrL%dP%5.4lf.dat",L,prob);
	fp2=fopen(Kspace_corr,"w");
	sprintf(avgvalue,"avg_valueL%dP%5.4lf.dat",L,prob);
	
	lat_init();
	for(ms=0;ms<T;ms++)
		evolve();
	for(p=1;p<=NUM;p++)
	{
		for(t=0;t<R/2;t++)
		{
			evolve();
			rho[t]=	1.0*(lasth+lastv)/(1.0*N);	
			Q[t]=1.0*(lasth-lastv)/(1.0*(lasth+lastv));
		}
		complex_1d_rhodft();
		complex_1d_orderdft();
		print_corr(p);
	}
	
	for(i=0;i<R;i++)
	{
		fft_rho[i]=(double)(fft_rho[i])/(1.0*NUM);
		fft_Q[i]=(double)(fft_Q[i])/(1.0*NUM);
		fprintf(fp2,"%d\t%e\t%e\n",i,fft_rho[i],fft_Q[i]);
	}

	fclose(fp2);
}


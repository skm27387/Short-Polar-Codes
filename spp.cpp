#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXn	7
#define MAXN	128
#define MAXNL	20000
#define MAXPL	20
#define MAXDBL	((double)(1<<30))
#define MAXNS	(2*MAXN-1)
#define TRUE	1
#define FALSE	0
#define MINSUM_EN

/* Comparator for qsort usage */
struct str
{
    double value;
    int index;
};

int n;
int R;
int N;
int K;
int nL;
int pL;
int f[MAXN];
int P[MAXN];
int w[MAXPL];

/* Decoder Memory */
double LLR[MAXNL][MAXn+1][MAXN]; 	/* beliefs in nL decoders*/
int ucap[MAXNL][MAXn+1][MAXN];		/* decisions in nL decoders */
double PML[MAXNL];					/* Path Metrics */
int vcap[MAXNL][MAXN];
int uhat[MAXNL][MAXN];

int cmp(const void *a, const void *b)
{
    struct str *a1 = (struct str *)a;
    struct str *a2 = (struct str *)b;
    if ((*a1).value > (*a2).value)
        return 1;
    else if ((*a1).value < (*a2).value)
        return -1;
    else
    {
		if ((*a1).index > (*a2).index)
			return 1;
		else
			return -1;
	}
}

/* Arikan's PAC Code */
void initPAC()
{
	n = 7;
	R = 0.5;
	N = (1 << n);
	K = floor(N * R);
	nL = 10000;

	FILE *fp = fopen("fRM_128.txt", "r");
	for (int i = 0; i < N; i++)
	{
		fscanf(fp, "%d", &(f[i]));
	}
	fclose(fp);
	
	pL = 7;
	
	for (int i = 0; i < N; i++)
	{
		P[i] = 1;
	}
	
	w[0] = 1;
	w[1] = 0;
	w[2] = 1;
	w[3] = 1;
	w[4] = 0;
	w[5] = 1;
	w[6] = 1;
}

void initSPC()
{
	n = 7;
	R = 0.5;
	N = (1 << n);
	K = floor(N * R);
	nL = 20000;

	FILE *fp = fopen("fRM_128.txt", "r");
	for (int i = 0; i < N; i++)
	{
		fscanf(fp, "%d", &(f[i]));
	}
	fclose(fp);
	
	pL = 11;
	
	for (int i = 0; i < N; i++)
	{
		if (f[i] == 0)
		{
			P[i] = 1;
		}
	}
	
	w[0] = 1;
	w[1] = 0;
	w[2] = 1;
	w[3] = 1;
	w[4] = 1;
	w[5] = 1;
	w[6] = 0;
	w[7] = 0;
	w[8] = 1;
	w[9] = 1;
	w[10] = 1;
}

void inity(double *y)
{
	for (int i = 0; i < N; i++)
	{
		y[i] = 1;
	}
}

void polarTransform(int *u)
{
	int m = 1; /* Number of bits combined */
	for (int d = (n-1); d >= 0; d--)
	{
		for (int i = 0; i < N;)
		{
			for (int j = i; j < (i + m); j++)
			{
				int a = u[j];
				int b = u[j + m];
				u[j] = a^b;
				u[j + m] = b;
			}
			i += (m << 1);
		}
		m <<= 1;
	}
}

int sum(int *arr, int len)
{
	int s = 0;
	for (int i = 0; i < len; i++)
	{
		s += arr[i];
	}
	return s;
}

void genPolarCwl()
{
	int cww[MAXNL] = {0}; 
	/* Weight of codeword and their frequency is denoted by uwc*/
	int uwc[MAXN + 1] = {0};
	
	for (int i = 0; i < nL; i++)
	{
		int *u = &(uhat[i][0]);
		polarTransform(u);
		cww[i] = sum(u, N);
		uwc[cww[i]]++;
	}
	
	for (int i = 0; i <= N; i++)
	{
		if (uwc[i] > 0)
		{
			printf("Weight %d Frequency %d\n", i, uwc[i]);
		}
	}
}

void initDecoderMemory()
{
	for (int i = 0; i < MAXNL; i++)
	{
		PML[i] = MAXDBL;
		for (int j = 0; j <= MAXn; j++)
		{
			for (int k = 0; k < MAXN; k++)
			{
				LLR[i][j][k] = 0;
				ucap[i][j][k] = 0;
			}
		}
	}
	
	for (int i = 0; i < MAXNL; i++)
	{
		for (int k = 0; k < MAXN; k++)
		{
			uhat[i][k] = 0;
			vcap[i][k] = 0;
		}
	}
}

void initSPCDecodeLLRPM(double *y)
{
	PML[0] = 0;
	for (int i = 0; i < nL; i++)
	{
		for (int k = 0; k < N; k++)
		{
			LLR[i][0][k] = y[k];
		}
	}
}

double f_func(double a, double b)
{
	#if defined(MINSUM_EN)
	/* Min-Sum method */
	double sign_a = (a >= 0) ? 1 : -1;
	double sign_b = (b >= 0) ? 1 : -1;
	double abs_a = fabs(a);
	double abs_b = fabs(b);
	double minabs = (abs_a <= abs_b) ? abs_a : abs_b;
	double out = sign_a * sign_b * minabs;
	#else
	double out = 2 * atanh(tanh(a/2) * tanh(b/2));
	#endif
	return out;
}

double g_func(double u, double a, double b)
{
	double out = b + ((1 - 2 * u) * a);
	return out;
}

void processLeftChild(int *node, int *depth)
{
	int temp = (1 << (n - *depth));
	/* Find which source nodes to process */
	int spS = temp * *node;
	int epS = temp * (*node + 1);
	/* Find start and end for a and b*/
	int nS = *node;
	int dS = *depth;
	int spA = spS;
	int spB = spA + (temp >> 1);

	/* 
	 * Find which destination nodes to put result
	 * In this case, it is left child 
	 */
	*node = (*node << 1);
	*depth = *depth + 1;
	temp >>= 1;
	int spD = temp * (*node);
	int epD = temp * (*node + 1);
	int dD = *depth;
	
	/* Do it for all current running lists */
	for (int i = 0; i < nL; i++)
	{
		int rA = spA;
		int rB = spB;
		for (int k = spD; k < epD; k++)
		{
			LLR[i][dD][k] = f_func(LLR[i][dS][rA], LLR[i][dS][rB]);
			rA++;
			rB++;
		}
	}
}

void processRightChild(int *node, int *depth)
{
	int temp = (1 << (n - *depth));
	/* Find which source nodes to process */
	int spS = temp * *node;
	int epS = temp * (*node + 1);
	/* Find start and end for a and b*/
	int nS = *node;
	int dS = *depth;
	int spA = spS;
	int spB = spA + (temp >> 1);
	
	/* 
	 * For g operation we need decision of bit taken by left
	 * child already
	 */
	int lnode = *node << 1;
	int ldepth = *depth + 1;
	int ltemp = temp >> 1;
	int spLc = ltemp * lnode;
	int epLc = ltemp * (lnode + 1);
	
	/* 
	 * Find which destination nodes to put result
	 * In this case, it is right child 
	 */
	*node = (*node << 1) + 1;
	*depth = *depth + 1;
	temp >>= 1;
	int spD = temp * (*node);
	int epD = temp * (*node + 1);
	int dD = *depth;
	
	/* Do it for all current running lists */
	for (int i = 0; i < nL; i++)
	{
		int rU = spLc;
		int rA = spA;
		int rB = spB;
		for (int k = spD; k < epD; k++)
		{
			LLR[i][dD][k] = g_func(ucap[i][ldepth][rU], LLR[i][dS][rA], LLR[i][dS][rB]);
			rU++;
			rA++;
			rB++;
		}
	}
}

void getPartialSums(int *node, int *depth)
{
	int temp = (1 << (n - *depth));
	int lnode = *node << 1;
	int rnode = lnode + 1;
	int cdepth = *depth + 1;
	int ctemp = temp >> 1;
	
	int spLc = ctemp * lnode;
	int epLc = ctemp * (lnode + 1);
	int spRc = ctemp * rnode;
	int epRc = ctemp * (rnode + 1);
	
	int spD = temp * (*node);
	int epD = temp * (*node + 1);
	
	/* Do it for all the lists */
	for (int i = 0; i < nL; i++)
	{
		int rLc = spLc;
		int rRc = spRc;
		for (int k = spD; k < (spD + ctemp); k++)
		{
			ucap[i][*depth][k] = ucap[i][cdepth][rLc] ^ ucap[i][cdepth][rRc];
			ucap[i][*depth][k + ctemp] = ucap[i][cdepth][rRc];
			rLc++;
			rRc++;
		}
	}
	
	*node = *node >> 1;
	*depth = *depth - 1;
}

int precode(int *v, int node)
{
	if (P[node] == 1)
	{
		int sp = node;
		int ep = sp - pL + 1;
		if (ep < 0)
		{
			ep = 0;
		}
		
		/* Now precode using w */
		int s = 0;
		int wIndex = 0;
		for (int i = sp; i >= ep; i--)
		{
			s = s ^ (w[wIndex] & v[i]);
			wIndex++;
		}
		return s;
	}
	else
	{
		return v[node];
	}
}

double calcPhi(double pm_old, double l, double u)
{
	double pm = pm_old;
	int sign_l = (l >= 0) ? 1 : -1;
	if (u != (0.5 * (1 - sign_l)))
	{
		pm = pm + fabs(l);
	}
	return pm;
}

void doFrozen(int node)
{
	/* For all lists, update the path metric, vcap, ucap */
	for (int i = 0; i < nL; i++)
	{
		double dm = LLR[i][n][node];
		vcap[i][node] = 0;
		ucap[i][n][node] = precode(&(vcap[i][0]), node);
		PML[i] = calcPhi(PML[i], dm, ucap[i][n][node]);
	}
}

/* Below data structures are used for duplication */
double LLR2[2 * MAXNL][MAXn+1][MAXN];
int ucap2[2 * MAXNL][MAXn+1][MAXN];
int vcap2[2 * MAXNL][MAXN];
double PM2[2 * MAXNL] = {0};

void dataDuplicate(int node)
{
	/* Duplication function */
	for (int i = 0; i < nL; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			for (int k = 0; k < N; k++)
			{
				LLR2[i][j][k] = LLR[i][j][k];
				LLR2[i + nL][j][k] = LLR[i][j][k];
				ucap2[i][j][k] = ucap[i][j][k];
				ucap2[i + nL][j][k] = ucap[i][j][k];
			}
		}
		
		for (int k = 0; k < N; k++)
		{
			vcap2[i][k] = vcap[i][k];
			vcap2[i + nL][k] = vcap[i][k];
		}
		
		PM2[i] = PML[i];
		PM2[i + nL] = PML[i];
	}
}

void calcDuplicatePMs(int node)
{
	for (int i = 0; i < nL; i++)
	{
		double dm = LLR[i][n][node];
		vcap2[i][node] = 0;
		vcap2[i + nL][node] = 1;
		ucap2[i][n][node] = precode(&(vcap2[i][0]), node);
		ucap2[i + nL][n][node] = precode(&(vcap2[i + nL][0]), node);
		PM2[i] = calcPhi(PM2[i], dm, ucap2[i][n][node]);
		PM2[i + nL] = calcPhi(PM2[i + nL], dm, ucap2[i + nL][n][node]);
	}
}

struct str pathMetrics[2 * MAXNL];
void sortPM(double *PM2)
{
	for (int i = 0; i < (nL << 1); i++)
	{
		pathMetrics[i].value = PM2[i];
		pathMetrics[i].index = i;
	}
	qsort(pathMetrics, (nL << 1), sizeof(pathMetrics[0]), cmp);
}

void doNonFrozen(int node)
{
	/* Data duplication */
	dataDuplicate(node);
	/* Calculate duplicate path metrics */
	calcDuplicatePMs(node);	
	/* Sorting of Path Metrics */
	sortPM(PM2);
	/* Retain nL paths only */
	for (int i = 0; i < nL; i++)
	{
		PML[i] = pathMetrics[i].value;
		int sI = pathMetrics[i].index;
		for (int j = 0; j <= n; j++)
		{
			for (int k = 0; k < N; k++)
			{
				LLR[i][j][k] = LLR2[sI][j][k];
				ucap[i][j][k] = ucap2[sI][j][k];
			}
		}
		
		for(int k = 0; k < N; k++)
		{
			vcap[i][k] = vcap2[sI][k];
		}
	}
}

void selectUHat()
{
	for (int i = 0; i < nL; i++)
	{
		for (int k = 0; k < N; k++)
		{
			uhat[i][k] = ucap[i][n][k];
		}
	}
}

void SPCDecodeMain()
{
	/* Define node and depth */
	int node = 0;
	int depth = 0;
	/* Node State Vector 0 : L, 1 : R, 2 : U*/
	int npos = 0;
	int ns[MAXNS] = {0};
	
	/* The loop runs till all bits are decoded */
	int done = FALSE;
	while (done == FALSE)
	{
		/* depth = n means we have reached a leaf */
		if (depth == n)
		{
			if (f[node] == 0)
			{
				/* We know that a zero can bring us to this state */
				doFrozen(node);
			}
			else
			{
				doNonFrozen(node);
			}
			
			/* Decide where to go next */
			if (node == (N - 1))
			{
				done = TRUE;
			}
			else
			{
				node >>= 1;
				depth--;
			}
		}
		else /* We are at an intermediate stage */
		{
			/* Position of the node in the node state vector */
			npos = ((1 << depth) - 1) + node;
			
			/* Step L and go to left child */
			if (ns[npos] == 0)
			{	
				processLeftChild(&node, &depth);
				ns[npos] = 1;
			}
			else
			{
				/* Step R and go to right child */
				if (ns[npos] == 1)
				{
					processRightChild(&node, &depth);
					ns[npos] = 2;
				}
				/* Step U and go to parent */
				else
				{
					getPartialSums(&node, &depth);
				}
			}
		}
	}
	
	/* Select uhat */
	for (int i = 0; i < nL; i++)
	{
		for (int k = 0; k < N; k++)
		{
			uhat[i][k] = ucap[i][n][k];
		}
	}
}

int uhat_test[MAXNL][MAXN];
void genPolarCwl2()
{
	int cww[MAXNL] = {0}; 
	/* Weight of codeword and their frequency is denoted by uwc*/
	int uwc[MAXN + 1] = {0};
	
	for (int i = 0; i < nL; i++)
	{
		int *u = &(uhat_test[i][0]);
		polarTransform(u);
		cww[i] = sum(u, N);
		uwc[cww[i]]++;
	}
	
	for (int i = 0; i <= N; i++)
	{
		if (uwc[i] > 0)
		{
			printf("Weight %d Frequency %d\n", i, uwc[i]);
		}
	}
}

void genUhatTest()
{
	for (int i = 0; i < nL; i++)
	{
		for (int k = 0; k < nL; k++)
		{
			uhat_test[i][k] = precode(&(vcap[i][0]), k);
		}
	}
}

void WeightAnalysis()
{
	double y[MAXN] = {0};
	//initPAC();
	initSPC();
	inity(y);
	initDecoderMemory();
	initSPCDecodeLLRPM(y);	
	SPCDecodeMain();
	genPolarCwl();
}

/* Main test function */
int main()
{
	WeightAnalysis();
	return 0;
}

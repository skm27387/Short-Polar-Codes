#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXn	7
#define MAXN	128
#define MAXNL	20000
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
double R;
int N;
int K;
int nL;
int pL;
int f[MAXN];

/* Decoder Memory */
double LLR[MAXNL][MAXn+1][MAXN]; 	/* beliefs in nL decoders*/
int ucap[MAXNL][MAXn+1][MAXN];		/* decisions in nL decoders */
double PML[MAXNL];					/* Path Metrics */

/* Below data structures are used for duplication */
double LLR2[2 * MAXNL][MAXn+1][MAXN];
int ucap2[2 * MAXNL][MAXn+1][MAXN];
int vcap2[2 * MAXNL][MAXN];
double PM2[2 * MAXNL] = {0};
int uhat[MAXNL][MAXN];

int nrOrder[1024] = {
0, 1, 2, 4, 8, 16, 32, 3, 5, 64, 9, 6, 17, 10, 18, 128, 12, 33, 65, 20, 256, 34, 24, 36, 7, 129, 66, 512, 11, 40, 68, 130,
19, 13, 48, 14, 72, 257, 21, 132, 35, 258, 26, 513, 80, 37, 25, 22, 136, 260, 264, 38, 514, 96, 67, 41, 144, 28, 69, 42,
516, 49, 74, 272, 160, 520, 288, 528, 192, 544, 70, 44, 131, 81, 50, 73, 15, 320, 133, 52, 23, 134, 384, 76, 137, 82, 56, 27,
97, 39, 259, 84, 138, 145, 261, 29, 43, 98, 515, 88, 140, 30, 146, 71, 262, 265, 161, 576, 45, 100, 640, 51, 148, 46, 75, 266, 273, 517, 104, 162,
53, 193, 152, 77, 164, 768, 268, 274, 518, 54, 83, 57, 521, 112, 135, 78, 289, 194, 85, 276, 522, 58, 168, 139, 99, 86, 60, 280, 89, 290, 529, 524,
196, 141, 101, 147, 176, 142, 530, 321, 31, 200, 90, 545, 292, 322, 532, 263, 149, 102, 105, 304, 296, 163, 92, 47, 267, 385, 546, 324, 208, 386, 150, 153,
165, 106, 55, 328, 536, 577, 548, 113, 154, 79, 269, 108, 578, 224, 166, 519, 552, 195, 270, 641, 523, 275, 580, 291, 59, 169, 560, 114, 277, 156, 87, 197,
116, 170, 61, 531, 525, 642, 281, 278, 526, 177, 293, 388, 91, 584, 769, 198, 172, 120, 201, 336, 62, 282, 143, 103, 178, 294, 93, 644, 202, 592, 323, 392,
297, 770, 107, 180, 151, 209, 284, 648, 94, 204, 298, 400, 608, 352, 325, 533, 155, 210, 305, 547, 300, 109, 184, 534, 537, 115, 167, 225, 326, 306, 772, 157,
656, 329, 110, 117, 212, 171, 776, 330, 226, 549, 538, 387, 308, 216, 416, 271, 279, 158, 337, 550, 672, 118, 332, 579, 540, 389, 173, 121, 553, 199, 784, 179,
228, 338, 312, 704, 390, 174, 554, 581, 393, 283, 122, 448, 353, 561, 203, 63, 340, 394, 527, 582, 556, 181, 295, 285, 232, 124, 205, 182, 643, 562, 286, 585,
299, 354, 211, 401, 185, 396, 344, 586, 645, 593, 535, 240, 206, 95, 327, 564, 800, 402, 356, 307, 301, 417, 213, 568, 832, 588, 186, 646, 404, 227, 896, 594,
418, 302, 649, 771, 360, 539, 111, 331, 214, 309, 188, 449, 217, 408, 609, 596, 551, 650, 229, 159, 420, 310, 541, 773, 610, 657, 333, 119, 600, 339, 218, 368,
652, 230, 391, 313, 450, 542, 334, 233, 555, 774, 175, 123, 658, 612, 341, 777, 220, 314, 424, 395, 673, 583, 355, 287, 183, 234, 125, 557, 660, 616, 342, 316,
241, 778, 563, 345, 452, 397, 403, 207, 674, 558, 785, 432, 357, 187, 236, 664, 624, 587, 780, 705, 126, 242, 565, 398, 346, 456, 358, 405, 303, 569, 244, 595,
189, 566, 676, 361, 706, 589, 215, 786, 647, 348, 419, 406, 464, 680, 801, 362, 590, 409, 570, 788, 597, 572, 219, 311, 708, 598, 601, 651, 421, 792, 802, 611,
602, 410, 231, 688, 653, 248, 369, 190, 364, 654, 659, 335, 480, 315, 221, 370, 613, 422, 425, 451, 614, 543, 235, 412, 343, 372, 775, 317, 222, 426, 453, 237,
559, 833, 804, 712, 834, 661, 808, 779, 617, 604, 433, 720, 816, 836, 347, 897, 243, 662, 454, 318, 675, 618, 898, 781, 376, 428, 665, 736, 567, 840, 625, 238,
359, 457, 399, 787, 591, 678, 434, 677, 349, 245, 458, 666, 620, 363, 127, 191, 782, 407, 436, 626, 571, 465, 681, 246, 707, 350, 599, 668, 790, 460, 249, 682,
573, 411, 803, 789, 709, 365, 440, 628, 689, 374, 423, 466, 793, 250, 371, 481, 574, 413, 603, 366, 468, 655, 900, 805, 615, 684, 710, 429, 794, 252, 373, 605,
848, 690, 713, 632, 482, 806, 427, 904, 414, 223, 663, 692, 835, 619, 472, 455, 796, 809, 714, 721, 837, 716, 864, 810, 606, 912, 722, 696, 377, 435, 817, 319,
621, 812, 484, 430, 838, 667, 488, 239, 378, 459, 622, 627, 437, 380, 818, 461, 496, 669, 679, 724, 841, 629, 351, 467, 438, 737, 251, 462, 442, 441, 469, 247,
683, 842, 738, 899, 670, 783, 849, 820, 728, 928, 791, 367, 901, 630, 685, 844, 633, 711, 253, 691, 824, 902, 686, 740, 850, 375, 444, 470, 483, 415, 485, 905,
795, 473, 634, 744, 852, 960, 865, 693, 797, 906, 715, 807, 474, 636, 694, 254, 717, 575, 913, 798, 811, 379, 697, 431, 607, 489, 866, 723, 486, 908, 718, 813,
476, 856, 839, 725, 698, 914, 752, 868, 819, 814, 439, 929, 490, 623, 671, 739, 916, 463, 843, 381, 497, 930, 821, 726, 961, 872, 492, 631, 729, 700, 443, 741,
845, 920, 382, 822, 851, 730, 498, 880, 742, 445, 471, 635, 932, 687, 903, 825, 500, 846, 745, 826, 732, 446, 962, 936, 475, 853, 867, 637, 907, 487, 695, 746,
828, 753, 854, 857, 504, 799, 255, 964, 909, 719, 477, 915, 638, 748, 944, 869, 491, 699, 754, 858, 478, 968, 383, 910, 815, 976, 870, 917, 727, 493, 873, 701,
931, 756, 860, 499, 731, 823, 922, 874, 918, 502, 933, 743, 760, 881, 494, 702, 921, 501, 876, 847, 992, 447, 733, 827, 934, 882, 937, 963, 747, 505, 855, 924,
734, 829, 965, 938, 884, 506, 749, 945, 966, 755, 859, 940, 830, 911, 871, 639, 888, 479, 946, 750, 969, 508, 861, 757, 970, 919, 875, 862, 758, 948, 977, 923,
972, 761, 877, 952, 495, 703, 935, 978, 883, 762, 503, 925, 878, 735, 993, 885, 939, 994, 980, 926, 764, 941, 967, 886, 831, 947, 507, 889, 984, 751, 942, 996,
971, 890, 509, 949, 973, 1000, 892, 950, 863, 759, 1008, 510, 979, 953, 763, 974, 954, 879, 981, 982, 927, 995, 765, 956, 887, 985, 997, 986, 943, 891, 998, 766,
511, 988, 1001, 951, 1002, 893, 975, 894, 1009, 955, 1004, 1010, 957, 983, 958, 987, 1012, 999, 1016, 767, 989, 1003, 990, 1005, 959, 1011, 1013, 895, 1006, 1014, 1017, 1018,
991, 1020, 1007, 1015, 1019, 1021, 1022, 1023
};

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

void initPolar()
{
	n = 7;
	R = 0.5;
	N = (1 << n);
	K = (int)floor((double)N * R);
	nL = 2;

	/* Get the polar frozen array */
	for (int i = 0; i < N; i++)
	{
		f[i] = 0;
	}
	int j = 0;
	int Q1[MAXN] = {0};
	for (int i = 0; i < 1024; i++)
	{
		if (nrOrder[i] < N)
		{
			Q1[j] = nrOrder[i];
			j++;
			if (j == N)
			{
				break;
			}
		}
	}
	for (int i = N-K; i < N; i++)
	{
		f[Q1[i]] = 1;
	}
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
	int minHammingWeight = N;
	for (int i = 0; i < nL; i++)
	{
		int *u = &(uhat[i][0]);
		polarTransform(u);
		cww[i] = sum(u, N);
		if (cww[i] < N)
		{
			minHammingWeight = cww[i];
		}
	}
	printf("The minimum hamming weight of this polar code is %d\n", minHammingWeight);
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
		ucap[i][n][node] = 0;
		PML[i] = calcPhi(PML[i], dm, ucap[i][n][node]);
	}
}

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
		
		PM2[i] = PML[i];
		PM2[i + nL] = PML[i];
	}
}

void calcDuplicatePMs(int node)
{
	for (int i = 0; i < nL; i++)
	{
		double dm = LLR[i][n][node];
		ucap2[i][n][node] = 0;
		ucap2[i + nL][n][node] = 1;
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

void WeightAnalysis()
{
	double y[MAXN] = {0};
	initPolar();
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#define euclidDist(X1, Y1, X2, Y2) (sqrt((((X1)-(X2))*((X1)-(X2))+((Y1)-(Y2))*((Y1)-(Y2))))) // Calculates the Euclidean distance between two points
// Input	:	X1 (X-coordinate of point #1)
//				Y1 (Y-coordinate of point #1)
//				X2 (X-coordinate of point #2)
//				Y2 (Y-coordinate of point #2)
// Output	:	(Euclidean distance between point #1 and #2)

double bruteClosestPair(
	unsigned L, unsigned R, // Current leftmost and rightmost indices
	double* X, double* Y, // Arrays storing the X and Y coordinates of each point ordered by index
	unsigned* Xid, // Array storing the indices of the points in the non-decreasing-order of X-coordinates
	unsigned* Yid, // Initially empty array used to store the indices of the points in the non-increasing-order of Y-coordinates
	unsigned* pt1, unsigned* pt2 // Indices of closest pair(pt1, pt2) of points
);
//	Input:	L, R, X[], Y[], Xid[], Yid[]
//	Output:	pt1, pt2( the distance between (X[pt1],Y[pt1])
//					and (X[pt2],Y[pt2]) is the closest ),
//			Euclidean distance between pt1 and pt2
void mergeSort(
	unsigned L, unsigned R, // Current leftmost and rightmost indices
	double* X, // Array storing the X-coordinates of each point ordered by index
	unsigned* Xid,	//	Array initially storing the indices of points in increasing order to be sorted in
					//	non-decreasing-order of the corresponding X-coordinates
	unsigned* TMP // Temporary array to store the sorted indices
);
//	Input:	L, R, X[], Xid[], TMP[]
//	Output:	Xid[L...R] sorted in non-decreasing order
void mergeXid(
	unsigned L, unsigned M, unsigned R, // Current leftmost, middle, and rightmost indices
	double* X,	//	Array storing the X-coordinates of each point ordered by index
	unsigned* Xid, // Arrays storing the partially sorted array of X indices in non-decreasing-order of X-coordinates
	unsigned* TMP // Temporary array to store the sorted indices
);
//	Input:	L, M, R, X[], Xid[], TMP[]
//	Output:	Xid[L...M] and Xid[M+1...R] merged in non-decreasing order
void mergeYid(
	unsigned L, unsigned M, unsigned R, // Current leftmost, middle, and rightmost indices
	double* Y, // Array storing the Y-coordinates of each point ordered by index
	unsigned* Yid, // Array storing the partially sorted array of Y indices in non-decreasing-order of Y-coordinates
	unsigned* TMP // Temporary array to store the sorted indices
);
//	Input:	L, M, R, Y[], Yid[], TMP[]
//	Output:	Yid[L...M] and Yid[M+1...R] merged in non-increasing order
double combine(
	unsigned L, unsigned M, unsigned R, // Current leftmost, middle, and rightmost indices
	double dLR, // Distance between the closest pair of points in the left or the right divisions
	double* X, double* Y, // Arrays storing the X and Y coordinates of each point ordered by index
	unsigned* Xid, // Array storing the indices of the points in the non-decreasing-order of X-coordinates
	unsigned* pt1, unsigned* pt2, // Indices of closest pair(pt1, pt2) of points
	unsigned* TMP // Temporary array to store the sorted indices
);
//	Input:	L, M, R, dLR, X[], Y[], Xid[], TMP[]
//	Output:	pt1, pt2( the distance between (X[pt1],Y[pt1])
//					and (X[pt2],Y[pt2]) is the closest ),
//			Euclidean distance between pt1 and pt2

void sortXid(double* X, unsigned* Xid, unsigned* TMP, unsigned N) {
	if(N == 0) return;
	mergeSort(0, N-1, X, Xid, TMP);
}

double closestPairDC(unsigned L, unsigned R, unsigned* pt1, unsigned* pt2, double* X, double* Y, unsigned* Xid, unsigned* Yid, unsigned* TMP, unsigned THR) {
	double dMin, dL, dR;
	unsigned M, ptL1, ptL2, ptR1, ptR2;
	*pt1 = *pt2 = 0; // Initialize pt1 and pt2 to 0 to account for a single point input
	if (R - L + 1 <= THR)
		// Use brute-force algorithm for number points under threshold
		dMin = bruteClosestPair(L, R, X, Y, Xid, Yid, pt1, pt2);
	else {
		M = (L + R) / 2;

		// Apply D&C algorithm
		dL = closestPairDC(L, M, &ptL1, &ptL2, X, Y, Xid, Yid, TMP, THR);
		dR = closestPairDC(M + 1, R, &ptR1, &ptR2, X, Y, Xid, Yid, TMP, THR);
		if (dL < dR) {
			dMin = dL;
			*pt1 = ptL1;
			*pt2 = ptL2;
		}
		else {
			dMin = dR;
			*pt1 = ptR1;
			*pt2 = ptR2;
		}

		mergeYid(L, M, R, Y, Yid, TMP); // Merge Yid[L...M] and Yid[M+1...R] in non-increasing order of corresponding Y

		dMin = combine(L, M, R, dMin, X, Y, Xid, pt1, pt2, TMP); // This step finds the closest pair between the points of left and right division 
	}

	return dMin;
}

double bruteClosestPair(unsigned L, unsigned R, double* X, double* Y, unsigned* Xid, unsigned* Yid, unsigned* pt1, unsigned* pt2) {
	unsigned i, j, idxTemp, ptsToComp;
	double dMin, dTemp;

	ptsToComp = R - L + 1; // Number of points to compare

	for (i = L; i <= R; i++) {
		Yid[i] = Xid[i];
	}
	//Bubble sort Yid[L...R] in non-increasing order
	for (i = 0; i < ptsToComp - 1; i++) {
		for (j = L+i; j < R; j++) {
			if (Y[Yid[j]] < Y[Yid[j+1]]) {
				idxTemp = Yid[j];
				Yid[j] = Yid[j+1];
				Yid[j+1] = idxTemp;
			}
		}
	}

	dMin = DBL_MAX;
	// Using the sorted Yid[L...R], compare every pair of points in Yid[L...R]
	for (i = L+1; i <= R; i++) {
		for (j = i; j <= R; j++) {
			if (Y[Yid[i-1]] - Y[Yid[j]] >= dMin) break; //	When this condition is met, it is unnecessary to compare more points
														//	of the lower Y-coordinates
			dTemp = euclidDist(X[Yid[i-1]], Y[Yid[i-1]], X[Yid[j]], Y[Yid[j]]);
			if (dTemp < dMin) {
				dMin = dTemp;
				// Assign new closest pair of points
				*pt1 = Yid[j];
				*pt2 = Yid[i-1];
			}
		}
	}

	return dMin;
}

void mergeSort(unsigned L, unsigned R, double* X, unsigned* Xid, unsigned* TMP) {
	unsigned M;
	if (L < R) {
		M = (L + R) / 2;
		// Using D&C algorithm
		mergeSort(L, M, X, Xid, TMP);
		mergeSort(M+1, R, X, Xid, TMP);
		mergeXid(L, M, R, X, Xid, TMP);
	}
}

void mergeXid(unsigned L, unsigned M, unsigned R, double* X, unsigned* Xid, unsigned* TMP) {
	unsigned i, j, tmpLen;
	i = L;
	j = M + 1;
	tmpLen = 0;
	// Xid[L...M] and Xid[M+1...R] are each already sorted in non-decreasing order
	// Store all indices of Xid[L...M] and Xid[M+1...R] in TMP array in non-decreasing order
	while (i <= M && j <= R) {
		if (X[Xid[i]] < X[Xid[j]])
			TMP[tmpLen++] = Xid[i++];
		else
			TMP[tmpLen++] = Xid[j++];
	}
	// Store the leftover indices
	if (i > M) {
		while (j <= R) {
			TMP[tmpLen++] = Xid[j++];
		}
	}
	else {
		while (i <= M) {
			TMP[tmpLen++] = Xid[i++];
		}
	}
	// Replace the indices of the array Xid[L...R] with sorted indices in TMP
	for (i = L, j = 0; j < tmpLen; i++, j++) {
		Xid[i] = TMP[j];
	}
}

void mergeYid(unsigned L, unsigned M, unsigned R, double* Y, unsigned* Yid, unsigned* TMP) {
	unsigned i, j, tmpLen;
	i = L;
	j = M + 1;
	tmpLen = 0;
	// Yid[L...M] and Yid[M+1...R] are each already sorted in non-increasing order
	// Store all indices of Yid[L...M] and Yid[M+1...R] in TMP array in non-increasing order
	while (i <= M && j <= R) {
		if (Y[Yid[i]] > Y[Yid[j]])
			TMP[tmpLen++] = Yid[i++];
		else
			TMP[tmpLen++] = Yid[j++];
	}
	// Store leftover indices
	if (i > M) {
		while (j <= R) {
			TMP[tmpLen++] = Yid[j++];
		}
	}
	else {
		while (i <= M) {
			TMP[tmpLen++] = Yid[i++];
		}
	}
	// Replace the indices of the array Yid[L...R] with sorted indices in TMP
	for (i = L, j = 0; j < tmpLen; i++, j++) {
		Yid[i] = TMP[j];
	}
}

double combine(unsigned L, unsigned M, unsigned R, double dLR, double* X, double* Y, unsigned* Xid, unsigned* pt1, unsigned* pt2, unsigned* TMP) {
	unsigned i, j, tmpLen, idxTemp, rightLend, leftRend;
	double dTemp;
	rightLend = M+1;
	leftRend = M;
	i = M;
	j = M + 1;
	tmpLen = 0;

	// Store points with its X-coordinates that satisfy
	// X[Xid[M+1]] <= X[p] <= X[Xid[M+1]] + dLR
	while (X[Xid[rightLend]] + dLR >= X[Xid[j]]) {
		TMP[tmpLen++] = Xid[j];
		if (j == R) break; // When index reaches the rightmost end, exit loop
		j++;
	}
	// Store points with its X-coordinates that satisfy
	// X[Xid[M]] - dLR <= X[p] <= X[Xid[M]]
	while (X[Xid[leftRend]] - dLR <= X[Xid[i]]) {
		TMP[tmpLen++] = Xid[i];
		if (i == L) break; // When index reaches the leftmost end, exit loop
		i--;
	}
	//Sort indices in TMP in the order of non-increasing order of Y-coordinates
	for (i = 0; i < tmpLen-1; i++) {
		for (j = 0; j < tmpLen-i-1; j++) {
			if (Y[TMP[j]] < Y[TMP[j + 1]]) {
				idxTemp = TMP[j];
				TMP[j] = TMP[j + 1];
				TMP[j + 1] = idxTemp;
			}
		}
	}
	// Find the closest pair of points arount the center of Xid[L...R]
	for (i = 0; i < tmpLen-1; i++) {
		for (j = i + 1; j < tmpLen; j++) {
			if (Y[TMP[i]] - Y[TMP[j]] >= dLR) break;	// When this condition is met, it is unnecessary to compare
														// to the points with smaller Y-coordinates
			dTemp = euclidDist(X[TMP[i]], Y[TMP[i]], X[TMP[j]], Y[TMP[j]]);
			if (dTemp < dLR) {
				dLR = dTemp;
				*pt1 = TMP[j];
				*pt2 = TMP[i];
			}
		}
	}
	return dLR;
}
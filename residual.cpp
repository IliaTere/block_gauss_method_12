#include "functions.hpp"

void add(double *a, double *b, int n, int m)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            a[i*m+j] = fabs(a[i*m+j]) + fabs(b[i*m+j]);
}

double find_residual(double *a, double *b, double *block, double *sum_array, int n, int m)
{
	int i = 0, p = 0, q = 0;
	int v = 0, h = 0;
	int r = 0, t = 0;
	int s = 0;
	int ah = 0;
	int z = 0;
	int v3 = 0, h3 = 0;
	double res = 0;
	double s00 = 0, s01 = 0, s02 = 0, s10 = 0, s20 = 0, s11 = 0, s12 = 0,
				 s21 = 0, s22 = 0, sum = 0;
	double *pa = nullptr, *pb = nullptr;
	int k = n / m;
	int l = n % m;
	int bl = (l != 0 ? k + 1 : k);

	for (int j = 0; j < bl; j++)
	{
		h = (j < k ? m : l);
		for (p = 0; p < h; p++)
			sum_array[p] = 0;

		for (i = 0; i < bl; i++)
		{
			v = (i < k ? m : l);
			for (p = 0; p < m * m; p++)
				block[p] = 0;

			for (s = 0; s < bl; s++)
			{
				ah = (s < k) ? m : l;
				pa = a + i * n * m + s * m;
				pb = b + s * n * m + j * m;
				v3 = v % 3;
				h3 = h % 3;
				for (r = 0; r < v3; r++)
				{
					for (t = 0; t < h3; t++)
					{
						sum = 0;
						for (q = 0; q < ah; q++)
							sum += pa[r * n + q] * pb[q * n + t];
						block[r * m + t] += sum;
					}
					for (; t < h; t += 3)
					{
						s00 = 0;
						s01 = 0;
						s02 = 0;
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s01 += pa[r * n + q] * pb[q * n + t + 1];
							s02 += pa[r * n + q] * pb[q * n + t + 2];
						}
						block[r * m + t] += s00;
						block[r * m + t + 1] += s01;
						block[r * m + t + 2] += s02;
					}
				}
				for (; r < v; r += 3)
				{
					for (t = 0; t < h3; t++)
					{
						s00 = 0;
						s10 = 0;
						s20 = 0;
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s10 += pa[(r + 1) * n + q] * pb[q * n + t];
							s20 += pa[(r + 2) * n + q] * pb[q * n + t];
						}
						block[r * m + t] += s00;
						block[(r + 1) * m + t] += s10;
						block[(r + 2) * m + t] += s20;
					}
					for (; t < h; t += 3)
					{
						s00 = 0;
						s01 = 0;
						s02 = 0;
						s10 = 0;
						s11 = 0;
						s12 = 0;
						s20 = 0;
						s21 = 0;
						s22 = 0;
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s01 += pa[r * n + q] * pb[q * n + t + 1];
							s02 += pa[r * n + q] * pb[q * n + t + 2];
							s10 += pa[(r + 1) * n + q] * pb[q * n + t];
							s11 += pa[(r + 1) * n + q] * pb[q * n + t + 1];
							s12 += pa[(r + 1) * n + q] * pb[q * n + t + 2];
							s20 += pa[(r + 2) * n + q] * pb[q * n + t];
							s21 += pa[(r + 2) * n + q] * pb[q * n + t + 1];
							s22 += pa[(r + 2) * n + q] * pb[q * n + t + 2];
						}
						block[r * m + t] += s00;
						block[r * m + t + 1] += s01;
						block[r * m + t + 2] += s02;
						block[(r + 1) * m + t] += s10;
						block[(r + 1) * m + t + 1] += s11;
						block[(r + 1) * m + t + 2] += s12;
						block[(r + 2) * m + t] += s20;
						block[(r + 2) * m + t + 1] += s21;
						block[(r + 2) * m + t + 2] += s22;
					}
				}
			}
			if (i == j)
				for (p = 0; p < v; p++)
					block[p * m + p] -= 1;

			for (p = 0; p < h; p++)
				for (z = 0; z < v; z++)
					sum_array[p] += fabs(block[z * m + p]);
		}

		for (i = 0; i < h; i++)
			if (sum_array[i] > res)
				res = sum_array[i];
	}
	return res;
}


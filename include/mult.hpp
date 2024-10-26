void mult(double *a, double *b, double *res, int m1, int m2, int m3, int m)
{
	// a m1 x m2
	// b m2 x m3
	int t = 0, q = 0, r = 0;
	int v = m1, h = m3, ah = m2;
	int v3 = v % 3, h3 = h % 3;
	double s00 = 0, s01 = 0, s02 = 0;
	double s10 = 0, s11 = 0, s12 = 0;
	double s20 = 0, s21 = 0, s22 = 0;

	for (r = 0; r < v; r++)
		for (t = 0; t < h; t++)
			res[r * m + t] = 0;

	for (r = 0; r < v3; r++)
	{
		for (t = 0; t < h3; t++)
		{
			double sum;
			sum = 0;
			for (q = 0; q < ah; q++)
				sum += a[r * m + q] * b[q * m + t];
			res[r * m + t] += sum;
		}
		for (; t < h; t += 3)
		{
			s00 = 0;
			s01 = 0;
			s02 = 0;
			for (q = 0; q < ah; q++)
			{
				s00 += a[r * m + q] * b[q * m + t];
				s01 += a[r * m + q] * b[q * m + t + 1];
				s02 += a[r * m + q] * b[q * m + t + 2];
			}
			res[r * m + t] += s00;
			res[r * m + t + 1] += s01;
			res[r * m + t + 2] += s02;
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
				s00 += a[r * m + q] * b[q * m + t];
				s10 += a[(r + 1) * m + q] * b[q * m + t];
				s20 += a[(r + 2) * m + q] * b[q * m + t];
			}
			res[r * m + t] += s00;
			res[(r + 1) * m + t] += s10;
			res[(r + 2) * m + t] += s20;
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
				s00 += a[r * m + q] * b[q * m + t];
				s01 += a[r * m + q] * b[q * m + t + 1];
				s02 += a[r * m + q] * b[q * m + t + 2];
				s10 += a[(r + 1) * m + q] * b[q * m + t];
				s11 += a[(r + 1) * m + q] * b[q * m + t + 1];
				s12 += a[(r + 1) * m + q] * b[q * m + t + 2];
				s20 += a[(r + 2) * m + q] * b[q * m + t];
				s21 += a[(r + 2) * m + q] * b[q * m + t + 1];
				s22 += a[(r + 2) * m + q] * b[q * m + t + 2];
			}
			res[r * m + t] += s00;
			res[r * m + t + 1] += s01;
			res[r * m + t + 2] += s02;
			res[(r + 1) * m + t] += s10;
			res[(r + 1) * m + t + 1] += s11;
			res[(r + 1) * m + t + 2] += s12;
			res[(r + 2) * m + t] += s20;
			res[(r + 2) * m + t + 1] += s21;
			res[(r + 2) * m + t + 2] += s22;
		}
	}
}
#define EPS 1e-16
#define CANNOT_SOLVE -1
#define IRREVERSIBLE -2
#define SUCCESS 0

int gauss_classic_row(double *matrix, double *inverse_matrix, int *index, int n, double matrix_norm, int row_ind)
{
	int max_index = 0;
	double tmp_ = 0, max = 0;
	int index_0 = 0;
	int index_1 = 0;
	int swap = 0;

	for (int i = 0; i < n; i++) // Присоединенная матрица
	{
		index_0 = i * row_ind;
		for (int j = 0; j < n; j++)
			i == j ? inverse_matrix[index_0 + j] = 1. : inverse_matrix[index_0 + j] = 0.;
	}

	for (int i = 0; i < n; i++)
		index[i] = i;

	for (int i = 0; i < n; i++) // Прямой ход метода Гаусса
	{
		index_0 = i * row_ind;
		max = fabs(matrix[index_0 + i]); // нахождение главного элемента
		max_index = i;

		for (int j = i + 1; j < n; j++)
			if (max < fabs(matrix[index_0 + j]))
			{
				max = fabs(matrix[index_0 + j]);
				max_index = j;
			}

		swap = index[i];
		index[i] = index[max_index];
		index[max_index] = swap;

		for (int j = 0; j < n; j++) // переставляем столбцы местами
		{
			index_1 = j * row_ind;
			tmp_ = matrix[index_1 + i];
			matrix[index_1 + i] = matrix[index_1 + max_index];
			matrix[index_1 + max_index] = tmp_;
		}

		if (fabs(matrix[index_0 + i]) < matrix_norm * EPS) // если элемент нулевой, метод неприменим
		{
			return IRREVERSIBLE;
		}

		// tmp_ = matrix[index_0 + i];
		tmp_ = 1 / matrix[index_0 + i];
		matrix[index_0 + i] = 1.0;
		for (int j = i + 1; j < n; j++) // умножаем i строку на обратный элемент
			matrix[index_0 + j] *= tmp_;

		for (int j = 0; j < n; j++) // присоединенная матрица. умножаем i строку на обратный элемент
			inverse_matrix[index_0 + j] *= tmp_;

		for (int j = i + 1; j < n; j++)
		{
			tmp_ = matrix[j * row_ind + i];
			for (int k = i; k < n; k++) // вычитание строки
				matrix[j * row_ind + k] -= matrix[index_0 + k] * tmp_;
			for (int k = 0; k < n; k++)
				inverse_matrix[j * row_ind + k] -= inverse_matrix[index_0 + k] * tmp_; // присоединенная матрица. вычитание строки умноженной на число
		}
	}

	for (int i = 0; i < n; i++) // Обратный ход
		for (int j = n - 1; j >= 0; j--)
		{
			index_0 = j * row_ind;
			tmp_ = inverse_matrix[index_0 + i];
			for (int k = j + 1; k < n; k++)
				tmp_ -= matrix[index_0 + k] * inverse_matrix[k * row_ind + i];
			inverse_matrix[index_0 + i] = tmp_;
		}

	for (int i = 0; i < n; i++)
	{
		index_0 = i * row_ind;
		index_1 = index[i] * row_ind;
		for (int j = 0; j < n; j++)
			matrix[index_1 + j] = inverse_matrix[index_0 + j];
	}

	for (int i = 0; i < n; i++)
	{
		index_0 = i * row_ind;
		for (int j = 0; j < n; j++)
			inverse_matrix[index_0 + j] = matrix[index_0 + j];
	}

	return SUCCESS;
}

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
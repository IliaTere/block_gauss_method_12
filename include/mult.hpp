#define EPS 1e-16

int gauss_classic_row(double *matrix, double *inverse_matrix, int *index, int n, double matrix_norm, int row_ind)
{
    int max_col_index = 0;
    double temp = 0, max_abs_value = 0;
    int row_offset = 0;
    int col_offset = 0;
    int temp_index = 0;

    // Инициализация обратной матрицы
    for (int row = 0; row < n; row++)
    {
        row_offset = row * row_ind;
        for (int col = 0; col < n; col++)
            inverse_matrix[row_offset + col] = (row == col) ? 1.0 : 0.0;
    }

    // Инициализация индексов
    for (int i = 0; i < n; i++)
        index[i] = i;

    // Прямой ход
    for (int row = 0; row < n; row++)
    {
        row_offset = row * row_ind;
        max_abs_value = fabs(matrix[row_offset + row]);
        max_col_index = row;
	
        // Поиск максимального элемента в столбце
        for (int col = row + 1; col < n; col++)
            if (max_abs_value < fabs(matrix[row_offset + col]))
            {
                max_abs_value = fabs(matrix[row_offset + col]);
                max_col_index = col;
            }

        // Перестановка индексов
        temp_index = index[row];
        index[row] = index[max_col_index];
        index[max_col_index] = temp_index;

        // Перестановка строк в матрице
        for (int col = 0; col < n; col++)
        {
            col_offset = col * row_ind;
            temp = matrix[col_offset + row];
            matrix[col_offset + row] = matrix[col_offset + max_col_index];
            matrix[col_offset + max_col_index] = temp;
        }

        // Проверка на обратимость
        if (fabs(matrix[row_offset + row]) < matrix_norm * EPS)
        {
            return -2;
        }

        // Нормализация строки
        temp = 1 / matrix[row_offset + row];
        matrix[row_offset + row] = 1.0;
        for (int col = row + 1; col < n; col++)
            matrix[row_offset + col] *= temp;

        for (int col = 0; col < n; col++)
            inverse_matrix[row_offset + col] *= temp;

        // Вычитание строки
        for (int next_row = row + 1; next_row < n; next_row++)
        {
            temp = matrix[next_row * row_ind + row];
            for (int col = row; col < n; col++) // Вычитание строки
                matrix[next_row * row_ind + col] -= matrix[row_offset + col] * temp;
            for (int col = 0; col < n; col++)
                inverse_matrix[next_row * row_ind + col] -= inverse_matrix[row_offset + col] * temp;
        }
    }

    // Обратный ход
    for (int row = 0; row < n; row++) // Обратный ход
        for (int col = n - 1; col >= 0; col--)
        {
            row_offset = col * row_ind;
            temp = inverse_matrix[row_offset + row];
            for (int next_col = col + 1; next_col < n; next_col++)
                temp -= matrix[row_offset + next_col] * inverse_matrix[next_col * row_ind + row];
            inverse_matrix[row_offset + row] = temp;
        }

    for (int row = 0; row < n; row++)
    {
        row_offset = row * row_ind;
        col_offset = index[row] * row_ind;
        for (int col = 0; col < n; col++)
            matrix[col_offset + col] = inverse_matrix[row_offset + col];
    }

    for (int row = 0; row < n; row++)
    {
        row_offset = row * row_ind;
        for (int col = 0; col < n; col++)
            inverse_matrix[row_offset + col] = matrix[row_offset + col];
    }

    return 0;
}
void mult(double *a, double *b, double *res, int m1, int m2, int m3, int m)
{
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

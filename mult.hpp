#define EPS 1e-16
#define UNUSED(x) (void)(x)
int inverse_matrix(double *matrix, double *inverse_matrix, int *index, int n, double matrix_norm, int row_ind)
{
    int max_col_index = 0;
    double temp = 0, max_abs_value = 0;
    int row_offset = 0;
    int col_offset = 0;
    int temp_index = 0;
    int i, row, col, next_row = 0;
    for (row = 0; row < n; row++)
    {
        row_offset = row * row_ind;
        for (col = 0; col < n; col++)
            if(row==col)
                inverse_matrix[row_offset + col] = 1.0;
            else
                inverse_matrix[row_offset + col] = 0.0;
    }

    for (i = 0; i < n; i++)
        index[i] = i;

    // Прямой ход
    for (row = 0; row < n; row++)
    {
        row_offset = row * row_ind;
        max_abs_value = fabs(matrix[row_offset + row]);
        max_col_index = row;
	
        // Поиск максимального элемента в столбце
        for (col = row + 1; col < n; col++)
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
        for (col = 0; col < n; col++)
        {
            col_offset = col * row_ind;
            temp = matrix[col_offset + row];
            matrix[col_offset + row] = matrix[col_offset + max_col_index];
            matrix[col_offset + max_col_index] = temp;
        }

        // Проверка на обратимость
        if (fabs(matrix[row_offset + row]) < matrix_norm * EPS)
        {
            return -1;
        }

        // Нормализация строки
        temp = 1 / matrix[row_offset + row];
        matrix[row_offset + row] = 1.0;
        for (col = row + 1; col < n; col++)
            matrix[row_offset + col] = matrix[row_offset + col] * temp;

        for (col = 0; col < n; col++)
            inverse_matrix[row_offset + col] = inverse_matrix[row_offset + col] * temp;

        // Вычитание строки
        for (next_row = row + 1; next_row < n; next_row++)
        {
            temp = matrix[next_row * row_ind + row];
            for (col = row; col < n; col++) // Вычитание строки
                matrix[next_row * row_ind + col] = matrix[next_row * row_ind + col] - matrix[row_offset + col] * temp;
            for (col = 0; col < n; col++)
                inverse_matrix[next_row * row_ind + col] = inverse_matrix[next_row * row_ind + col] - inverse_matrix[row_offset + col] * temp;
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
void mult(double *a, double *b, double *res, int m1, int m2, int m3, int m, double norm) {
    int count_b = 0;
    UNUSED(m1);
    UNUSED(m2);
    UNUSED(m3);


    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            res[i * m + j] = 0.0;
            double abs_value = fabs(b[i * m + j]);
            if (1e+250 * norm < abs_value || abs_value < 1e-250 * norm) {
                b[i * m + j] = 0.0;
                ++count_b;
            }
        }
    }
    if (count_b == m * m) {
        return;
    }

    for (int k = 0; k < m; ++k) {
        for (int i = 2; i < m; i += 3) {
            double temp = a[i * m + k];
            double temp_1 = a[(i - 1) * m + k];
            double temp_2 = a[(i - 2) * m + k];

            for (int j = 2; j < m; j += 3) {
                double b_j = b[k * m + j];
                double b_jm1 = b[k * m + j - 1];
                double b_jm2 = b[k * m + j - 2];

                res[i * m + j] += b_j * temp;
                res[(i - 1) * m + j] += b_j * temp_1;
                res[(i - 2) * m + j] += b_j * temp_2;

                res[i * m + (j - 1)] += b_jm1 * temp;
                res[(i - 1) * m + (j - 1)] += b_jm1 * temp_1;
                res[(i - 2) * m + (j - 1)] += b_jm1 * temp_2;

                res[i * m + (j - 2)] += b_jm2 * temp;
                res[(i - 1) * m + (j - 2)] += b_jm2 * temp_1;
                res[(i - 2) * m + (j - 2)] += b_jm2 * temp_2;
            }

            for (int j = (m / 3) * 3; j < m; ++j) {
                double b_j = b[k * m + j];
                res[i * m + j] += b_j * temp;
                res[(i - 1) * m + j] += b_j * temp_1;
                res[(i - 2) * m + j] += b_j * temp_2;
            }
        }

        for (int i = (m / 3) * 3; i < m; ++i) {
            double temp = a[i * m + k];
            for (int j = 2; j < m; j += 3) {
                double b_j = b[k * m + j];
                double b_jm1 = b[k * m + j - 1];
                double b_jm2 = b[k * m + j - 2];

                res[i * m + j] += b_j * temp;
                res[i * m + (j - 1)] += b_jm1 * temp;
                res[i * m + (j - 2)] += b_jm2 * temp;
            }

            for (int j = (m / 3) * 3; j < m; ++j) {
                double b_j = b[k * m + j];
                res[i * m + j] += b_j * temp;
            }
        }
    }
}

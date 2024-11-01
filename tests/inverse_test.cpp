#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#define EPS 1e-16
#define CANNOT_SOLVE -1
#define IRREVERSIBLE -2
#define SUCCESS 0

int treug(double * a, double * b, int n, double norma, double* c) {
    int i;
    int j;
    int k;
    int t;
    double p;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i * n + j] = 0;
        }
        for (int i = 0; i < n; i++) {
            b[i * n + i] = 1;
        }
    }

    for (i = 0; i < n; i++) {
        t = -1;
        for (j = i; j < n; j++)
            if (a[j * n + i] >  5e-15 * norma || a[j * n + i] < -5e-15 * norma) {
                //printf("%e\n", a[j * n + i]);
                t = 1;
            }
        if (t == -1) {
            return -1;
        }

        p = a[i * n + i];
        t = i;

        for (j = i; j < n; j++) {
            if (fabs(a[j * n + i]) > fabs(p)) {
                p = a[j * n + i];
                t = j;
            }
        }

        for (j = 0; j < n; j++) c[j] = a[t * n + j];
        for (j = 0; j < n; j++) a[t * n + j] = a[i * n + j];
        for (j = 0; j < n; j++) a[i * n + j] = c[j];

        for (j = 0; j < n; j++) c[j] = b[t * n + j];
        for (j = 0; j < n; j++) b[t * n + j] = b[i * n + j];
        for (j = 0; j < n; j++) b[i * n + j] = c[j];

        p = a[i * n + i];

        for (j = 0; j < n; j++) {
            a[i * n + j] = a[i * n + j] / p;
            b[i * n + j] = b[i * n + j] / p;
        }

        for (j = i + 1; j < n; j++) {
            p = a[j * n + i];
            for (k = 0; k < n; k++) {
                a[j * n + k] = a[j * n + k] - a[i * n + k] * p;
                b[j * n + k] = b[j * n + k] - b[i * n + k] * p;
            }
        }

    }
    return 1;
}

void diag(double * a, double * b, int n) {
    int i, j, k;
    double p;
    for (i = n - 1; i > 0; i--) {
        for (j = 0; j < i; j++) {
            p = a[j * n + i];
            for (k = 0; k < n; k++) {
                a[j * n + k] = a[j * n + k] - a[i * n + k] * p;
                b[j * n + k] = b[j * n + k] - b[i * n + k] * p;
            }
        }
    }
}
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
// Прототипы функций
bool inverse_matrixx(double *matrix, double *inverse_matrix, int n, int m, double norm)
{
     int i;
    int j;
    int k;
    int t;
    double p;
    double *c;
    double *a;
    double *b;

    a = new double[n * n];
    b = new double[n * n];
    c = new double[n];

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            a[i * n + j] = matrix[i * m + j];
            b[i * n + j] = (i == j ? 1 : 0);
        }
    }

    // Проверяем на каждом шаге, что матрица сводится к треугольному виду.
    // Иначе ранк не максимален и нет обратной матрицы.
    for (i = 0; i < n; i++)
    {
        t = -1;

        for (j = i; j < n; j++)
        {
            if (fabs(a[j * n + i]) > 5e-15 * norm)
            {
                t = 1;
            }
        }

        if (t == -1)
        {
            delete[] a;
            delete[] b;
            delete[] c;
            return false;
        }

        p = a[i * n + i];
        t = i;

        // Метод Гаусса с выбором главного элемента по строке
        for (j = i; j < n; j++)
        {
            if (fabs(a[j * n + i]) > fabs(p))
            {
                p = a[j * n + i];
                t = j;
            }
        }

        for (j = 0; j < n; j++)
        {
            c[j] = a[t * n + j];
        }
        for (j = 0; j < n; j++)
        {
            a[t * n + j] = a[i * n + j];
        }
        for (j = 0; j < n; j++)
        {
            a[i * n + j] = c[j];
        }

        for (j = 0; j < n; j++)
        {
            c[j] = b[t * n + j];
        }
        for (j = 0; j < n; j++)
        {
            b[t * n + j] = b[i * n + j];
        }
        for (j = 0; j < n; j++)
        {
            b[i * n + j] = c[j];
        }

        p = a[i * n + i];

        for (j = 0; j < n; j++)
        {
            a[i * n + j] = a[i * n + j] / p;
            b[i * n + j] = b[i * n + j] / p;
        }

        for (j = i + 1; j < n; j++)
        {
            p = a[j * n + i];
            for (k = 0; k < n; k++)
            {
                a[j * n + k] = a[j * n + k] - a[i * n + k] * p;
                b[j * n + k] = b[j * n + k] - b[i * n + k] * p;
            }
        }
    }

    for (i = n - 1; i > 0; i--)
    {
        for (j = 0; j < i; j++)
        {
            p = a[j * n + i];
            for (k = 0; k < n; k++)
            {
                a[j * n + k] = a[j * n + k] - a[i * n + k] * p;
                b[j * n + k] = b[j * n + k] - b[i * n + k] * p;
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            inverse_matrix[i * m + j] = b[i * n + j];
        }
    }

    delete[] a;
    delete[] b;
    delete[] c;
    return true;
}
bool inverse_matrix(double*a, double* b, int n, double norma, double* c)
{
    int t = treug(a, b, n, norma, c);
        if (t == -1) {
            return false;
        }
        diag(a, b, n);
        return true;
}

// Функция для сравнения двух матриц
bool compare_matrices(double* matrix1, double* matrix2, int n, double tolerance) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (fabs(matrix1[i * n + j] - matrix2[i * n + j]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}
void PrintDouble(double* matrix, int n, int r) {
	for(int i = 0; i < std::min(n, r); i++) {
        for(int j = 0; j < std::min(n, r); j++) {
            printf("%10.3e ", matrix[i*n+j]);
        }
        std::cout << std::endl;
    }
	printf("\n");
}
double norma(double* block, int m) {
    int i,j;
    double norm_m=0, temp;
    for(i=0;i<m;i++)  
    {
        temp=0;
        for(j=0;j<m;j++)
            temp+= fabs(block[j*m + i]);
        if(temp>norm_m)
            norm_m=temp;
    }
    return norm_m;
} 

// Тестовая функция
void test_inverse_matrix() {
    const int n = 3; // Размер матрицы
    // Норма матрицы
    const double tolerance = 1e-9; // Допустимая погрешность

    double matrix[n * n];
    double inverse_matrix1[n * n];
    double inverse_matrix2[n * n];
    int index[n];
    double c[n];

    // Заполняем матрицу случайными числами
    std::srand(std::time(nullptr));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i * n + j] = static_cast<double>(std::rand()% 100) ;
        }
    }
    PrintDouble(matrix, n, n);
    double matrix_norm = norma(matrix, n); 
    // Вычисляем обратную матрицу с помощью первой функции
    int result1 = gauss_classic_row(matrix, inverse_matrix1, index, n, matrix_norm, n);

    // Вычисляем обратную матрицу с помощью второй функции
    bool result2 = inverse_matrix(matrix, inverse_matrix2, n, matrix_norm, c);
    PrintDouble(inverse_matrix1, n, n);
    printf("Мое решение\n");
    PrintDouble(inverse_matrix2, n, n);
    std::cout << result1 << std::endl;
    // Проверяем, что обе функции вернули корректные результаты
    if (result1 == SUCCESS && result2) {
        // Сравниваем результаты
        if (compare_matrices(inverse_matrix1, inverse_matrix2, n, tolerance)) {
            std::cout << "Тест пройден: обе функции вернули одинаковые результаты." << std::endl;
        } else {
            std::cout << "Тест не пройден: функции вернули разные результаты." << std::endl;
        }
    } else {
        std::cout << "Тест не пройден: одна из функций не смогла вычислить обратную матрицу." << std::endl;
    }
}

int main() {
    test_inverse_matrix();
    return 0;
}
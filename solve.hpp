#include <iostream>
#include <math.h>
#include <float.h>
#include "multiplication.hpp"
#define UNUSED(x) (void)(x)

void get_block (double *matr, double *block , int n, int m, int i , int j )
{
    int ii=0,jj=0;
    int k = n / m;
    int l = n - k * m;
    int v = (j < k ? m : l), h = (i < k ? m : l);
    //std::cout << "v= " << v << " h= " << h << std::endl;
    //block = new double[v * h];
    double *ma_bl = matr + i * n * m+ j * m;
    for (ii = 0; ii < h; ii++)
    {
        for (jj = 0; jj < v; jj++)
        {
            block[ii *m+jj]=ma_bl[ii *n+jj];
        }   
    }
}
void put_block (double *matr, double *block , int n, int m, int i , int j )
{
int ii=0,jj=0;
int k = n / m;
int l = n - k * m;
int v = (j < k ? m : l), h = (i < k ? m : l);
double *ma_bl = matr + i * n * m+ j * m;
for (ii = 0; ii < h; ii++)
{
    for (jj = 0; jj < v; jj++)
    {
        ma_bl[ii *n+jj]=block[ii *m+jj];
    }   
}
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

bool inverseMatrix(double* matrix, int size, double* inverseMatrix) {
    // Создаем временную матрицу с расширенным размером
    int tempSize = 2*size;
    double* tempMatrix = new double[tempSize*tempSize];

    // Заполняем временную матрицу расширенными данными
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            tempMatrix[i*tempSize + j] = matrix[i*size + j];
            if (i == j) {
                tempMatrix[i*tempSize + j + size] = 1.0;
            } else {
                tempMatrix[i*tempSize + j + size] = 0.0;
            }
        }
    }

    // Приводим временную матрицу к треугольному виду
    for (int i = 0; i < size; i++) {
        if (tempMatrix[i*tempSize + i] == 0.0) {
            // Если главный элемент равен 0, меняем строки местами
            int swapRow = i + 1;
            while (swapRow < size && tempMatrix[swapRow*tempSize + i] == 0.0) {
                swapRow++;
            }
            if (swapRow == size) {
                // Если все элементы в столбце равны 0, матрица необратима
                delete[] tempMatrix;
                return false;
            }
            for (int j = 0; j < tempSize; j++) {
                std::swap(tempMatrix[i*tempSize + j], tempMatrix[swapRow*tempSize + j]);
            }
        }

        double pivot = tempMatrix[i*tempSize + i];
        for (int j = 0; j < tempSize; j++) {
            tempMatrix[i*tempSize + j] /= pivot;
        }

        for (int k = 0; k < size; k++) {
            if (k != i) {
                double factor = tempMatrix[k*tempSize + i];
                for (int j = 0; j < tempSize; j++) {
                    tempMatrix[k*tempSize + j] -= factor * tempMatrix[i*tempSize + j];
                }
            }
        }
    }

    // Копируем обратную матрицу из временной матрицы
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            inverseMatrix[i*size + j] = tempMatrix[i*tempSize + j + size];
        }
    }

    delete[] tempMatrix;
    return true;
}

int findmax(double* matr, double* block, int n, int m,int l, int j) {
    int k = -1;
    bool t;
    double tmp;
    double max_n = DBL_MAX;
    double* inverse = new double[m * m];
    for (int i = l; i < n/m; i++) {
        get_block(matr, block, n ,m , i, j);
        t = inverseMatrix(block, m, inverse);
        if (!t) {
            std::cout << "Обратной матрицы нет  " << std::endl;
            continue;
        }
        tmp = norma(inverse, m);
        std::cout << "Норма " << i << "ой= " << tmp << std::endl;
        if (tmp < max_n) {
            max_n = tmp;
            k = i;
        }
        
    }
    delete [] inverse;
    return k;
}

int solve(int n, int m, double* matrix, double* x) {
    UNUSED(n);
    UNUSED(m);
    UNUSED(matrix);
    UNUSED(x);
    return 0;
}

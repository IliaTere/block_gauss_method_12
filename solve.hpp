#include <iostream>
#include <math.h>
#include <float.h>
#include "mult.hpp"
//#include "reader.hpp"
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

bool inverseMatrix(double* matrix, int n, double* inverseMatrix) {
    // Создаем временную матрицу с расширенным размером
    int tempSize = 2*n;
    double* solution = new double[tempSize*tempSize];

    // Заполняем временную матрицу расширенными данными
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            solution[i*tempSize + j] = matrix[i*n + j];
            if (i == j) {
                solution[i*tempSize + j + n] = 1.0;
            } else {
                solution[i*tempSize + j + n] = 0.0;
            }
        }
    }

    // Приводим временную матрицу к треугольному виду
    for (int i = 0; i < n; i++) {
        if (solution[i*tempSize + i] == 0.0) {
            // Если главный элемент равен 0, меняем строки местами
            int swapRow = i + 1;
            while (swapRow < n && solution[swapRow*tempSize + i] == 0.0) {
                swapRow++;
            }
            if (swapRow == n) {
                // Если все элементы в столбце равны 0, матрица необратима
                delete[] solution;
                return false;
            }
            for (int j = 0; j < tempSize; j++) {
                std::swap(solution[i*tempSize + j], solution[swapRow*tempSize + j]);
            }
        }

        double pivot = solution[i*tempSize + i];
        for (int j = 0; j < tempSize; j++) {
            solution[i*tempSize + j] /= pivot;
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = solution[k*tempSize + i];
                for (int j = 0; j < tempSize; j++) {
                    solution[k*tempSize + j] -= factor * solution[i*tempSize + j];
                }
            }
        }
    }

    // Копируем обратную матрицу из временной матрицы
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverseMatrix[i*n + j] = solution[i*tempSize + j + n];
        }
    }

    delete[] solution;
    return true;
}
// l - 
int findmax(double* matr, double* block, int n, int m,int l, int j) { // TODO:  подавать на вход inverse
    int k = -1;
    bool t;
    double tmp;
    double max_n = DBL_MAX;
    double* inverse = new double[m * m];
    for (int i = l; i < n/m; i++) {
        get_block(matr, block, n ,m , i, j);
        t = inverseMatrix(block, m, inverse);
        if (!t) {
            //std::cout << "Обратной матрицы нет  " << std::endl;
            continue;
        }
        tmp = norma(inverse, m);
        //std::cout << "Норма " << i << "ой= " << tmp << std::endl;
        if (tmp < max_n) {
            max_n = tmp;
            k = i;
        }
        
    }
    delete [] inverse;
    return k;
}

void swap_rows(double* matr, int k, int l, int n, int m) {
    double* tmp = new double[m*n];
    for (int i = 0; i < m*n; i++) {
        tmp[i] = matr[k*m * n + i];
    }
    for (int i = 0; i < m*n; i++) {
        matr[m*k * n + i] = matr[m*l * n + i];
    }
    for (int i = 0; i < m*n; i++) {
        matr[m*l * n + i] = tmp[i];
    }
    delete[] tmp;
}

void subtraction(double* a, double* b, int m) {
    for (int i = 0; i < m*m; i++) {
        a[i] -= b[i];
    }

    
}

int solve(int n, int m, double* matr, double* block, double* solution, double* inverse, double* tmp, double* block1, double* block2) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                solution[i*n + j] = 1.0;
            } else {
                solution[i*n + j] = 0.0;
            }
        }
    }
    inverseMatrix(matr, n, inverse);
    PrintDouble(inverse, n, n);
    int t = -1; // ToDo: можем итерироваться до size который менятеся в зависимости от n, m
    for (int i = 0; i < n/m ; i++) {
        t = findmax(matr, block, n, m, i, i);
        //std::cout << "t= "<< t  << "i= "<< i<< std::endl;
        if (t == -1) {
            std::cout << "Необратимая матрица"  << std::endl;
        }
        if (t != i) {
        swap_rows(matr, i, t, n, m);
        swap_rows(solution, i, t, n, m);
        }
        get_block(matr, block, n, m, i, i);
        
        inverseMatrix(block, m, inverse);
        // printf("--------Обратная матрица\n");
        // PrintDouble(block, m ,m);
        // PrintDouble(inverse, m ,m);
        // printf("\n");
        // printf("Норма= %f\n", norma(inverse, m));
        
        for (int j=i; j < n/m; j++) { // Заполняем solution с i+1 эл-та
            get_block(matr, block, n, m, i, j);
            mult(inverse, block, tmp, m, m, m ,m);
            put_block(matr, tmp, n,m, i, j);
            get_block(solution, block, n, m, i, j);
            mult(inverse, block, tmp, m, m, m ,m);
            put_block(solution, tmp, n,m, i, j);
        }                                       // Правильно +
        
        for (int k=i+1; k < n/m - i; k++) {
            for (int f=i+1; f < n/m - i; f++) {
                get_block(matr,block, n, m, i+1, i); // TODO: посмотреть можно ли брать этот блок раньше
                get_block(matr,block1, n, m, k, f);
                get_block(matr,block2, n,m, k-1, f);
                mult(block, block2, tmp, m , m, m, m);
                // printf("Умножаем \n");
                // PrintDouble(block, m ,m);
                // PrintDouble(block2, m ,m);
                // printf("Получаем и вычитаем из  \n");
                // PrintDouble(tmp, m ,m);
                // PrintDouble(block, m ,m);
                subtraction(block1, tmp, m);
                // printf("Получаем и записываем на (%d,%d)\n", k ,f);
                // PrintDouble(block1, m, m);
                put_block(matr, block1, n, m, k, f);

                get_block(solution,block, n, m, i+1, i+1); 
                get_block(solution,block1, n, m, k, f);
                get_block(solution,block2, n,m, k-1, f);
                mult(block, block2, tmp, m , m, m, m);
                subtraction(block1, tmp, m);
                put_block(solution, block1, n, m, k, f);
            }
        }
        
    }
    // Обратный ход
        for (int i = n/m-1; i > 0; i--) {
            for (int j = i-1; j > -1; j--) {
                get_block(matr,block,n,m,j,i);
                get_block(solution,block1,n,m,j,i);
                mult(block1, block, tmp, m,m,m,m);
                subtraction(block1, tmp, m);
                put_block(solution, block1,n , m, j , i);
                printf("(%d, %d)\n", i, j);
                PrintDouble(block1, m, m);   
            }
        }

    printf("-----------------------ОТВЕТ\n");
    
    PrintDouble(matr, n ,n);
    PrintDouble(solution, n , n);
    return 0;
}

#pragma once
#include <iostream>
#include <math.h>
#include <float.h>
#include "mult.hpp"
#include <iomanip>
#include <algorithm>
#include <immintrin.h>
#include <cstring>
#define UNUSED(x) (void)(x)

using namespace std;
void get_block(double *matrix, double *block, int n, int m, int i_block, int j_block)
{
    int k = n / m;
    int l = n % m;
    int h = (i_block < k ? m : l);
    int v = (j_block < k ? m : l);
    
    std::memset(block, 0, m * m * sizeof(double));
    for (int i = 0; i < h; ++i) {
        std::copy_n(&matrix[n * (i_block * m + i) + j_block * m], v, &block[i * m]);
    }
    // if(block[0] > 1e+300 || block[0] < 1e-300) {
    //     printf("(%2d,%2d| %lf)", i_block, j_block, block[0]);
    // }
}
void put_block(double *matrix, double *block, int n, int m, int i_block, int j_block)
{
	int k = n / m;
	int l = n % m;
	int h = (i_block < k ? m : l);
	int v = (j_block < k ? m : l);
	int tmp;
	for (int i = 0; i < h; i++)
	{
		tmp = i * m;
		for (int j = 0; j < v; j++)
			matrix[n * (i_block * m + i) + j_block * m + j] = block[tmp + j];
	}
}

double norma(double *matrix, int n)
{
	double norm = 0;
    double tmp_sum = 0;
	for (int i = 0; i < n; i++)
	{
		tmp_sum = 0;
		for (int j = 0; j < n; j++)
			tmp_sum = tmp_sum + fabs(matrix[j * n + i]);
		if (tmp_sum - norm > 0)
			norm = tmp_sum;
	}
	return norm;
}

double ravno(double x, double norma)
{
    if(fabs(x)<1e-7*norma){
        return 0.0;}
    else {return x;}
}

void ravnoBlock(double *block, int m, double norma)
{
    int i,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            block[j*m + i] = ravno(block[j*m + i], norma);
        }
    }
}
void pcord(int i, int j)
{
    printf("(%d, %d)\n", i, j);
}

void swap_rows(double* matr, int k, int l, int n, int m) {
    int start_k = k * m * n;
    int start_l = l * m * n;

    // Обмен блочными строками
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            std::swap(matr[start_k + i * n + j], matr[start_l + i * n + j]);
        }
    }
}

void subtraction(double *a, double *b, int i_block, int j_block, int k, int m, int l)
{
	int h = (i_block < k ? m : l);
	int v = (j_block < k ? m : l);
	int tmp;
	for (int i = 0; i < h; i++)
	{
		tmp = i * m;
		for (int j = 0; j < v; j++)
			a[tmp + j] = a[tmp + j] - b[tmp + j];
	}
}
void initialize_matrix(double* solution, int n) {
    std::memset(solution, 0, n * n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        solution[i * n + i] = 1.0;
    }
}

void copyMatrix(double* in, double* out, int current_size, int req_size) {
    for (int i = 0; i < current_size*current_size; i++)
    {
        out[i] = 0;
    }
    for (int i = 0; i < req_size; ++i) {
        for (int j = 0; j < req_size; ++j) {
            out[i * req_size + j] = in[i * current_size + j];
        }
    }
}
void setZeroBlock(double *matrix, double *block, int i_block, int j_block, int n, int m)
{
    for (int row = 0; row < m; row++)
    {
        for (int col = 0; col < m; col++)
        {
            if (row == col && i_block == j_block)
            {
                block[row * m + col] = 1.0;
            }
            else
            {
                block[row * m + col] = 0.0;
            }
        }
    }

    put_block(matrix, block, n, m, i_block, j_block);
}


int solve(int n, int m, double* matr, double* block, double* solution, double* inverse, double* tmp, double* block1, /*double* block2,*/ double matrix_norm, int *block_index) {
    int s, ss, i, j;
    int k = n/m;
    int l = n%m;
    int bl = (l==0?k:k+1);
    double buffer_norma, min_norma = 0;
    int min_index = 0;

    for (int p=0; p<k; p++) {
        min_norma = 0;
        min_index = -1;
        for(j=p; j < k; j++)
        {
            buffer_norma = 0;
            get_block(matr, block, n, m, j, p);
            if (inverse_matrix(block, inverse, block_index, m, matrix_norm, m) != -1)
            {
                buffer_norma = norma(inverse, m);
                if (min_norma > buffer_norma || min_index == -1)
                {
                    min_norma = buffer_norma;
                    min_index = j;
                }
            }
        }
        if (min_index == -1) {
            printf("Метод не применим\n");
            return -1;
        }
        if (min_index != p) {
            swap_rows(matr, min_index, p, n, m);
            swap_rows(solution, min_index, p, n, m);
            // printf("После перестановки\n");
        }

        get_block(matr, block, n, m, p, p);
        // PrintDouble(block, m , m);
        inverse_matrix(block, inverse, block_index, m, matrix_norm, m);
        get_block(matr, block, n, m, p, p);
        setZeroBlock(matr, block, p, p, n , m);
        for (s = p+1; s < bl; s++) 
        {
            get_block(matr, block, n, m, p, s);

            // PrintDouble(block, m, m);
            // PrintDouble(inverse, m,m);

            mult(inverse, block, tmp, m, m, s==k?l:m, m, matrix_norm);

            // PrintDouble(tmp, m, m);

            put_block(matr, tmp, n, m, p, s);
        }
        for (s = 0; s < bl; s++) 
        {
            get_block(solution, block, n, m, p, s);
            mult(inverse, block, tmp, m, m, m, m, matrix_norm);
            put_block(solution, tmp, n, m, p, s);
        }
        // PrintDouble(matr, n , n);
        for (s = p+1; s < bl; s++)
        {
            get_block(matr, block, n , m, s, p);
            get_block(matr, block1, n, m, s, p);
            setZeroBlock(matr, block1, s, p, n, m);
            for (ss=p+1; ss<bl; ss++) // Тут с p+1
            {
                get_block(matr, block1, n , m, p, ss);
                
                // printf("Начало умножения\n");
                // PrintDouble(block1, m , m);
                // PrintDouble(block, m , m);
                
                mult(block, block1, tmp, m, m, m, m, matrix_norm);
                
                // PrintDouble(tmp, m , m);
                // printf("Конец умножения, начало вычитания\n");

                get_block(matr, block1, n ,m, s, ss);
                
                // pcord(s, ss);
                // PrintDouble(block1, m ,m);
                // printf(" - \n");
                // PrintDouble(tmp, m, m);
                
                subtraction(block1, tmp, s, ss, k, m, l);
                
                // PrintDouble(block1, m, m);
                
                put_block(matr, block1, n, m, s, ss);
            }
            for (ss=0; ss<bl; ss++) // Тут с p+1
            {
                get_block(solution, block1, n , m, p, ss);
                mult(block, block1, tmp, m , m, m, m, matrix_norm);
                get_block(solution, block1, n ,m, s, ss);
                subtraction(block1, tmp, s, ss, k, m, l);
                put_block(solution, block1, n, m, s, ss);
                
            }
        }
    }
    // PrintDouble(matr, n, n);
    // PrintDouble(solution, n, n);
    // printf("До if ###################\n");
    if (l!=0)
    {
        get_block(matr, block, n ,m, k, k);
        // PrintDouble(block, m, m);
        if (inverse_matrix(block, inverse, block_index, l, matrix_norm, m) == -1)
        {
            return -1;
        }
        setZeroBlock(matr, block1, k, k, n, m);
        for (j = 0; j <= k; j++) {
            get_block(solution, block, n, m, k, j);
            // PrintDouble(block, m, m);
            // PrintDouble(inverse, m, m);
            mult(inverse, block, tmp, l, l, j==k?l:m, m, matrix_norm);
            
            put_block(solution, tmp, n ,m, k, j);
        }
    }
    for ( i = bl - 1; i > 0; i--)
    {
        for ( j = i - 1; j >= 0; j--)
        {
            get_block(matr, inverse, n, m, j, i);
            for ( ss = 0; ss < bl; ss++)
            {
                get_block(solution, block1, n, m, i, ss);
                mult(inverse, block1, block, m, m, m, m, matrix_norm);
                get_block(solution, block1, n, m, j, ss);
                subtraction(block1, block, j, ss, k, m, l);
                put_block(solution, block1, n, m, j, ss);
            }
        }
    }
    return 0;
}

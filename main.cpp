#include <iostream>
#include <string>
#include <fstream>
#include "reader.hpp"
#include "solve.hpp"
#include "formula.hpp"
#include <cstring>
bool isNumber(std::string& str)
    {
    std::string::iterator it = std::begin(str);
    while (it != str.end() && std::isdigit(*it)) {
        it++;
    }
    return !str.empty() && it == str.end();
    }

int main(int argc, char **argv)
{
    
    if (argc < 5 || argc > 6) {
        std::cout << "error: To many(few) arguments";
        return -1;
    }
    for (int i = 1; i < 5; i++) {
        std::string str(argv[i]);
        if (isNumber(str) == false) {
            std::cout << "error: Invalid argument: " << str << std::endl;
            return -4;
        }
    }
    std::string s1(argv[1]);
    std::string s2(argv[2]);
    std::string s3(argv[3]);
    int n = std::stoi(s1); // Размерность матрицы
    int m = std::stoi(s2); // Размерность блока
    int r = std::stoi(s3); // Кол-во выводимых значений в матрице
	double* matr = new double[n*n];
    double* x = new double[n*n];
    if (strcmp(argv[4],"0") == 0) {
        if (argc != 6) {
            std::cout << "error: File not found" << std::endl;
            delete[] matr;
            delete[] x;
            return -2;
        }
        std::string name(argv[5]);
		int t = read_ff(name , matr , n*n);
        if(t != 0) {
            delete[] matr;
            delete[] x;
            return -5;
        }
        std::cout << std::endl;
        PrintDouble(matr,n,r);
        double* block = new double[m*m];
        double* inverse = new double[m*m];
        //get_block(matr, block, n, m, 1 , 0);
        int k = findmax(matr, block, n, m, 0, 0);
        std::cout << "--------------------------------" << " Номер строки в которой норма минимальная " << k << std::endl;
        //PrintDouble(block, m, m);
        delete[] block;
        delete[] inverse;
        return solve(n,m,matr,x);
    }
    if ((strcmp(argv[4],"0") != 0)) 
    {
        if (argc > 5) {
            std::cout << "error: To many(few) arguments" << std::endl;
            delete[] matr;
            delete[] x;
            return -6;
        }
        std::string tmp(argv[4]);
        int s = stoi(tmp);
        if (s < 1 || s > 4) {
            std::cerr << "error: Pasametr s is not a valid number" << std::endl;
            delete[] matr;
            delete[] x;
            return -7;
        }
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                matr[n*i + j] = formula(s,n,i+1,j+1);
            }
        }
        PrintDouble(matr, n, r);
        double* block = new double[m*m];
        double* inverse = new double[m*m];
        get_block(matr, block, n, m, 1 , 0);
        int k = findmax(matr, block, n, m, 0, 0);
        std::cout << "--------------------------------" << " Номер строки в которой норма минимальная " << k << std::endl;
        //PrintDouble(block, m, m);
        delete[] block;
        delete[] inverse;
        return solve(n,m,matr,x);
    }
    delete[] matr;
    delete[] x;
    return 0;
}

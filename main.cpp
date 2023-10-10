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
        std::cout << "error: To much arguments";
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
    std::vector<double> matrix;
    std::vector<double> x;
    if (strcmp(argv[4],"0") == 0) {
        if (argc != 6) {
            std::cout << "error: File not found" << std::endl;
            return -2;
        }
        std::string name(argv[5]);
        int t = read_doubles_from_file(name, matrix, n*n);
        if(t != 0) {
            return -5;
        }
        printMatrix(matrix, n, r);
        return solve(n,m,matrix,x);
    }
    if ((strcmp(argv[4],"0") != 0)) 
    {
        if (argc > 5) {
            std::cout << "error: To much arguments" << std::endl;
            return -6;
        }
        std::string tmp(argv[4]);
        int s = stoi(tmp);
        if (s < 1 || s > 4) {
            std::cerr << "error: Pasametr s is not a valid number" << std::endl;
            return -7;
        }
        for(int i = 1; i < n+1; i++) {
            for(int j = 1; j < n+1; j++) {
                matrix.push_back(formula(s, n, i, j));
            }
        }
        printMatrix(matrix,n,r);
        return solve(n,m,matrix,x);
    }

    return 0;
}

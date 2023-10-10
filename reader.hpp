#include <iostream>
#include <cstdio>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

bool is_double(const std::string& s) {
    if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+') && (s[0] != '.'))) return false;
    char* p;
    strtod(s.c_str(), &p);
    return (*p == 0);
}

 int read_doubles_from_file(const std::string& filename, std::vector<double>& result, int n) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return -1;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            if (is_double(token)) {
                result.push_back(std::stod(token));
            } else {
                std::cerr << "Invalid double value: " << token << std::endl;
                return -2;
            }
        }
    }
    if (result.size() != (size_t)n) {
        std::cout << "Invalid lenght" << std::endl;
        return -3;
    }
    return 0;
}

void printMatrix(const std::vector<double> matrix, int n, int r) {
    for(int i = 0; i < std::min(n, r); i++) {
        for(int j = 0; j < std::min(n, r); j++) {
            printf("%10.3e ", matrix[i*n+j]);
        }
        std::cout << std::endl;
    }
}

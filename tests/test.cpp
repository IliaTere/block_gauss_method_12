#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <Eigen/Dense>
#include"v1_mult.hpp"
#include"mult.hpp"

// Функция для проверки корректности работы mult и измерения времени
bool test_mult(int n, int m) {
    // Создаем случайные матрицы A и B
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = dis(gen);
            B(i, j) = dis(gen);
        }
    }

    // Создаем матрицы для результатов
    Eigen::MatrixXd C_eigen = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd C_custom = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd C_custom2 = Eigen::MatrixXd::Zero(n, n);

    // Преобразуем Eigen-матрицы в массивы для передачи в mult
    double *a = new double[n * n];
    double *b = new double[n * n];
    double *c = new double[n * n];
    double *c2 = new double[n * n];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i * n + j] = A(i, j);
            b[i * n + j] = B(i, j);
        }
    }

    // Измеряем время работы умножения через Eigen
    auto start_eigen = std::chrono::high_resolution_clock::now();
    C_eigen = A * B;
    auto end_eigen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_eigen = end_eigen - start_eigen;
    std::cout << "Eigen multiplication time: " << duration_eigen.count() << " seconds" << std::endl;

    // Измеряем время работы вашей функции mult
    auto start_custom = std::chrono::high_resolution_clock::now();
    mult1(a, b, c, n, m);
    auto end_custom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_custom = end_custom - start_custom;
    std::cout << "Custom multiplication time: " << duration_custom.count() << " seconds" << std::endl;

    // Измеряем время работы другой функции mult
    auto start_custom2 = std::chrono::high_resolution_clock::now();
    mult(a, b, c2, n, n, n, n, 1.0); // Предположим, что норма равна 1.0
    auto end_custom2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_custom2 = end_custom2 - start_custom2;
    std::cout << "Custom2 multiplication time: " << duration_custom2.count() << " seconds" << std::endl;

    // Преобразуем результат обратно в Eigen-матрицу
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C_custom(i, j) = c[i * n + j];
            C_custom2(i, j) = c2[i * n + j];
        }
    }

    // Выводим первые 5 элементов каждого решения
    // std::cout << "First 5 elements of Eigen result:" << std::endl;
    // for (int i = 0; i < 5; ++i) {
    //     std::cout << C_eigen(i, 0) << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "First 5 elements of Custom result:" << std::endl;
    // for (int i = 0; i < 5; ++i) {
    //     std::cout << C_custom(i, 0) << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "First 5 elements of Custom2 result:" << std::endl;
    // for (int i = 0; i < 5; ++i) {
    //     std::cout << C_custom2(i, 0) << " ";
    // }
    // std::cout << std::endl;

    // Сравниваем результаты
    bool success = true;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (std::abs(C_eigen(i, j) - C_custom(i, j)) > 1e-6) {
                success = false;
                std::cout << "Error at (" << i << ", " << j << ") with Custom: "
                          << C_eigen(i, j) << " != " << C_custom(i, j) << std::endl;
            }
            if (std::abs(C_eigen(i, j) - C_custom2(i, j)) > 1e-6) {
                success = false;
                std::cout << "Error at (" << i << ", " << j << ") with Custom2: "
                          << C_eigen(i, j) << " != " << C_custom2(i, j) << std::endl;
            }
        }
    }

    // Освобождаем память
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] c2;

    return success;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <n> <m>" << std::endl;
        return 1;
    }

    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);

    if (test_mult(n, m)) {
        std::cout << "Test passed!" << std::endl;
    } else {
        std::cout << "Test failed!" << std::endl;
    }

    return 0;
}
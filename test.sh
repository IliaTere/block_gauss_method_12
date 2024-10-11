#!/bin/bash

# Очищаем файл результатов перед началом
> results.txt

for input_formula in {1..4} # изменяем размер матрицы
do
    for block_size in {1..10} # изменяем размер блока
    do
        for matrix_size in {10..20} # изменяем размер вывода
        do
            echo "Running with matrix_size=$matrix_size, block_size=$block_size, input_formula=$input_formula" >> results.txt
            ./a.out $matrix_size $block_size 0 $input_formula >> results.txt 2>&1
            echo "----------------------------------------" >> results.txt
        done
    done
done
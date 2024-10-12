#!/bin/bash

# Проверяем, что переданы два аргумента
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file1> <file2>"
    exit 1
fi

file1="$1"
file2="$2"
output_file="differences.txt"

# Очищаем файл с отличиями перед началом
> "$output_file"

# Сравниваем строки в двух файлах с помощью 'comm'
comm -3 <(sort "$file1") <(sort "$file2") > "$output_file"

echo "Differences have been written to $output_file"
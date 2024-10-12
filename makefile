# Переменные
CC = g++
CFLAGS = -O3 -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=2 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-property-attribute-mismatch

# Исходные файлы
SRCS = src/main.cpp
OBJS = $(SRCS:.cpp=.o)

# Имя исполняемого файла
TARGET = my_program

# Правило по умолчанию
all: $(TARGET)

# Правило для сборки исполняемого файла
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

# Правило для компиляции отдельных .cpp файлов в .o файлы
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Правило для очистки
clean:
	rm -f $(OBJS) $(TARGET)

# Правило для перекомпиляции всего проекта
rebuild: clean all

# Правило для запуска программы
run: $(TARGET)
	./$(TARGET)
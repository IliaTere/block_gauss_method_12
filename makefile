# Переменные
CC = g++
CFLAGS = -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
CFLAGS_NO_UNUSED = -O0 -fstack-protector-all -g -W -Wall -Wextra -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=2 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-property-attribute-mismatch

# Исходные файлы
SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)

# Имя исполняемого файла
TARGET = a.out

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

# Правило для компиляции без флагов -Wunused
no_unused:
	$(CC) $(CFLAGS_NO_UNUSED) -c main.cpp -o main.o
	$(CC) $(CFLAGS_NO_UNUSED) -o $(TARGET) main.o
# Правило для запуска тестов
test: $(TARGET)
	rm tests/$(TARGET)
	cp $(TARGET) tests/$(TARGET)
	cd tests && ./test_small.sh

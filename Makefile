CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -march=native -Wall -Wextra -pedantic
TARGET := nprime
SRC := src/nprime.cpp

.PHONY: all clean test bench bench-scale

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

test: $(TARGET)
	./$(TARGET) --self-test

bench: $(TARGET)
	./$(TARGET) --bench 1000000000 10000
	./$(TARGET) --bench 999999487 10000

bench-scale: $(TARGET)
	./$(TARGET) --bench 1000000000 10000
	./$(TARGET) --bench 9999999999 10000
	./$(TARGET) --bench 10000000000 10000
	./$(TARGET) --bench 20000000000 1

clean:
	rm -f $(TARGET)

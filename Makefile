CXX = g++
CXXFLAGS = -O3 -std=c++17 -IInclude
OMPFLAGS = -O3 -march=native -flto -std=c++17 -fopenmp -IInclude  -static-libstdc++

SRC = Src/Main.cpp $(shell find Lib -name '*.cpp')
TARGET = Bin/bucket-partitioned-MDS

all: $(TARGET)

$(TARGET): $(SRC)
	@mkdir -p Bin
	$(CXX) $(OMPFLAGS) $(SRC) -o $(TARGET)
	@echo "✅ Build successful: $(TARGET)"

clean:
	rm -f $(TARGET)
	@echo "🧹 Cleaned build artifacts."

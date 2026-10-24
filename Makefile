CXX      := g++
CXXFLAGS := -O3 -std=c++17 -I include/ -I $(mkEigenInc)

BUILD_DIR := build
TARGET    := $(BUILD_DIR)/main
SRCS      := src/main.cpp

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRCS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(BUILD_DIR)
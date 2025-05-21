CXX = /opt/homebrew/bin/g++-14

CXXFLAGS = -std=c++17 -O3
TARGET = build/test-numblib.x
SRCS = main.cpp

OBJS = $(patsubst %.cpp, build/%.o, $(SRCS))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

build/%.o: %.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean
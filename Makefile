# Compiler
CXX = g++

# Compiler flags
# -std=c++20 (from your comment)
# -Wall (Enable all warnings)
# -g (Include debug symbols)
CXXFLAGS = -std=c++20 -Wall -g -O2

# Linker flags
LDFLAGS =

# Source files
SOURCES = main.cpp circuit.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE = main

# Default target: build the executable
all: $(EXECUTABLE)

# Rule to link the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

# Rule to compile .cpp files into .o files
%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Specific rule for main.o (it only depends on circuit.h)
main.o: main.cpp circuit.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Specific rule for circuit.o (it depends on circuit.h)
circuit.o: circuit.cpp circuit.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean target: remove compiled files
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

# Phony targets
.PHONY: all clean
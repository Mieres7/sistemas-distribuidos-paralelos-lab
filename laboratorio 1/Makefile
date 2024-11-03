# Nombre del compilador
CXX = g++

# Flags de compilación
CXXFLAGS = -Wall -O2 -fopenmp

# Librerías necesarias
LIBS = -lcfitsio -lm

# Objetivo
TARGET = main

# Reglas
all: $(TARGET)

$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.cpp $(LIBS)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET) -i elipsesimple1.fits -d 4 -u 4 -b 128 -a 50 -r 0.4

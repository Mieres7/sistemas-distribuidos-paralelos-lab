# Nombre del compilador
CXX = g++

# Flags de compilación
CXXFLAGS = -Wall -g

# Librerías necesarias
LIBS = -lcfitsio

# Objetivo
TARGET = main

# Reglas
all: $(TARGET)

$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.cpp $(LIBS)

clean:
	rm -f $(TARGET)

IN=src/complex.cpp src/polynomial.cpp src/newton_poly.cpp src/main.cpp
CC=g++-13 -std=c++20 -lm -O3
OUT=main

default:
	$(CC) -o $(OUT) $(IN) && ./$(OUT)

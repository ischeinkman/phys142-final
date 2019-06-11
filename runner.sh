#!/bin/sh
clang -Wall -o out/main main.c hermite_polynomial.c -lm 
./out/main

g++ -march=core2 -ffast-math -pthread -std=c++11 -DSFMT_MEXP=607 -I ./SFMT-src-1.4.1/ -I ./../eigen-3.3.7/ -O3 -o MCMR SFMT.c main.cpp 

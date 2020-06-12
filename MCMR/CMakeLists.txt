cmake_minimum_required(VERSION 3.6)
project(TopPPR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    SFMT-src-1.4.1/SFMT-alti.h
    SFMT-src-1.4.1/SFMT-common.h
    SFMT-src-1.4.1/SFMT-params.h
    SFMT-src-1.4.1/SFMT-params11213.h
    SFMT-src-1.4.1/SFMT-params1279.h
    SFMT-src-1.4.1/SFMT-params132049.h
    SFMT-src-1.4.1/SFMT-params19937.h
    SFMT-src-1.4.1/SFMT-params216091.h
    SFMT-src-1.4.1/SFMT-params2281.h
    SFMT-src-1.4.1/SFMT-params4253.h
    SFMT-src-1.4.1/SFMT-params44497.h
    SFMT-src-1.4.1/SFMT-params607.h
    SFMT-src-1.4.1/SFMT-params86243.h
    SFMT-src-1.4.1/SFMT-sse2-msc.h
    SFMT-src-1.4.1/SFMT-sse2.h
    SFMT-src-1.4.1/SFMT.c
    SFMT-src-1.4.1/SFMT.h
    SFMT-src-1.4.1/test.c
    alias.h
    Graph.h
    ppr.h
    main.cpp
    Random.h weightedGraph.h manifoldRanking.h graphAlias.h algo.h SFMT-src-1.4.1/heap.h)

add_executable(TopPPR ${SOURCE_FILES})
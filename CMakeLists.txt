cmake_minimum_required(VERSION 3.26)
project(BigNumsLibrary)

set(CMAKE_CXX_STANDARD 20)

add_executable(BigNumsLibrary main.cpp
        BigNum.cpp
        BigNum.h
        BigNumError.cpp
        BigNumError.h)

cmake_minimum_required(VERSION 3.15)
project(simulador)

find_package(pybind11 REQUIRED)

pybind11_add_module(simulador
    simulador.cpp
    binding.cpp
)
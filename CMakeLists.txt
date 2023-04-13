cmake_minimum_required(VERSION 3.24)
project(NumericMethods)

set(CMAKE_CXX_STANDARD 17)

set(UTILS
        week01_tasks/utils.hpp
        week01_tasks/utils.cpp week03_tasks/optimization.cpp)
add_executable(kahan week01_tasks/kahan.cpp ${UTILS} week03_tasks/optimization.cpp)
add_executable(poly week01_tasks/poly.cpp ${UTILS} week03_tasks/optimization.cpp)
add_executable(pairwise_sum week01_tasks/pairwise_sum.cpp ${UTILS} week03_tasks/optimization.cpp)
add_executable(statistics week01_tasks/statistics.cpp week01_tasks/statistics.hpp week03_tasks/optimization.cpp)
add_executable(length week01_tasks/vector_length.cpp week03_tasks/optimization.cpp)
add_executable(opt week03_tasks/optimization.cpp)
add_executable(tank week03_tasks/tank.cpp)
add_executable(global_opt week04_tasks/global_optimization.cpp)
add_executable(equation_solving week04_tasks/equation_solving.cpp)
add_executable(comp_geom week15_tasks/computational_geometry.cpp)
add_executable(integration week07_tasks/integration.cpp week07_tasks/integration.cpp)
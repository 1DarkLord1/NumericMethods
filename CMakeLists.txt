cmake_minimum_required(VERSION 3.24)
project(NumericMethods)

set(CMAKE_CXX_STANDARD 17)

set(UTILS
        week01_tasks/utils.hpp
        week01_tasks/utils.cpp)
add_executable(kahan week01_tasks/kahan.cpp ${UTILS})
add_executable(poly week01_tasks/poly.cpp ${UTILS})
add_executable(pairwise_sum week01_tasks/pairwise_sum.cpp ${UTILS})
add_executable(statistics week01_tasks/statistics.cpp week01_tasks/statistics.hpp)
add_executable(length week01_tasks/vector_length.cpp)
cmake_minimum_required(VERSION 2.6)

project(HEmpute LANGUAGES CXX)

add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-std=c++17)
add_compile_options(-march=native)
add_compile_options(-g)


add_executable("hefoo" "src/main_hempute.cpp" "src/param.cpp" "src/utils.cpp" "src/thread.cpp" "src/utils_data.cpp" "src/HEmpute_bfv.cpp" "src/HEmpute_ckks.cpp" "src/TestHEmpute.cpp")


# Import Microsoft SEAL
find_package(SEAL 3.4.5 EXACT REQUIRED
    # Providing a path so this can be built without installing Microsoft SEAL
    PATHS external/seal/native/src/cmake
)

TARGET_LINK_LIBRARIES("hefoo" SEAL::seal)

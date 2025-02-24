cmake_minimum_required(VERSION 3.25)

set(PRJ_NAME onednn_kernel)
project(${PRJ_NAME})

MESSAGE(STATUS "== Start: ${PRJ_NAME}")

find_package(dnnl REQUIRED PATHS ${CMAKE_CURRENT_SOURCE_DIR}/oneDNN/build/install/lib/cmake/dnnl/)

file(GLOB SRC 
	"./*.cpp")

add_executable(${PRJ_NAME}
	${SRC}
)

target_link_libraries(${PRJ_NAME} PRIVATE DNNL::dnnl)
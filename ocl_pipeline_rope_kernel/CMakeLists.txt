cmake_minimum_required(VERSION 3.20)
set(PRJ_NAME ocl_pipeline_rope_kernel)

project(${PRJ_NAME} LANGUAGES C CXX)

find_package(OpenCL REQUIRED)

file(GLOB SRC 
	"src/*.cpp")

add_executable(${PRJ_NAME} ${SRC})
target_link_libraries(${PRJ_NAME} PRIVATE OpenCL::OpenCL)
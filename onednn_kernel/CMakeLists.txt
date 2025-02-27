cmake_minimum_required(VERSION 3.25)

if(GPU_RUNTIME STREQUAL "SYCL")
	set(PRJ_NAME onednn_kernel_sycl)
elseif (GPU_RUNTIME STREQUAL "OCL")
	set(PRJ_NAME onednn_kernel_ocl)
else()
	MESSAGE(FATAL_ERROR "Please specific GPU_RUNTIME")
endif()

project(${PRJ_NAME})

# Define the macro
macro(DEFINE_CPP_MACRO MACRO_NAME MACRO_VALUE)
  target_compile_definitions(${PRJ_NAME} PRIVATE ${MACRO_NAME}=${MACRO_VALUE})
endmacro()

MESSAGE(STATUS "== Start: ${PRJ_NAME}")
MESSAGE(STATUS "GPU_RUNTIME = ${GPU_RUNTIME}")
if(GPU_RUNTIME STREQUAL "SYCL")
    find_package(dnnl REQUIRED PATHS ${CMAKE_CURRENT_SOURCE_DIR}/onednn_sycl/oneDNN/build/install/lib/cmake/dnnl/)
	find_package(IntelSYCL REQUIRED)
elseif (GPU_RUNTIME STREQUAL "OCL")
	find_package(dnnl REQUIRED PATHS ${CMAKE_CURRENT_SOURCE_DIR}/onednn_ocl/oneDNN/build/install/lib/cmake/dnnl/)
endif()

file(GLOB SRC 
	"./*.cpp")

add_executable(${PRJ_NAME}
	${SRC}
)

add_sycl_to_target(TARGET ${PRJ_NAME}  SOURCES ${SRC} )

target_link_libraries(${PRJ_NAME} PRIVATE DNNL::dnnl)

if(GPU_RUNTIME STREQUAL "SYCL")
	DEFINE_CPP_MACRO(GPU_RUNTIME_SYCL ${GPU_RUNTIME})
endif()
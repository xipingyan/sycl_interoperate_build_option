# HelloOpenCL

Verify OpenCL rope kernel, compare result with SYCL runtime.

# How to run

    source /opt/intel/oneapi/setvars.sh 
    cd sycl_interoperate_build_option/ocl_pipeline_rope_kernel
    mkdir build && cd build
    cmake ..
    make -j20
# Sycl pipeline

#### build

    cd sycl_interoperate_build_option/sycl_pipeline_rope_kernel
    mkdir build && cd build
    source /opt/intel/oneapi/setvars.sh 
    cmake -DCMAKE_CXX_COMPILER=icpx ..

    SYCL_UR_TRACE=3
    ./sycl_pipeline_rope_kernel
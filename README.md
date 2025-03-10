# sycl_interoperate_build_option
Just compare the performance between SYCL pipeline and OCL pipeline if they have some OpenCL C kernel source.

# Run OpenCL pipeline

    source /opt/intel/oneapi/setvars.sh 
    cd sycl_interoperate_build_option/ocl_pipeline_sdpa_kernel/
    mkdir build && cd build
    cmake ..
    make -j20
    ../run.sh 

    == Infer 147, time = 1208 micro sec.
    == Infer 148, time = 1211 micro sec.
    == Infer 149, time = 1219 micro sec.

# Run SYCL pipeline

    source /opt/intel/oneapi/setvars.sh 
    cd sycl_interoperate_build_option/sycl_pipeline_sdpa_kernel/
    mkdir build && cd build
    cmake -DCMAKE_CXX_COMPILER=icpx ..
    make -j20
    ../run.sh 

    == Infer 147, time = 1840 micro sec.
    == Infer 148, time = 1841 micro sec.
    == Infer 149, time = 1844 micro sec.


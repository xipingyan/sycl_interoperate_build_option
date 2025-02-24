# How to build oneDNN (OCL runtime)

    cd sycl_interoperate_build_option/onednn_kernel/onednn_ocl/
    
    git clone https://github.com/oneapi-src/oneDNN.git
    cd oneDNN
    git checkout e8901c56a3358883ac34fb15f1c5adc7695ee987

    mkdir build && cd build
    source /opt/intel/oneapi/setvars.sh
    cmake -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx -DDNNL_GPU_RUNTIME=SYCL -DDNNL_CPU_RUNTIME=SYCL -DCMAKE_INSTALL_PREFIX=install ..
    make -j20 && make install

# Build and Run APP

    cd sycl_interoperate_build_option/onednn_kernel/onednn_sycl
    mkdir build && cd build
    source /opt/intel/oneapi/setvars.sh
    export LD_LIBRARY_PATH=../oneDNN/build/install/lib:$LD_LIBRARY_PATH
    cmake -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx ..
    make -j20

    export LD_LIBRARY_PATH=../oneDNN/build/install/lib:$LD_LIBRARY_PATH
    ./onednn_kernel
    
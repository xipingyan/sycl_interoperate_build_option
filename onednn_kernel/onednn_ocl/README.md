# How to build oneDNN (OCL runtime)

    cd sycl_interoperate_build_option/onednn_kernel/onednn_ocl/
    git clone https://github.com/oneapi-src/oneDNN.git
    cd oneDNN
    git checkout e8901c56a3358883ac34fb15f1c5adc7695ee987
    mkdir build && cd build
    source /opt/intel/oneapi/setvars.sh
    cmake -DDNNL_GPU_RUNTIME=OCL -DCMAKE_INSTALL_PREFIX=install ..

# Build and Run APP

    cd sycl_interoperate_build_option/onednn_kernel/onednn_ocl
    mkdir build && cd build
    cmake ..
    ./onednn_kernel

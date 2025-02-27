# ONEDNN GEMM(MatMul)

Verify performance oneDNN GEMM for SYCL and OpenCL runtime.

``onednn_ocl: `` Build oneDNN with OpenCL runtime; <br>
``onednn_sycl: `` Build oneDNN with SYCL runtime;  <br>

# Build app

#### SYCL RUNTIME and Run APP

    cd sycl_interoperate_build_option/onednn_kernel/
    build.sh sycl
    run.sh sycl

#### OCL RUNTIME

    cd sycl_interoperate_build_option/onednn_kernel/
    build.sh ocl
    run.sh ocl


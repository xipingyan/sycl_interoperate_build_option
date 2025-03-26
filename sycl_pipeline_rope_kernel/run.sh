source /opt/intel/oneapi/setvars.sh

export OV_DEVICE='GPU'

# Test: model_gemm.py
# logs_dir=vtune_log_dir_sycl_pipeline_rope_kernel
# rm -rf $logs_dir
# mkdir -p $logs_dir
#OV_GPU_Verbose=4 ONEDNN_VERBOSE=1 

# SYCL_KERNEL=1 
# PERFORMANCE=1 ./sycl_pipeline_rope_kernel
# PERFORMANCE=1  onetrace --chrome-call-logging --chrome-device-timeline ./sycl_pipeline_rope_kernel

# Test compare OpenCL/L0 benchend's performance.
# =====================================================================================================
# PERFORMANCE=1 SYCL_KERNEL=1 ./sycl_pipeline_rope_kernel
# PERFORMANCE=1 SYCL_KERNEL=1 ONEAPI_DEVICE_SELECTOR=OPENCL:GPU ./sycl_pipeline_rope_kernel

 KERNEL_LOOP=100000 PERFORMANCE=1 SYCL_KERNEL=1 ./sycl_pipeline_rope_kernel
 KERNEL_LOOP=100000 PERFORMANCE=1 SYCL_KERNEL=1 ONEAPI_DEVICE_SELECTOR=OPENCL:GPU ./sycl_pipeline_rope_kernel

  KERNEL_LOOP=100000 PERFORMANCE=1 ./sycl_pipeline_rope_kernel
 KERNEL_LOOP=100000 PERFORMANCE=1  ONEAPI_DEVICE_SELECTOR=OPENCL:GPU ./sycl_pipeline_rope_kernel


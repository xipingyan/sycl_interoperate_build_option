source /opt/intel/oneapi/setvars.sh

export OV_DEVICE='GPU'

# Test: model_gemm.py
# logs_dir=vtune_log_dir_sycl_pipeline_rope_kernel
# rm -rf $logs_dir
# mkdir -p $logs_dir
# OV_GPU_Verbose=4 ONEDNN_VERBOSE=1 

PERFORMANCE=1 ./ocl_pipeline_sdpa_kernel

# onetrace profiling
# PERFORMANCE=1  onetrace --chrome-call-logging --chrome-device-timeline ./ocl_pipeline_sdpa_kernel
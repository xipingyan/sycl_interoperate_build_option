source /opt/intel/oneapi/setvars.sh

export OV_DEVICE='GPU'

# Test: model_gemm.py
logs_dir=vd_sycl_sdpa_kernel
rm -rf $logs_dir
mkdir -p $logs_dir
# OV_GPU_Verbose=4 ONEDNN_VERBOSE=1 

PERFORMANCE=1 vtune -collect gpu-hotspots -r $logs_dir -knob gpu-sampling-interval=0.1 -- ./sycl_pipeline_sdpa_kernel
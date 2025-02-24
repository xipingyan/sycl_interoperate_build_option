
source /opt/intel/oneapi/setvars.sh

# export LD_LIBRARY_PATH=/mnt/xiping/gpu_profiling/sycl_interoperate_build_option/onednn_kernel/onednn_sycl/oneDNN/build/install/lib:$LD_LIBRARY_PATH

logs_dir=vtune_log_dir_sycl_onednn

rm -rf $logs_dir
mkdir -p $logs_dir
# OV_GPU_Verbose=4 ONEDNN_VERBOSE=1 
vtune -collect gpu-hotspots -r $logs_dir -knob gpu-sampling-interval=0.1 -- ./onednn_kernel
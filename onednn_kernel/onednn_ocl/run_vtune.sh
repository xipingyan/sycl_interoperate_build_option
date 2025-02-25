
#source /opt/intel/oneapi/vtune/latest/vtune-vars.sh
source /opt/intel/oneapi/vtune/2025.1_internal/vtune-vars.sh

export LD_LIBRARY_PATH=/mnt/xiping/gpu_profiling/sycl_interoperate_build_option/onednn_kernel/onednn_ocl/oneDNN/build/install/lib:$LD_LIBRARY_PATH

logs_dir=vtune_log_dir_ocl_onednn

rm -rf $logs_dir
mkdir -p $logs_dir
# OV_GPU_Verbose=4 ONEDNN_VERBOSE=1 
vtune -collect gpu-hotspots -r $logs_dir -knob gpu-sampling-interval=0.1 -- ./onednn_kernel


source /opt/intel/oneapi/vtune/latest/vtune-vars.sh

logs_dir=vtune_log_dir_ocl_onednn

rm -rf $logs_dir
mkdir -p $logs_dir
# OV_GPU_Verbose=4 ONEDNN_VERBOSE=1 
vtune -collect gpu-hotspots -r $logs_dir -knob gpu-sampling-interval=0.1 -- ./onednn_gemm
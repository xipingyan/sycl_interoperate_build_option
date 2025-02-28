#!/bin/bash
cd $(dirname "$0")

if [[ -z $1 ]];
then 
    echo "Please specific vtune profiling sycl or ocl. For exmaple:"
    echo "$ vtune.sh sycl"
    echo "$ vtune.sh ocl"
    exit
fi
build_option=$1

loop_num=100
if [ ! -z $2 ];
then 
    loop_num=$2
fi
echo "== loop_num=${loop_num}"

source /opt/intel/oneapi/setvars.sh
source /opt/intel/oneapi/vtune/2025.1_internal/vtune-vars.sh

option_sycl="sycl"
if [ "$build_option" == "$option_sycl" ];
then
    logs_dir=vtune_dir_sycl_onednn
    rm -rf $logs_dir
    mkdir -p $logs_dir

    export LD_LIBRARY_PATH=`pwd`/onednn_sycl/oneDNN/build/install/lib:$LD_LIBRARY_PATH
    # ONEDNN_VERBOSE=1 SYCL_UR_TRACE=3 
    vtune -collect gpu-hotspots -r $logs_dir -knob gpu-sampling-interval=0.1 -- ./build_sycl/onednn_kernel_sycl $loop_num
else
    logs_dir=vtune_dir_ocl_onednn
    rm -rf $logs_dir
    mkdir -p $logs_dir

    export LD_LIBRARY_PATH=`pwd`/onednn_ocl/oneDNN/build/install/lib:$LD_LIBRARY_PATH
    # ONEDNN_VERBOSE=1
    vtune -collect gpu-hotspots -r $logs_dir -knob gpu-sampling-interval=0.1 -- ./build_ocl/onednn_kernel_ocl $loop_num
fi
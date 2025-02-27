#!/bin/bash
cd $(dirname "$0")

if [[ -z $1 ]];
then 
    echo "Please specific run sycl or ocl. For exmaple:"
    echo "$ run.sh sycl"
    echo "$ run.sh ocl"
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

option_sycl="sycl"
if [ "$build_option" == "$option_sycl" ];
then
    export LD_LIBRARY_PATH=`pwd`/onednn_sycl/oneDNN/build/install/lib:$LD_LIBRARY_PATH
    # ONEDNN_VERBOSE=1 SYCL_UR_TRACE=3 
    ./build_sycl/onednn_kernel_sycl $loop_num
else
    export LD_LIBRARY_PATH=`pwd`/onednn_ocl/oneDNN/build/install/lib:$LD_LIBRARY_PATH
    ONEDNN_VERBOSE=1 ./build_ocl/onednn_kernel_ocl  $loop_num
fi
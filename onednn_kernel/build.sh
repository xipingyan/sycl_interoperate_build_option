#!/bin/bash
cd $(dirname "$0")

if [[ -z $1 ]];
then 
    echo "Please specific build sycl or ocl. For exmaple:"
    echo "$ build.sh sycl"
    echo "$ build.sh ocl"
    exit
fi

build_option=$1
echo "Start to cmake. build_option = $build_option"

source /opt/intel/oneapi/setvars.sh

option_sycl="sycl"
if [ "$build_option" == "$option_sycl" ];
then
    mkdir -p build_sycl && cd build_sycl
    export LD_LIBRARY_PATH=../onednn_sycl/oneDNN/build/install/lib:$LD_LIBRARY_PATH
    cmake -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx -DGPU_RUNTIME=SYCL ..
else
    mkdir -p build_ocl && cd build_ocl
    export LD_LIBRARY_PATH=../onednn_ocl/oneDNN/build/install/lib:$LD_LIBRARY_PATH
    cmake -DGPU_RUNTIME=OCL ..
fi

echo "Start to make."
make -j 32
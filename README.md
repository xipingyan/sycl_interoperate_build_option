# sycl_interoperate_build_option
Verify build options work or not for SYCL interoperate with OpenCL kernel source.

# Build

    Build ocl_pipeline_rope_kernel based on ocl_pipeline_rope_kernel/README.md
    Build sycl_pipeline_rope_kernel based on sycl_pipeline_rope_kernel/README.md

# Run OpenCL

    <!-- Get OpenCL pipeline result -->
    source /opt/intel/oneapi/setvars.sh 
    cd ocl_pipeline_rope_kernel/build/

    USE_OPTION=1 ./ocl_pipeline_rope_kernel result_ocl_option.log
    ./ocl_pipeline_rope_kernel result_ocl.log

#### Compare result

    cd sycl_interoperate_build_option

    python cmp_2_fn.py ./ocl_pipeline_rope_kernel/build/result_ocl.log ./ocl_pipeline_rope_kernel/build/result_ocl_option.log 

    ...
    @- Line-5351 0.079956
    #+ Line-5351 0.079895

    @- Line-5354 2.851562
    #+ Line-5354 2.849609

    Total diff: 410 / 5377

# Run SYCL

    cd sycl_interoperate_build_option/sycl_pipeline_rope_kernel/build
    USE_OPTION=1 ./sycl_pipeline_rope_kernel > result_sycl_option.log
    ./sycl_pipeline_rope_kernel > result_sycl.log

    <!-- Convert and sort print output -->
    python ../../cvt.py result_sycl_option.log ./result_sycl_option_cvted.log
    python ../../cvt.py result_sycl.log ./result_sycl_cvted.log

#### Compare sycl result

    python ../../cmp_2_fn.py ./result_sycl_option_cvted.log result_sycl_cvted.log

    Total diff: 0 / 5377         Means: enable and not enable build option have same result.

# Compare OpenCL and Sycl Result

    cd sycl_interoperate_build_option

    <!-- Sycl without option VS Ocl with option -->
    python cmp_2_fn.py sycl_pipeline_rope_kernel/build/result_sycl_cvted.log ocl_pipeline_rope_kernel/build/result_ocl_option.log 

    @- Line-5351 0.079956
    #+ Line-5351 0.079895

    @- Line-5354 2.851562
    #+ Line-5354 2.849609

    Total diff: 410 / 5377

    <!-- Sycl without option VS Ocl without option -->
    python cmp_2_fn.py sycl_pipeline_rope_kernel/build/result_sycl_cvted.log ocl_pipeline_rope_kernel/build/result_ocl.log 

    Total diff: 0 / 5377  Means: they are same.

``Conclusion``:  So I guess build option doesn't work.

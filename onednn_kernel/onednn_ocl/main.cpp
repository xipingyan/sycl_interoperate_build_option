
#include <fstream>
#include <thread>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstring>

#include "oneapi/dnnl/dnnl.hpp"
#include "oneapi/dnnl/dnnl_types.h"

inline dnnl::memory::dim product(const dnnl::memory::dims &dims) {
    return std::accumulate(dims.begin(), dims.end(), (dnnl::memory::dim)1,
            std::multiplies<dnnl::memory::dim>());
}

void test_onednn_gemm_ocl(bool runtime_ocl,
                          bool test_performance)
{
    dnnl::engine::kind engine_kind = dnnl::engine::kind::gpu;

    // Create execution dnnl::engine.
    dnnl::engine engine(engine_kind, 0);

    // Create dnnl::stream.
    dnnl::stream engine_stream(engine);

    // onednn_verbose,v1,primitive,exec,gpu:0,inner_product,jit:gemm:any,
    // forward_inference,
    // src:f16::blocked:abcd::f0 
    // wei:f16:a:blocked:abcd::f0 
    // bia:undef::undef::: 
    // dst:f32::blocked:ab::f0,
    // attr-scratchpad:user,,
    // mb8ic896ih1iw1oc151936,
    // 21.5701

    dnnl::memory::dims src_dims = {8, 896};
    dnnl::memory::dims wei_dims = {896, 65536/2}; // 151936, 4096
    dnnl::memory::dims dst_dims = {8, 65536/2};

    // Allocate buffers.
    std::vector<ushort> src_data(product(src_dims));
    std::vector<ushort> weights_data(product(wei_dims));
    std::vector<ushort> dst_data(product(dst_dims));

    // Initialize src, weights, and dst tensors.
    std::generate(src_data.begin(), src_data.end(), []() {
        static int i = 0;
        return (ushort)(std::cos(i++ / 10.f) * 10);
    });
    std::generate(weights_data.begin(), weights_data.end(), []() {
        static int i = 0;
        return (ushort)(std::sin(i++ * 2.f) * 10);
    });

    dnnl::memory::desc src_desc = dnnl::memory::desc(src_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::any);
    dnnl::memory::desc weights_desc = dnnl::memory::desc(wei_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::any);
    dnnl::memory::desc dst_desc = dnnl::memory::desc(dst_dims, dnnl::memory::data_type::f32, dnnl::memory::format_tag::any);

    dnnl::memory::desc usr_src_desc = dnnl::memory::desc(src_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::ab);
    dnnl::memory::desc usr_weights_desc = dnnl::memory::desc(wei_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::ab);
    dnnl::memory::desc usr_dst_desc = dnnl::memory::desc(dst_dims, dnnl::memory::data_type::f32, dnnl::memory::format_tag::ab);

    std::cout << "== prepare primitive_desc." << std::endl;
    auto mm_pd = dnnl::matmul::primitive_desc(engine, src_desc, weights_desc, dst_desc);
    std::cout << "== construct matmul primitive." << std::endl;
    auto mm_prim = dnnl::matmul(mm_pd);

    std::cout << "== prepare user memory." << std::endl;
    auto user_src_mem = dnnl::memory(usr_src_desc, engine);
    auto user_weights_mem = dnnl::memory(usr_weights_desc, engine);
    auto user_dst_mem = dnnl::memory(usr_dst_desc, engine);

    void *mapped_ptr = user_src_mem.map_data();
    if (mapped_ptr)
        std::memcpy(mapped_ptr, src_data.data(), user_src_mem.get_desc().get_size());
    user_src_mem.unmap_data(mapped_ptr);
    mapped_ptr = user_weights_mem.map_data();
    if (mapped_ptr)
        std::memcpy(mapped_ptr, weights_data.data(), user_weights_mem.get_desc().get_size());
    user_weights_mem.unmap_data(mapped_ptr);

    std::cout << "== user to mm memory." << std::endl;
    auto mm_src_mem = user_src_mem;
    if (mm_pd.src_desc() != user_src_mem.get_desc())
    {
        mm_src_mem = dnnl::memory(mm_pd.src_desc(), engine);
        dnnl::reorder(user_src_mem, mm_src_mem)
            .execute(engine_stream, user_src_mem, mm_src_mem);
    }
    auto mm_weights_mem = user_weights_mem;
    if (mm_pd.weights_desc() != user_weights_mem.get_desc())
    {
        mm_weights_mem = dnnl::memory(mm_pd.weights_desc(), engine);
        dnnl::reorder(user_weights_mem, mm_weights_mem)
            .execute(engine_stream, user_weights_mem, mm_weights_mem);
    }
    auto mm_dst_mem = user_dst_mem;
    if (mm_pd.dst_desc() != user_dst_mem.get_desc())
    {
        mm_dst_mem = dnnl::memory(mm_pd.dst_desc(), engine);
        dnnl::reorder(user_dst_mem, mm_dst_mem)
            .execute(engine_stream, user_dst_mem, mm_dst_mem);
    }

    std::cout << "== prepare mm_args." << std::endl;
    std::unordered_map<int, dnnl::memory> mm_args;
    mm_args.insert({DNNL_ARG_SRC, mm_src_mem});
    mm_args.insert({DNNL_ARG_WEIGHTS, mm_weights_mem});
    mm_args.insert({DNNL_ARG_DST, mm_dst_mem});

    std::cout << "== statrt execute." << std::endl;
    for (int i = 0; i < 10; i++)
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        mm_prim.execute(engine_stream, mm_args);
        engine_stream.wait();
        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "  == " << i << ", tm = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " micro sec" << std::endl;
    }
}

int main(int argc, char** argv)
{
    std::cout << "== $ onednn_gemm -h, show help message and exit app." << std::endl;
    if (argc >= 2 && std::string(argv[1]) == std::string("-h")) {       
        std::cout << "  argv[1]:  ocl: means test OpenCL runtime; sycl: means SYCL runtime. Default 1" << std::endl;
        std::cout << "  argv[2]:  acc: means test accuracy. Default test performance." << std::endl;
        exit(-1);
    }

    bool runtime_ocl = true;
    bool test_performance = true;
    if (argc >= 2 && std::string(argv[1]) == std::string("sycl"))
    {
        runtime_ocl = false;
    }
    if (argc >= 3 && std::string(argv[2]) == std::string("acc"))
    {
        test_performance = false;
    }
    std::cout << "  == runtime_ocl = " << runtime_ocl << std::endl;
    std::cout << "  == test_performance = " << test_performance << std::endl;

    std::cout << "== Start test onednn gemm. " << std::endl;
    test_onednn_gemm_ocl(runtime_ocl, test_performance);
    
    return 0;
}
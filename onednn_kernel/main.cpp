
#include <fstream>
#include <thread>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstring>

#include "oneapi/dnnl/dnnl.hpp"
#include "oneapi/dnnl/dnnl_debug.h"

#ifdef GPU_RUNTIME_SYCL
#include <sycl/sycl.hpp>
#include <sycl/ext/oneapi/backend/level_zero.hpp>
namespace syclex = sycl::ext::oneapi::experimental;
#include "oneapi/dnnl/dnnl_sycl.hpp"
#endif


#define M 8
#define K 896
#define N 151936 // 151936, 65536/2

inline dnnl::memory::dim product(const dnnl::memory::dims &dims) {
    return std::accumulate(dims.begin(), dims.end(), (dnnl::memory::dim)1,
            std::multiplies<dnnl::memory::dim>());
}

void test_onednn_kernel(dnnl::engine engine, dnnl::stream engine_stream, int loop_num)
{
    dnnl::memory::dims src_dims = {M, K};
    dnnl::memory::dims wei_dims = {N, K};
    dnnl::memory::dims dst_dims = {M, N};

    // Allocate buffers.
    std::vector<ushort> src_data(product(src_dims));
    std::vector<ushort> weights_data(product(wei_dims));
    std::vector<float> dst_data(product(dst_dims));

    // Initialize src, weights, and dst tensors.
    std::generate(src_data.begin(), src_data.end(), []()
                  {
        static int i = 0;
        return (ushort)(std::cos(i++ / 10.f) * 10); });
    std::generate(weights_data.begin(), weights_data.end(), []()
                  {
        static int i = 0;
        return (ushort)(std::sin(i++ * 2.f) * 10); });

    #define DEBUG_POS std::cout << "Line:" << __LINE__ << std::endl
    DEBUG_POS;
    dnnl::memory::desc src_desc = dnnl::memory::desc(src_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::any);
    dnnl::memory::desc weights_desc = dnnl::memory::desc(wei_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::any);
    dnnl::memory::desc dst_desc = dnnl::memory::desc(dst_dims, dnnl::memory::data_type::f32, dnnl::memory::format_tag::any);
    DEBUG_POS;

    dnnl::memory::desc usr_src_desc = dnnl::memory::desc(src_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::ab);
    dnnl::memory::desc usr_weights_desc = dnnl::memory::desc(wei_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::ab);
    DEBUG_POS;
    dnnl::memory::desc usr_dst_desc = dnnl::memory::desc(dst_dims, dnnl::memory::data_type::f32, dnnl::memory::format_tag::ab);

    std::cout << "== prepare primitive_desc." << std::endl;
    auto ip_pd = dnnl::inner_product_forward::primitive_desc(engine, dnnl::prop_kind::forward_inference, src_desc, weights_desc, dst_desc);

    std::cout << "== construct inner_product_forward primitive." << std::endl;
    auto ip_prim = dnnl::inner_product_forward(ip_pd);

    std::cout << "== prepare user memory." << std::endl;
    auto user_src_mem = dnnl::memory(usr_src_desc, engine);
    auto user_weights_mem = dnnl::memory(usr_weights_desc, engine);
    auto user_dst_mem = dnnl::memory(usr_dst_desc, engine);

    // Copy data from host to device.
#ifdef GPU_RUNTIME_SYCL
    sycl::queue queue = dnnl::sycl_interop::get_queue(engine_stream);
    queue.memcpy(user_src_mem.map_data(), src_data.data(), user_src_mem.get_desc().get_size()).wait();
    queue.memcpy(user_weights_mem.map_data(), weights_data.data(), user_weights_mem.get_desc().get_size()).wait();
#else
#endif

    std::cout << "== user to mm memory." << std::endl;
    auto ip_src_mem = user_src_mem;
    if (ip_pd.src_desc() != user_src_mem.get_desc())
    {
        ip_src_mem = dnnl::memory(ip_pd.src_desc(), engine);
        dnnl::reorder(user_src_mem, ip_src_mem)
            .execute(engine_stream, user_src_mem, ip_src_mem);
    }
    auto ip_weights_mem = user_weights_mem;
    if (ip_pd.weights_desc() != user_weights_mem.get_desc())
    {
        ip_weights_mem = dnnl::memory(ip_pd.weights_desc(), engine);
        dnnl::reorder(user_weights_mem, ip_weights_mem)
            .execute(engine_stream, user_weights_mem, ip_weights_mem);
    }
    auto ip_dst_mem = user_dst_mem;
    if (ip_pd.dst_desc() != user_dst_mem.get_desc())
    {
        ip_dst_mem = dnnl::memory(ip_pd.dst_desc(), engine);
        dnnl::reorder(user_dst_mem, ip_dst_mem)
            .execute(engine_stream, user_dst_mem, ip_dst_mem);
    }

    std::cout << "== prepare ip_args." << std::endl;
    std::unordered_map<int, dnnl::memory> ip_args;
    ip_args.insert({DNNL_ARG_SRC, ip_src_mem});
    ip_args.insert({DNNL_ARG_WEIGHTS, ip_weights_mem});
    ip_args.insert({DNNL_ARG_DST, ip_dst_mem});

    std::cout << "== statrt execute." << std::endl;
    for (int i = 0; i < loop_num; i++)
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        ip_prim.execute(engine_stream, ip_args);
        engine_stream.wait();
        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "  == " << i << ", tm = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " micro sec" << std::endl;
    }

    std::cout << "== copy dst buf." << std::endl;
#ifdef GPU_RUNTIME_SYCL
    queue.memcpy(dst_data.data(), ip_dst_mem.map_data(), ip_dst_mem.get_desc().get_size()).wait();
    std::string result_fn = "result_sycl.log";
#else
    std::string result_fn = "result_ocl.log";
#endif

    std::cout << "== Start to save result to " << result_fn << std::endl;
    std::filebuf fb;
    fb.open(result_fn.c_str(), std::ios::out);
    std::ostream os(&fb);
    os.precision(6);
    for (size_t i = 0; i < dst_data.size(); i++)
    {
        os << std::fixed << dst_data[i] << std::endl;
    }
    fb.close();
}

int main(int argc, char** argv)
{
    int loop_num = 100;
    if (argc == 2)
    {
        loop_num = std::atoi(argv[1]);
    }

#ifdef GPU_RUNTIME_SYCL
    auto queue = sycl::queue(sycl::gpu_selector_v);
    std::cout << "  == Using "
              << queue.get_device().get_info<sycl::info::device::name>()
              << ", Backend: " << queue.get_backend()
              << std::endl;

    auto engine = dnnl::sycl_interop::make_engine(queue.get_device(), queue.get_context());
    auto engine_stream = dnnl::sycl_interop::make_stream(engine, queue);
#else
    dnnl::engine::kind engine_kind = dnnl::engine::kind::gpu;

    // Create execution dnnl::engine.
    dnnl::engine engine(engine_kind, 0);

    // Create dnnl::stream.
    dnnl::stream engine_stream(engine);
#endif

    std::cout << "  == loop_num: " << loop_num << std::endl;
    test_onednn_kernel(engine, engine_stream, loop_num);

    return 0;
}
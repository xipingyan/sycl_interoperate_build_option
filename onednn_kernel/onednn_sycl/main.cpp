
#include <fstream>
#include <thread>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstring>

#include <sycl/sycl.hpp>
#include <sycl/ext/oneapi/backend/level_zero.hpp>
namespace syclex = sycl::ext::oneapi::experimental;

#include "oneapi/dnnl/dnnl.hpp"
#include "oneapi/dnnl/dnnl_debug.h"
#include "oneapi/dnnl/dnnl_sycl.hpp"


#define M 8
#define K 896
#define N 65536/2

inline dnnl::memory::dim product(const dnnl::memory::dims &dims) {
    return std::accumulate(dims.begin(), dims.end(), (dnnl::memory::dim)1,
            std::multiplies<dnnl::memory::dim>());
}

void test_onednn_kernel(dnnl::engine engine, dnnl::stream engine_stream, sycl::queue& queue)
{
    dnnl::memory::dims src_dims = {M, K};
    dnnl::memory::dims wei_dims = {K, N};
    dnnl::memory::dims dst_dims = {M, N};

    // Allocate buffers.
    std::vector<ushort> src_data(product(src_dims));
    std::vector<ushort> weights_data(product(wei_dims));
    std::vector<ushort> dst_data(product(dst_dims));

    // Initialize src, weights, and dst tensors.
    std::generate(src_data.begin(), src_data.end(), []()
                  {
        static int i = 0;
        return (ushort)(std::cos(i++ / 10.f) * 10); });
    std::generate(weights_data.begin(), weights_data.end(), []()
                  {
        static int i = 0;
        return (ushort)(std::sin(i++ * 2.f) * 10); });

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

    // Copy data from host to device.
    queue.memcpy(user_src_mem.map_data(), src_data.data(), user_src_mem.get_desc().get_size()).wait();
    queue.memcpy(user_weights_mem.map_data(), weights_data.data(), user_weights_mem.get_desc().get_size()).wait();

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
    // auto sycl_queue
    //                     = dnnl::sycl_interop::get_queue(dnnl::stream(eng));
    //             sycl_queue.memcpy(dst_ptr, handle, size).wait();
    auto queue = sycl::queue(sycl::gpu_selector_v);
    std::cout << "  == Using "
              << queue.get_device().get_info<sycl::info::device::name>()
              << ", Backend: " << queue.get_backend()
              << std::endl;

    auto engine = dnnl::sycl_interop::make_engine(queue.get_device(), queue.get_context());
    auto engine_stream = dnnl::sycl_interop::make_stream(engine, queue);

    test_onednn_kernel(engine, engine_stream, queue);

    return 0;
}
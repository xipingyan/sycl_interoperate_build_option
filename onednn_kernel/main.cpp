
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


typedef unsigned short ushort;
typedef unsigned int uint;

inline uint as_uint(const float x)
{
    return *(uint *)&x;
}
inline float as_float(const uint x)
{
    return *(float *)&x;
}

inline float half_to_float(const ushort x)
{                                                                                                                                                        // IEEE-754 16-bit floating-point format (without infinity): 1-5-10, exp-15, +-131008.0, +-6.1035156E-5, +-5.9604645E-8, 3.311 digits
    const uint e = (x & 0x7C00) >> 10;                                                                                                                   // exponent
    const uint m = (x & 0x03FF) << 13;                                                                                                                   // mantissa
    const uint v = as_uint((float)m) >> 23;                                                                                                              // evil log2 bit hack to count leading zeros in denormalized format
    return as_float((x & 0x8000) << 16 | (e != 0) * ((e + 112) << 23 | m) | ((e == 0) & (m != 0)) * ((v - 37) << 23 | ((m << (150 - v)) & 0x007FE000))); // sign : normalized : denormalized
}

inline ushort float_to_half(const float x)
{                                                                                                                                                                                       // IEEE-754 16-bit floating-point format (without infinity): 1-5-10, exp-15, +-131008.0, +-6.1035156E-5, +-5.9604645E-8, 3.311 digits
    const uint b = as_uint(x) + 0x00001000;                                                                                                                                             // round-to-nearest-even: add last bit after truncated mantissa
    const uint e = (b & 0x7F800000) >> 23;                                                                                                                                              // exponent
    const uint m = b & 0x007FFFFF;                                                                                                                                                      // mantissa; in line below: 0x007FF000 = 0x00800000-0x00001000 = decimal indicator flag - initial rounding
    return (b & 0x80000000) >> 16 | (e > 112) * ((((e - 112) << 10) & 0x7C00) | m >> 13) | ((e < 113) & (e > 101)) * ((((0x007FF000 + m) >> (125 - e)) + 1) >> 1) | (e > 143) * 0x7FFF; // sign : normalized : denormalized : saturate
}

inline dnnl::memory::dim product(const dnnl::memory::dims &dims) {
    return std::accumulate(dims.begin(), dims.end(), (dnnl::memory::dim)1,
            std::multiplies<dnnl::memory::dim>());
}

inline void init_inputs(std::vector<ushort> &src_data, std::vector<ushort> &weights_data)
{
    std::generate(src_data.begin(), src_data.end(), []()
                  {
        static int i = 0;
        return float_to_half(std::cos(i++ / 3.2f)); });
    std::generate(weights_data.begin(), weights_data.end(), []()
                  {
        static int i = 0;
        return float_to_half(std::sin(i++ * 3.2f)); });

    std::cout << "weights_data[20] = " << half_to_float(weights_data[20]) << std::endl;
}

dnnl::memory reorder_memory(dnnl::memory &user_mem, const dnnl::memory::desc& pd, dnnl::engine &engine, dnnl::stream &engine_stream)
{
    auto mem = dnnl::memory(pd, engine);
    dnnl::reorder(user_mem, mem)
        .execute(engine_stream, user_mem, mem);
    return mem;
}

void test_onednn_kernel(dnnl::engine engine, dnnl::stream engine_stream, int loop_num)
{
    dnnl::memory::dims src_dims = {M, K, 1, 1};
    dnnl::memory::dims wei_dims = {N, K, 1, 1};
    dnnl::memory::dims dst_dims = {M, N};

    // Allocate buffers.
    std::vector<ushort> src_data(product(src_dims));
    std::vector<ushort> weights_data(product(wei_dims));
    std::vector<float> dst_data(product(dst_dims));

    // Initialize src, weights, and dst tensors.
    init_inputs(src_data, weights_data);

    // onednn_verbose,v1,primitive,exec,gpu:0,inner_product,jit:gemm:any,forward_inference,
    // src:f16::blocked:abcd::f0 
    // wei:f16:a:blocked:abcd::f0 bia:undef::undef::: 
    // dst:f32::blocked:ab::f0,attr-scratchpad:user,,
    // mb8ic896ih1iw1oc151936,21.5701

    dnnl::memory::desc src_desc = dnnl::memory::desc(src_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::any);
    dnnl::memory::desc weights_desc = dnnl::memory::desc(wei_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::any);
    dnnl::memory::desc dst_desc = dnnl::memory::desc(dst_dims, dnnl::memory::data_type::f32, dnnl::memory::format_tag::any);

    dnnl::memory::desc usr_src_desc = dnnl::memory::desc(src_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::abcd);
    dnnl::memory::desc usr_weights_desc = dnnl::memory::desc(wei_dims, dnnl::memory::data_type::f16, dnnl::memory::format_tag::abcd);
    dnnl::memory::desc usr_dst_desc = dnnl::memory::desc(dst_dims, dnnl::memory::data_type::f32, dnnl::memory::format_tag::ab);

    std::cout << "== prepare primitive_desc." << std::endl;
    dnnl::primitive_attr prim_attr;
    prim_attr.set_scratchpad_mode(dnnl::scratchpad_mode::user);
    auto ip_pd = dnnl::inner_product_forward::primitive_desc(engine,
                                                             dnnl::prop_kind::forward_inference,
                                                             src_desc, weights_desc, dst_desc, prim_attr);

    std::cout << "== construct inner_product_forward primitive." << std::endl;
    auto ip_prim = dnnl::inner_product_forward(ip_pd);

    std::cout << "== prepare user memory." << std::endl;
    auto user_src_mem = dnnl::memory(usr_src_desc, engine);
    auto user_weights_mem = dnnl::memory(usr_weights_desc, engine);
    auto user_dst_mem = dnnl::memory(usr_dst_desc, engine);

    // Copy data from user.
#ifdef GPU_RUNTIME_SYCL
    sycl::queue queue = dnnl::sycl_interop::get_queue(engine_stream);
#endif
    void *mapped_ptr = user_src_mem.map_data();
    if (mapped_ptr)
        std::memcpy(mapped_ptr, src_data.data(), user_src_mem.get_desc().get_size());
    user_src_mem.unmap_data(mapped_ptr);
    mapped_ptr = user_weights_mem.map_data();
    if (mapped_ptr)
        std::memcpy(mapped_ptr, weights_data.data(), user_weights_mem.get_desc().get_size());
    user_weights_mem.unmap_data(mapped_ptr);

    std::cout << "== primitive desc reorder." << std::endl;
    auto ip_src_mem = reorder_memory(user_src_mem, ip_pd.src_desc(), engine, engine_stream);
    auto ip_weights_mem = reorder_memory(user_weights_mem, ip_pd.weights_desc(), engine, engine_stream);

#ifdef GPU_RUNTIME_SYCL
    auto* dst_usm_host = sycl::malloc_host(user_dst_mem.get_desc().get_size(),
                                          dnnl::sycl_interop::get_context(engine));
    auto* dst_usm_dev = sycl::malloc_device(user_dst_mem.get_desc().get_size(),
                                           dnnl::sycl_interop::get_device(engine), dnnl::sycl_interop::get_context(engine));
    auto ip_dst_mem = dnnl::sycl_interop::make_memory(
        ip_pd.dst_desc(), engine, dnnl::sycl_interop::memory_kind::usm, dst_usm_dev);
#else
    auto ip_dst_mem = reorder_memory(user_dst_mem, ip_pd.dst_desc(), engine, engine_stream);
#endif

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
    float* dst_host = (float*)dst_usm_host;
    queue.memcpy(dst_usm_host, dst_usm_dev, ip_dst_mem.get_desc().get_size()).wait();
    std::string result_fn = "result_sycl.log";
#else
    float* dst_host = dst_data.data();
    mapped_ptr = ip_dst_mem.map_data();
    if (mapped_ptr)
        std::memcpy(dst_host, mapped_ptr, ip_dst_mem.get_desc().get_size());
    ip_dst_mem.unmap_data(mapped_ptr);
    std::string result_fn = "result_ocl.log";
#endif

    std::cout << "== Start to save result to " << result_fn << std::endl;
    std::filebuf fb;
    fb.open(result_fn.c_str(), std::ios::out);
    std::ostream os(&fb);
    os.precision(5);
    // for (size_t i = 0; i < 3; i++)
    //     std::cout << "dst_host[i]=" << dst_host[i] << std::endl;
    for (size_t i = 0; i < dst_data.size(); i++)
    {
        os << std::fixed << dst_host[i] << std::endl;
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
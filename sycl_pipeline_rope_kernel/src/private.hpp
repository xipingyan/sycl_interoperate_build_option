#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <sycl/sycl.hpp>

// Just verify if opencl and sycl have light different result.
int test_sycl_olc_interoperate_l0_backend_rope_ref();

struct DumpData
{
    std::string format;
    std::vector<int> shape;
    std::vector<float> data;
    std::string to_string();
    sycl::half* to_half(sycl::queue queue);
    int* to_int(sycl::queue queue);
};
DumpData load_dump_data(std::string fn);
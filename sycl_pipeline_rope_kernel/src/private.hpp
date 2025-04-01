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
    sycl::half *to_half(sycl::queue queue);
    int *to_int(sycl::queue queue);
};
DumpData load_dump_data(std::string fn);

inline int64_t calc_kernel_tm(sycl::event event)
{
    auto start_time = event.get_profiling_info<sycl::info::event_profiling::command_start>();
    auto end_time = event.get_profiling_info<sycl::info::event_profiling::command_end>();

    // Calculate execution time in nanoseconds.
    auto execution_time_ns = end_time - start_time;

    // Convert to milliseconds for easier reading.
    auto execution_time_micro_sec = static_cast<int64_t>(execution_time_ns) / 1000.0;
    return execution_time_micro_sec;
}

class CStaticTM
{
private:
    CStaticTM() = default;
    CStaticTM(const CStaticTM &) = delete;
    CStaticTM &operator=(const CStaticTM &) = delete;

public:
    int64_t count = 0;
    int64_t sum = 0;

public:
    static std::shared_ptr<CStaticTM> create()
    {
        return std::shared_ptr<CStaticTM>(new CStaticTM());
    }

    void add(int64_t tm)
    {
        sum += tm;
        count++;
    }

    void print(std::string prefix)
    {
        std::cout << prefix << ", Mean: [" << sum << "/" << count << "] = " << (float)sum / count << " micro sec." << std::endl;
    }

    void print_host_kernel_diff(std::string prefix, std::shared_ptr<CStaticTM> kernel_stm)
    {
        auto kernel_mean = (float)kernel_stm->sum / kernel_stm->count;
        auto host_mean = (float)sum / count;
        std::cout << prefix << " " << host_mean - kernel_mean << " micro sec." << std::endl;
    }
};

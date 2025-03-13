#include "private.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#elif __has_include(<CL/sycl.hpp>)
#include <CL/sycl.hpp>
#else
#error "Unsupported compiler"
#endif

#include <sycl/ext/oneapi/backend/level_zero.hpp>
namespace syclex = sycl::ext::oneapi::experimental;
#include <sycl/ext/oneapi/properties/properties.hpp>

static bool get_env(std::string env)
{
	auto env_str = std::getenv(env.c_str());
	if (env_str && std::string("1") == env_str)
	{
		std::cout << "  == Get: " << env << " = 1" << std::endl;
		return true;
	}
	std::cout << "  == Get: " << env << " = 0" << std::endl;
	return false;
}

static std::string load_kernel(std::string kernel_fn)
{
	std::ifstream kernel_file(kernel_fn.c_str(), std::ios::in | std::ios::binary);
	if (kernel_file.is_open())
	{
		std::string ret;
		auto beg = kernel_file.tellg();
		kernel_file.seekg(0, std::ios::end);
		auto end = kernel_file.tellg();
		kernel_file.seekg(0, std::ios::beg);

		ret.resize((size_t)(end - beg));
		kernel_file.read(&ret[0], (size_t)(end - beg));

		return {std::move(ret)};
	}
	return std::string();
}

static sycl::event launchOpenCLKernelOnline(sycl::queue &q, std::string source,
											std::string func_name, std::vector<std::pair<void *, size_t>> &params,
											sycl::event &dep_event, bool test_performance)
{
	std::cout << "  == Start to kernel_bundle opencl source" << std::endl;
	sycl::kernel_bundle<sycl::bundle_state::ext_oneapi_source> kb_src =
		syclex::create_kernel_bundle_from_source(
			q.get_context(),
			syclex::source_language::opencl,
			source);

	// Compile and link the kernel from the source definition.
	std::cout << "  == Start to kernel_bundle kb_src" << std::endl;
	std::vector<std::string> option_flags = {"-cl-mad-enable", "-cl-std=CL2.0"};

	sycl::kernel_bundle<sycl::bundle_state::executable> kb_exe = syclex::build(kb_src, syclex::properties{syclex::build_options{option_flags}});

	// Get a "kernel" object representing the kernel defined in the
	// source string.
	std::cout << "  == Start to get sycl::kernel" << std::endl;
	sycl::kernel k = kb_exe.ext_oneapi_get_kernel(func_name);

	// global and local range.
	sycl::nd_range ndr = sycl::nd_range{sycl::range{256, 64, 28}, sycl::range{256, 1, 1}};
	// sycl::nd_range ndr = sycl::nd_range{sycl::range{28, 64, 256}, sycl::range{1, 1, 256}};

	std::cout << "  == ndr=["
			  << ndr.get_global_range()[0] << ", " << ndr.get_global_range()[1] << ", " << ndr.get_global_range()[2]
			  << "], local_range=[" << ndr.get_local_range()[0] << ", " << ndr.get_local_range()[1] << ", "
			  << ndr.get_local_range()[2] << "]" << std::endl;

	for (int i = 0; i < params.size(); i++)
	{
		std::cout << "    params[" << i << "] = " << params[i].second << ", addr = " << params[i].first << std::endl;
	}

	std::cout << "  == Start to submit" << std::endl;
	sycl::event ret_ev;
	size_t loop_num = test_performance ? 150 : 1;
	for (size_t i = 0; i < loop_num; i++)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		ret_ev = q.submit([&](sycl::handler &cgh)
						  {
                        cgh.depends_on(dep_event);
						// auto out = sycl::stream(1024, 768, cgh);

						for (int i = 0; i < params.size(); i++)
						{
							cgh.set_arg(i, params[i].first);
						}

						// Invoke the kernel over an nd-range.
						cgh.parallel_for(ndr, k); });
		ret_ev.wait();
		auto t2 = std::chrono::high_resolution_clock::now();
		if (test_performance)
			std::cout << "  == Infer " << i << ", time = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " micro sec." << std::endl;
	}
	return ret_ev;
}

int test_sycl_olc_interoperate_l0_backend()
{
	std::cout << "== Test: " << __FUNCTION__ << ", " << __FILE__ << ":" << __LINE__ << std::endl;

	bool test_performance = get_env("PERFORMANCE");
	bool test_sycl_kernel = get_env("SYCL_KERNEL");
	if (get_env("OCL_KERNEL"))
	{
		test_sycl_kernel = false;
	}
	std::cout << "  == test_sycl_kernel = " << test_sycl_kernel << std::endl;

	// It's hard to read for original cl file.
	// Convert to clean code via:
	// $ cpp original.cl > clean.cl
	std::string kernel_path = "../src/kernel_sdpa/";
	std::string kernel_source = load_kernel(kernel_path + "SYCL_LZ_program_1_bucket_0_part_67_8198298808319190670.cl");

	// std::cout << "  kernel_source = " << kernel_source << std::endl;

	auto queue = sycl::queue(sycl::gpu_selector_v);
	std::cout << "  == Using "
			  << queue.get_device().get_info<sycl::info::device::name>()
			  << ", Backend: " << queue.get_backend()
			  << std::endl;

	sycl::event ev;
	// std::vector<std::pair<sycl::buffer<uint8_t, 1>, bool>> params;
	std::vector<std::pair<void *, size_t>> params;

	//  gws=[28, 64, 256] lws=[1, 1, 256]
	// global_range=[256, 64, 28], local_range=[256, 1, 1]

	DumpData input0;
	input0.data = std::vector<float>({1, 28, 1, 1, 1, 1, 1024, 128,
									  1, 4, 1, 1, 1, 1, 1024, 128,
									  0, 37, 1, 4, 1, 1, 1, 1, 1024, 128,
									  0, 136, 1, 1, 1, 1, 1, 1, 1024, 1024,
									  1, 1, 1, 1, 1, 1, 1024, 1,
									  1, 28, 1, 1, 1, 1, 1024, 128});
	input0.format = "bfyx";
	input0.shape = {52};

	std::string param_name = "program1_network1_0_indirectsdpa___module.model.layers.0.self_attn_aten__scaled_dot_product_attention_ScaledDotProductAttention_";

	auto input1 = load_dump_data(kernel_path + param_name + "src0.txt");
	auto input2 = load_dump_data(kernel_path + param_name + "src1.txt");
	auto input3 = load_dump_data(kernel_path + param_name + "src2.txt");
	auto input4 = load_dump_data(kernel_path + param_name + "src3.txt");
	auto input5 = load_dump_data(kernel_path + param_name + "src4.txt");
	auto output_expected = load_dump_data(kernel_path + param_name + "dst0.txt");

	// auto buf_host_0 = input0.to_int_host(queue);
	auto buf_host_0 = input0.to_int_device(queue);
	auto buf_dev_1 = input1.to_half_device(queue);
	auto buf_dev_2 = input2.to_half_device(queue);
	auto buf_dev_3 = input3.to_half_device(queue);
	auto buf_dev_4 = input4.to_half_device(queue);
	auto buf_dev_5 = input5.to_half_device(queue);
	auto buf_dev_6 = sycl::malloc_device(sizeof(float), queue);
	auto buf_dev_7 = sycl::malloc_device(sizeof(float), queue);
	auto buf_dev_8 = sycl::malloc_device(sizeof(sycl::half), queue);
	queue.fill(buf_dev_6, 0, sizeof(float)).wait();
	queue.fill(buf_dev_7, 0, sizeof(float)).wait();
	queue.fill(buf_dev_8, 0, sizeof(sycl::half)).wait();

	// Device memory will trigger app crash.
	auto output_host = sycl::malloc_host<sycl::half>(output_expected.data.size(), queue);

	std::cout << "  == run OpenCL C kernel." << std::endl;
	params.push_back({buf_host_0, input0.data.size() * sizeof(int)});
	params.push_back({buf_dev_1, input1.data.size() * sizeof(sycl::half)});
	params.push_back({buf_dev_2, input2.data.size() * sizeof(sycl::half)});
	params.push_back({buf_dev_3, input3.data.size() * sizeof(sycl::half)});
	params.push_back({buf_dev_4, input4.data.size() * sizeof(sycl::half)});
	params.push_back({buf_dev_5, input5.data.size() * sizeof(sycl::half)});
	params.push_back({buf_dev_6, sizeof(float)});
	params.push_back({buf_dev_7, sizeof(float)});
	params.push_back({buf_dev_8, sizeof(sycl::half)});

	auto kernel_name = "sdpa_opt_multi_tokens_8660372428234100028_0_0__sa";
	auto ret_ev = launchOpenCLKernelOnline(queue, kernel_source, kernel_name, params, ev, test_performance);
	ret_ev.wait();

	// Device to Host
	// queue.memcpy(output_host, output_dev, output_expected.data.size() * sizeof(sycl::half)).wait();

	// Compare result.
	// ???????????????????????
	// I don't why I can't get result from output_buf(it is shared),
	// because If I replace to some simple kernels, I can get result from output_buf.
	if (!test_performance)
	{
		std::cout << "  == Compare input and output:" << std::endl;
		int diff_num = 0;
		for (size_t i = 0; i < output_expected.data.size(); i++)
		{
			if (fabs(output_expected.data[i] - output_host[i]) > 0.00001f)
			{
				diff_num++;
				std::cout << "    output_expected[" << i << "] = " << output_expected.data[i] << " VS output_host[" << i << "] = " << output_host[i] << std::endl;
			}
			if (diff_num > 5)
			{
				break;
			}
		}

		std::cout << std::endl;
		std::cout << "  == Input and Ouput are " << (diff_num == 0 ? "Same." : "Not Same.") << std::endl;
	}

	sycl::free(buf_host_0, queue);
	// sycl::free(buf_dev_0, queue);
	sycl::free(buf_dev_1, queue);
	sycl::free(buf_dev_2, queue);
	sycl::free(buf_dev_3, queue);
	sycl::free(buf_dev_4, queue);
	sycl::free(buf_dev_5, queue);
	sycl::free(output_host, queue);

	return 0;
}

int test_build_asm() {
#if 0
	std::string fun_name = "igc_check";
	const char *kernel_code = R""""(
        // Kernel name: igc_check
        kernel void igc_check() {
            __asm__ volatile(
                    ".decl AA0 v_type=G type=ud num_elts=1\n"
                    ".decl AA1 v_type=G type=ud num_elts=1\n"
                    ".implicit_PSEUDO_INPUT AA0 offset=256 size=4\n"
                    ".implicit_PSEUDO_INPUT AA1 offset=256 size=4\n"
                    "mov (M1_NM,1) AA0(0,0)<1> AA1(0,0)<0;1,0>\n"
            );
        }
        )"""";
	// global and local range.
	sycl::nd_range ndr = sycl::nd_range{sycl::range{256, 64, 64}, sycl::range{8, 8, 8}};
	std::vector<std::string> option_flags = {};
#else
	std::string fun_name = "is_local_block_io_supported";
	std::string kernel_code =
		"__attribute__((intel_reqd_sub_group_size(8)))"
		"__attribute__((reqd_work_group_size(8, 1, 1)))"
		"void kernel is_local_block_io_supported(global uchar* dst) {"
		"    uint lid = get_sub_group_local_id();"
		"    uchar val = (uchar)lid * 2;"
		"    __local uchar tmp_slm[8];"
		"    intel_sub_group_block_write_uc2(tmp_slm, (uchar2)(val));"
		"    barrier(CLK_LOCAL_MEM_FENCE);"
		"    uchar2 read = intel_sub_group_block_read_uc2(tmp_slm);"
		"    dst[lid] = read.s0 + 1;"
		"}";
	// global and local range.
	sycl::nd_range ndr = sycl::nd_range{sycl::range{256, 64, 64}, sycl::range{1, 1, 8}};

	std::vector<std::string> option_flags = {
		// "-cl-mad-enable", "-cl-std=CL2.0", 
		"-Dcl_intel_subgroup_local_block_io", 
		"-DLOCAL_BLOCK_IO_SUPPORTED=1"
	};
#endif

	auto queue = sycl::queue(sycl::gpu_selector_v);
	std::cout << "  == Using "
			  << queue.get_device().get_info<sycl::info::device::name>()
			  << ", Backend: " << queue.get_backend()
			  << std::endl;

	std::cout << "  == Start to kernel_bundle opencl source" << std::endl;
	sycl::kernel_bundle<sycl::bundle_state::ext_oneapi_source> kb_src =
		syclex::create_kernel_bundle_from_source(
			queue.get_context(),
			syclex::source_language::opencl,
			kernel_code);

	// Compile and link the kernel from the source definition.
	std::cout << "  == Start to kernel_bundle kb_src" << std::endl;

	sycl::kernel_bundle<sycl::bundle_state::executable> kb_exe = syclex::build(kb_src, syclex::properties{syclex::build_options{option_flags}});

	// Get a "kernel" object representing the kernel defined in the
	// source string.
	std::cout << "  == Start to get sycl::kernel" << std::endl;
	sycl::kernel k = kb_exe.ext_oneapi_get_kernel(fun_name);

	std::cout << "  == Start to submit" << std::endl;
	sycl::event ret_ev;
	bool test_performance = true;
	size_t loop_num = test_performance ? 5 : 1;
	for (size_t i = 0; i < loop_num; i++)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		ret_ev = queue.submit([&](sycl::handler &cgh)
						  {
						// Invoke the kernel over an nd-range.
						cgh.parallel_for(ndr, k); });
		ret_ev.wait();
		auto t2 = std::chrono::high_resolution_clock::now();
		if (test_performance)
			std::cout << "  == Infer " << i << ", time = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " micro sec." << std::endl;
	}
	return 0;
}

int main()
{
    std::cout << "== Debug MACRO tip:" << std::endl;
    std::cout << "  Default:        test accuracy only." << std::endl;
    std::cout << "  PERFORMANCE=1:  Test performance only." << std::endl;

	// return test_build_asm();
    test_sycl_olc_interoperate_l0_backend();

    return 0;
}
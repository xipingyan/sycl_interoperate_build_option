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
	sycl::nd_range ndr = sycl::nd_range{sycl::range{192, 14, 1}, sycl::range{192, 2, 1}};
	std::cout << "  == ndr=["
			  << ndr.get_global_range()[0] << ", " << ndr.get_global_range()[1] << ", " << ndr.get_global_range()[2]
			  << "], local_range=[" << ndr.get_local_range()[0] << ", " << ndr.get_local_range()[1] << ", "
			  << ndr.get_local_range()[2] << "]" << std::endl;

	for (int i = 0; i < params.size(); i++)
	{
		std::cout << "    params[" << i << "] = " << params[i].second << ", addr = " << params[i].first << std::endl;
	}

	auto profiling_host_tm = CStaticTM::create();
	auto profiling_kernel_tm = CStaticTM::create();

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
		{
			auto kernel_tm = calc_kernel_tm(ret_ev);
			auto tm = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
			// std::cout << "  == Infer " << i << ", time = " << tm << " micro sec." << std::endl;
			if (i >= 100)
			{
				profiling_host_tm->add(tm);
				profiling_kernel_tm->add(kernel_tm);

				if (loop_num == i + 1)
				{
					profiling_host_tm->print("--->host   tm:");
					profiling_kernel_tm->print("--->kernel tm:");
					profiling_host_tm->print_host_kernel_diff("---> host - kernel tm:", profiling_kernel_tm);
				}
			}
		}
	}
	return ret_ev;
}

// Rewrite OpenCL C kernel based on SYCL.
static sycl::event launchSyclKernel(sycl::queue &q, int *buf0, sycl::half *buf1, sycl::half *buf2, sycl::half *buf3, sycl::half *out_buf, bool test_performance)
{
	sycl::nd_range ndr = sycl::nd_range{sycl::range{1, 14, 192}, sycl::range{1, 2, 192}};
	auto *shape_info = buf0;
	auto *input = buf1;
	auto *cos = buf2;
	auto *sin = buf3;
	auto *output = out_buf;

	size_t kernel_loop_num = 1;
	if (std::getenv("KERNEL_LOOP"))
	{
		kernel_loop_num = std::atoi(std::getenv("KERNEL_LOOP"));
		std::cout << "  == Get: KERNEL_LOOP = " << kernel_loop_num << std::endl;
	}

	auto profiling_host_tm = CStaticTM::create();
	auto profiling_kernel_tm = CStaticTM::create();

	sycl::event ret_ev;
	size_t loop_num = test_performance ? 150 : 1;
	for (size_t i = 0; i < loop_num; i++)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		auto cur_event = q.submit([&](sycl::handler &cgh)
								  { 
					// auto out = sycl::stream(1024, 768, cgh);
					cgh.parallel_for<class rope_sycl>(ndr, [=](sycl::nd_item<3> itm)
													  {
															for (size_t kl = 0; kl < kernel_loop_num; kl++) 
															{
																const uint b = itm.get_global_id(0);
																const uint h = itm.get_global_id(1);
																const uint p = (uint)itm.get_global_id(2) / 32;
																const uint r = (uint)itm.get_global_id(2) % 32;
															
																uint input_idx = ((1 * 0) + (64 * (shape_info[8])) + ((64 * (14 + (shape_info[8] + shape_info[9]))) * 0) + ((64 * (14 + (shape_info[8] + shape_info[9])) * 1) * 0) + ((64 * (14 + (shape_info[8] + shape_info[9])) * 1 * 1 * 1 * 1) * 0) + ((64 * (14 + (shape_info[8] + shape_info[9])) * 1 * 1 * 1 * 1 * (shape_info[1] + 0)) * 0)) + (0) * 1 + (h) * 64 + (p) * (64 * (14 + (shape_info[8] + shape_info[9])) * 1 * 1 * 1 * 1) + (b) * (64 * (14 + (shape_info[8] + shape_info[9])) * 1 * 1 * 1 * 1 * (shape_info[1] + 0));

																uint cos_sin_b = b < (shape_info[10]) ? b : 0;
																uint cos_sin_h = h < 1 ? h : 0;
																uint cos_sin_p = p;
																cos_sin_p = cos_sin_p < (shape_info[16]) ? cos_sin_p : 0;

																uint cos_sin_idx = ((1 * 0) + (64 * 0) + ((64 * (shape_info[16] + 0)) * 0) + ((64 * (shape_info[16] + 0) * 1) * 0) + ((64 * (shape_info[16] + 0) * 1 * 1 * 1 * 1) * 0) + ((64 * (shape_info[16] + 0) * 1 * 1 * 1 * 1 * 1) * 0)) + (0) * 1 + (cos_sin_p) * 64 + (cos_sin_h) * (64 * (shape_info[16] + 0) * 1 * 1 * 1 * 1) + (cos_sin_b) * (64 * (shape_info[16] + 0) * 1 * 1 * 1 * 1 * 1);
																uint cos_idx = cos_sin_idx;
																uint sin_idx = cos_sin_idx;

																uint output_idx = ((1 * 0) + (64 * 0) + ((64 * (shape_info[32] + 0)) * 0) + ((64 * (shape_info[32] + 0) * 1) * 0) + ((64 * (shape_info[32] + 0) * 1 * 1 * 1 * 1) * 0) + ((64 * (shape_info[32] + 0) * 1 * 1 * 1 * 1 * 14) * 0)) + (0) * 1 + (p) * 64 + (h) * (64 * (shape_info[32] + 0) * 1 * 1 * 1 * 1) + (b) * (64 * (shape_info[32] + 0) * 1 * 1 * 1 * 1 * 14);
																sycl::half in1 = input[input_idx + r];
																sycl::half in2 = input[input_idx + 32 + r];
																output[output_idx + r] = cos[cos_idx + r] * in1 - sin[sin_idx + r] * in2;
																output[output_idx + 32 + r] = cos[cos_idx + 32 + r] * in2 +
																							sin[sin_idx + 32 + r] * in1;
																							// out << "output[" << output_idx << "]=" << output[output_idx + r] << sycl::endl;
															} 
														}); });
		cur_event.wait();
		auto t2 = std::chrono::high_resolution_clock::now();
		if (test_performance)
		{
			auto kernel_tm = calc_kernel_tm(cur_event);
			auto tm = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
			// std::cout << "  == Infer " << i << ", time = " << tm << " micro sec." << std::endl;
			if (i >= 100)
			{
				profiling_host_tm->add(tm);
				profiling_kernel_tm->add(kernel_tm);

				if (loop_num == i + 1)
				{
					profiling_host_tm->print("--->host   tm:");
					profiling_kernel_tm->print("--->kernel tm:");
					profiling_host_tm->print_host_kernel_diff("---> host - kernel tm:", profiling_kernel_tm);
				}
			}
		}
	}
	return ret_ev;
}

static sycl::event launchSyclKernel_expand_shape(sycl::queue &q, sycl::half *buf1, sycl::half *buf2, sycl::half *buf3, sycl::half *out_buf, bool test_performance)
{
	sycl::nd_range ndr = sycl::nd_range{sycl::range{1, 14, 192}, sycl::range{1, 2, 192}};
	auto *input = buf1;
	auto *cos = buf2;
	auto *sin = buf3;
	auto *output = out_buf;

	// 0, 1, 2, 3, 4, 5, 6,  7,  8, 9,10,11,12,13,14,15,16, 17, 18,19,20,21,22,23,24,25,26, 27,28,29,30,31,32, 33
	//{1, 6, 1, 1, 1, 1, 14, 64, 0, 4, 1, 1, 1, 1, 1, 1, 6, 64, 1, 1, 1, 1, 1, 1, 6, 64, 1, 14, 1, 1, 1, 1, 6, 64}

	sycl::event ret_ev;
	size_t loop_num = test_performance ? 15 : 1;
	for (size_t i = 0; i < loop_num; i++)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		q.submit([&](sycl::handler &cgh)
				 { 
					// auto out = sycl::stream(1024, 768, cgh);
					cgh.parallel_for<class rope_sycl_exp_shape>(ndr, [=](sycl::nd_item<3> itm)
													 {
										const uint b = itm.get_global_id(0);
										const uint h = itm.get_global_id(1);
										const uint p = (uint)itm.get_global_id(2) / 32;
										const uint r = (uint)itm.get_global_id(2) % 32;

										uint input_idx = (h) * 64 + (p) * 1152 + (b) * 6912;

										uint cos_sin_b = b < 1 ? b : 0;
										uint cos_sin_h = h < 1 ? h : 0;
										uint cos_sin_p = p;
										cos_sin_p = cos_sin_p < 6 ? cos_sin_p : 0;

										uint cos_sin_idx = (cos_sin_p) * 64 + (cos_sin_h) * 384 + (cos_sin_b) * 384;
										uint cos_idx = cos_sin_idx;
										uint sin_idx = cos_sin_idx;

										uint output_idx = (p) * 64 + (h) * 384 + (b) * 5376;
										sycl::half in1 = input[input_idx + r];
										sycl::half in2 = input[input_idx + 32 + r];
										output[output_idx + r] = cos[cos_idx + r] * in1 - sin[sin_idx + r] * in2;
										output[output_idx + 32 + r] = cos[cos_idx + 32 + r] * in2 +
																	  sin[sin_idx + 32 + r] * in1;
										// out << "output[" << output_idx << "]=" << output[output_idx + r] << sycl::endl;
																	   }); })
			.wait();
		auto t2 = std::chrono::high_resolution_clock::now();
		if (test_performance)
			std::cout << "  == Infer " << i << ", time = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " micro sec." << std::endl;
	}
	return ret_ev;
}

int test_sycl_olc_interoperate_l0_backend_rope_ref()
{
	std::cout << "  OCL_KERNEL=1:  [Default] SYCL pipeline + OpenCL C kernel." << std::endl;
	std::cout << "  SYCL_KERNEL=1: Test SYCL kernel in sycl pipeline.\n"
			  << std::endl;
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
	std::string kernel_path = "../src/kernel_rope_ref/";
	std::string kernel_source = test_performance ? load_kernel(kernel_path + "SYCL_LZ_program_1_bucket_0_part_53_8491392767821923070.cl") : load_kernel(kernel_path + "SYCL_LZ_program_1_bucket_0_part_53_8491392767821923070_print_result.cl");

	// std::cout << "  kernel_source = " << kernel_source << std::endl;

	// If you want to use OpenCL backend,
	// just set env: ONEAPI_DEVICE_SELECTOR=OPENCL:GPU
	sycl::queue queue = sycl::queue(sycl::gpu_selector_v, sycl::property::queue::enable_profiling());

	std::cout << "  == Using "
			  << queue.get_device().get_info<sycl::info::device::name>()
			  << ", Backend: " << queue.get_backend()
			  << std::endl;

	sycl::event ev;
	// std::vector<std::pair<sycl::buffer<uint8_t, 1>, bool>> params;
	std::vector<std::pair<void *, size_t>> params;

	// gws=[1, 14, 192] lws=[1, 2, 192]
	// input0: data_type=i32, {1, 6, 1, 1, 1, 1, 14, 64, 0, 4, 1, 1, 1, 1, 1, 1, 6, 64, 1, 1, 1, 1, 1, 1, 6, 64, 1, 14, 1, 1, 1, 1, 6, 64}
	// input1: data_type=f16;format=bfyx;shape=[1,6,14,64];program1_network1_0_rope___module.model.layers.1.self_attn_aten__add_Add_src0.txt
	// pad_l=[0, 0, 0, 0, 0, 0, 0, 0, 0, ];
	// pad_u=[0, 0, 4, 0, 0, 0, 0, 0, 0, ];
	// dyn_pad_dims=[000000100];
	// input2: data_type=f16;format=bfyx;shape=[1,1,6,64];
	// input3: data_type=f16;format=bfyx;shape=[1,1,6,64];
	// input4: data_type=f16;format=bfyx;shape=[1,14,6,64];

	DumpData input0;
	input0.data = std::vector<float>({1, 6, 1, 1, 1, 1, 14, 64, 0, 4, 1, 1, 1, 1, 1, 1, 6, 64, 1, 1, 1, 1, 1, 1, 6, 64, 1, 14, 1, 1, 1, 1, 6, 64});
	input0.format = "bfyx";
	input0.shape = {34};

	auto input1 = load_dump_data(kernel_path + "program1_network1_0_rope___module.model.layers.0.self_attn_aten__add_Add_src0.txt");
	auto input2 = load_dump_data(kernel_path + "program1_network1_0_rope___module.model.layers.0.self_attn_aten__add_Add_src1.txt");
	auto input3 = load_dump_data(kernel_path + "program1_network1_0_rope___module.model.layers.0.self_attn_aten__add_Add_src2.txt");
	auto output_expected = load_dump_data(kernel_path + "program1_network1_0_rope___module.model.layers.0.self_attn_aten__add_Add_dst0.txt");

	auto buf_host_0 = input0.to_int(queue);
	auto buf_host_1 = input1.to_half(queue);
	auto buf_host_2 = input2.to_half(queue);
	auto buf_host_3 = input3.to_half(queue);

	// Device memory will trigger app crash.
#if 1
	auto buf_dev_0 = sycl::malloc_device(input0.data.size() * sizeof(int), queue);
	auto buf_dev_1 = sycl::malloc_device(input1.data.size() * sizeof(sycl::half), queue);
	auto buf_dev_2 = sycl::malloc_device(input2.data.size() * sizeof(sycl::half), queue);
	auto buf_dev_3 = sycl::malloc_device(input3.data.size() * sizeof(sycl::half), queue);
	auto output_host = sycl::malloc_host<sycl::half>(output_expected.data.size(), queue);
	auto output_dev = sycl::malloc_device(output_expected.data.size() * sizeof(sycl::half), queue);
#else
	// aligned_alloc_device, aligned_alloc_shared
	// auto buf_dev_0 = sycl::aligned_alloc_shared(128, input0.data.size() * sizeof(int) + 64, queue);
	// auto buf_dev_1 = sycl::aligned_alloc_shared(128, input1.data.size() * sizeof(sycl::half) + 64, queue);
	// auto buf_dev_2 = sycl::aligned_alloc_shared(128, input2.data.size() * sizeof(sycl::half) + 64, queue);
	// auto buf_dev_3 = sycl::aligned_alloc_shared(128, input3.data.size() * sizeof(sycl::half) + 64, queue);
	// auto output_host = sycl::malloc_host<sycl::half>(output_expected.data.size(), queue);
	// auto output_dev = sycl::aligned_alloc_shared(128, output_expected.data.size() * sizeof(sycl::half) + 64, queue);

	auto buf_dev_0 = sycl::malloc_shared(input0.data.size() * sizeof(int) + 64, queue);
	auto buf_dev_1 = sycl::malloc_shared(input1.data.size() * sizeof(sycl::half) + 64, queue);
	auto buf_dev_2 = sycl::malloc_shared(input2.data.size() * sizeof(sycl::half) + 64, queue);
	auto buf_dev_3 = sycl::malloc_shared(input3.data.size() * sizeof(sycl::half) + 64, queue);
	auto output_host = sycl::malloc_host<sycl::half>(output_expected.data.size(), queue);
	auto output_dev = sycl::malloc_shared(output_expected.data.size() * sizeof(sycl::half) + 64, queue);
#endif

	// Host to Device
	queue.memcpy(buf_dev_0, buf_host_0, input0.data.size() * sizeof(int)).wait();
	queue.memcpy(buf_dev_1, buf_host_1, input1.data.size() * sizeof(sycl::half)).wait();
	queue.memcpy(buf_dev_2, buf_host_2, input2.data.size() * sizeof(sycl::half)).wait();
	queue.memcpy(buf_dev_3, buf_host_3, input3.data.size() * sizeof(sycl::half)).wait();

	if (test_sycl_kernel)
	{
		std::cout << "  == run sycl kernel." << std::endl;
		auto ret_ev = launchSyclKernel(queue, (int *)buf_dev_0, (sycl::half *)buf_dev_1,
									   (sycl::half *)buf_dev_2, (sycl::half *)buf_dev_3,
									   (sycl::half *)output_dev, test_performance);
		// auto ret_ev = launchSyclKernel_expand_shape(queue, buf1, buf2, buf3, output_buf, test_performance);
		// ret_ev.wait();
	}
	else
	{
		std::cout << "  == run OpenCL C kernel." << std::endl;
		params.push_back({buf_dev_0, input0.data.size() * sizeof(int)});
		params.push_back({buf_dev_1, input1.data.size() * sizeof(sycl::half)});
		params.push_back({buf_dev_2, input2.data.size() * sizeof(sycl::half)});
		params.push_back({buf_dev_3, input3.data.size() * sizeof(sycl::half)});
		params.push_back({output_dev, output_expected.data.size() * sizeof(sycl::half)});

		auto ret_ev = launchOpenCLKernelOnline(queue, kernel_source, "rope_ref_11982042700243959200_0_0__sa", params, ev, test_performance);
		ret_ev.wait();
	}

	// Device to Host
	queue.memcpy(output_host, output_dev, output_expected.data.size() * sizeof(sycl::half)).wait();

	// Compare result.
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
	sycl::free(buf_dev_0, queue);
	sycl::free(buf_host_1, queue);
	sycl::free(buf_dev_1, queue);
	sycl::free(buf_host_2, queue);
	sycl::free(buf_dev_2, queue);
	sycl::free(buf_host_3, queue);
	sycl::free(buf_dev_3, queue);
	sycl::free(output_host, queue);
	sycl::free(output_dev, queue);

	return 0;
}
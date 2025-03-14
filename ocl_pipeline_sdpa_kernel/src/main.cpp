// Reference:

#include <stdio.h>
#include <iostream>
#include <CL/opencl.hpp>
#include <stddef.h>
#include <stdint.h>
#include <fstream>
#include <iomanip>

#include "private.hpp"

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
	std::cout << "  == Fail: Can't open: " << kernel_fn << std::endl;
	exit(0);
	return std::string();
}

int test_build_asm(cl::Device& default_device, bool test_performance = true) {
	// const char *kernel_code = R""""(
    //     // Kernel name: igc_check
    //     kernel void igc_check() {
    //         __asm__ volatile(
    //                 ".decl AA0 v_type=G type=ud num_elts=1\n"
    //                 ".decl AA1 v_type=G type=ud num_elts=1\n"
    //                 ".implicit_PSEUDO_INPUT AA0 offset=256 size=4\n"
    //                 ".implicit_PSEUDO_INPUT AA1 offset=256 size=4\n"
    //                 "mov (M1_NM,1) AA0(0,0)<1> AA1(0,0)<0;1,0>\n"
    //         );
    //     }
    //     )"""";

	std::string kernel_code = load_kernel("../../sycl_pipeline_sdpa_kernel/src/kernel_sdpa_micro/sdpa_micro_prefill_8660372428234100028_0_0__sa.cl");

	std::cout << "== Create context" << std::endl;
	cl::Context context({default_device});

	std::cout << "== Create Sources" << std::endl;
	cl::Program::Sources sources;

	std::cout << "== Put kernel string to source." << std::endl;
	sources.push_back({kernel_code.c_str(), strlen(kernel_code.c_str())});

	std::cout << "== Construct program with source and context." << std::endl;
	cl::Program program(context, sources);

	std::string options = "-cl-mad-enable -cl-std=CL2.0 -cl-intel-256-GRF-per-thread -Dcl_intel_dot_accumulate -Dcl_intel_global_float_atomic -Dcl_intel_subgroup_matrix_multiply_accumulate -Dcl_intel_subgroup_split_matrix_multiply_accumulate";
	if (program.build({default_device}, options.c_str()) != CL_SUCCESS)
	{
		std::cout << " Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << "\n";
		exit(1);
	}

	// Construct kernel 1
	cl::vector<cl::Kernel> kernels;
	program.createKernels(&kernels);
	for (auto func_name : kernels)
	{
		auto kernel_name = func_name.getInfo<CL_KERNEL_FUNCTION_NAME>();
		std::cout << "  == Get kernel function name from  = " << kernel_name << std::endl;
	}

	std::cout << "== Create command queue" << std::endl;
	// create queue to which we will push commands for the device.
	cl::CommandQueue queue(context, default_device);

	cl::Kernel kernel_asm = cl::Kernel(program, "igc_check");

	auto kernel_name = kernel_asm.getInfo<CL_KERNEL_FUNCTION_NAME>();
	std::cout << "== Test get kernel name from cl::Kernel, kernel_name = " << kernel_name << std::endl;
	size_t loop_num = test_performance ? 150 : 1;
	for (size_t i = 0; i < loop_num; i++)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		queue.enqueueNDRangeKernel(kernel_asm, cl::NullRange, cl::NDRange(28, 64, 256), cl::NDRange(1, 1, 256));
		queue.finish();
		auto t2 = std::chrono::high_resolution_clock::now();
		if (test_performance)
			std::cout << "== Infer " << i << ", time = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " micro sec." << std::endl;
	}

	std::cout << "== Read result." << std::endl;

	return 0;
}

int main(int argc, char **argv)
{
	std::cout << "== Debug MACRO tip:" << std::endl;
	std::cout << "  Default:        test accuracy only." << std::endl;
	std::cout << "  PERFORMANCE=1:  Test performance only." << std::endl;
	std::string options = "-cl-mad-enable -cl-std=CL2.0";
	std::cout << "== Usage:" << std::endl;
	std::cout << "  ./ocl_pipeline_sdpa_kernel [result_fn]" << std::endl;

	bool test_performance = false;
	if (std::getenv("PERFORMANCE"))
	{
		test_performance = std::getenv("PERFORMANCE") == std::string("1");
	}

	std::cout << "== ocl_pipeline_sdpa_kernel." << std::endl;
	// get all platforms (drivers)
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);
	size_t selected_platform = 0;
	for (size_t i = 0; i < all_platforms.size(); i++)
	{
		std::string platname = all_platforms[i].getInfo<CL_PLATFORM_NAME>();
		if (platname.find("Graphics") != std::string::npos)
		{
			selected_platform = i;
			break;
		}
	}
	if (all_platforms.size() == 0)
	{
		std::cout << " No platforms found. Check OpenCL installation!\n";
		exit(1);
	}
	cl::Platform default_platform = all_platforms[selected_platform];
	std::cout << "Using platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";

	// get default device of the default platform
	std::vector<cl::Device> all_devices;
	default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
	if (all_devices.size() == 0)
	{
		std::cout << " No devices found. Check OpenCL installation!\n";
		exit(1);
	}
	cl::Device default_device = all_devices[0];
	std::cout << "Using device: " << default_device.getInfo<CL_DEVICE_NAME>() << "\n";

	return test_build_asm(default_device);

	std::cout << "== Create context" << std::endl;
	cl::Context context({default_device});

	std::cout << "== Create Sources" << std::endl;
	cl::Program::Sources sources;

	// kernel
	auto kernel_fn = "../../sycl_pipeline_sdpa_kernel/src/kernel_sdpa/SYCL_LZ_program_1_bucket_0_part_67_8198298808319190670.cl";
	std::cout << "== load kernel: " << kernel_fn << std::endl;
	std::string kernel_code = load_kernel(kernel_fn);

	std::cout << "== Put kernel string to source." << std::endl;
	sources.push_back({kernel_code.c_str(), kernel_code.length()});

	std::cout << "== Construct program with source and context." << std::endl;
	cl::Program program(context, sources);

	if (program.build({default_device}, options.c_str()) != CL_SUCCESS)
	{
		std::cout << " Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << "\n";
		exit(1);
	}

	// Construct kernel 1
	cl::vector<cl::Kernel> kernels;
	program.createKernels(&kernels);
	for (auto func_name : kernels)
	{
		auto kernel_name = func_name.getInfo<CL_KERNEL_FUNCTION_NAME>();
		std::cout << "  == Get kernel function name from  = " << kernel_name << std::endl;
	}

	DumpData input_0;
	input_0.data = std::vector<float>({1, 28, 1, 1, 1, 1, 1024, 128,
									   1, 4, 1, 1, 1, 1, 1024, 128,
									   0, 37, 1, 4, 1, 1, 1, 1, 1024, 128,
									   0, 136, 1, 1, 1, 1, 1, 1, 1024, 1024,
									   1, 1, 1, 1, 1, 1, 1024, 1,
									   1, 28, 1, 1, 1, 1, 1024, 128});
	input_0.format = "bfyx";
	input_0.shape = {52};

	// Dump data is from OpenVINO.
	std::string root_path = "../../sycl_pipeline_sdpa_kernel/src/kernel_sdpa/";
	std::string param_name = "program1_network1_0_indirectsdpa___module.model.layers.0.self_attn_aten__scaled_dot_product_attention_ScaledDotProductAttention_";

	auto input_1 = load_dump_data(root_path + param_name + "src0.txt");
	auto input_2 = load_dump_data(root_path + param_name + "src1.txt");
	auto input_3 = load_dump_data(root_path + param_name + "src2.txt");
	auto input_4 = load_dump_data(root_path + param_name + "src3.txt");
	auto input_5 = load_dump_data(root_path + param_name + "src4.txt");
	auto output_expected = load_dump_data(root_path + param_name + "dst0.txt");

	// create buffers on the device
	cl::Buffer inp_buf_0(context, CL_MEM_READ_ONLY, sizeof(int) * input_0.data.size());
	cl::Buffer inp_buf_1(context, CL_MEM_READ_ONLY, sizeof(cl_half16) * input_1.data.size());
	cl::Buffer inp_buf_2(context, CL_MEM_READ_ONLY, sizeof(cl_half16) * input_2.data.size());
	cl::Buffer inp_buf_3(context, CL_MEM_READ_ONLY, sizeof(cl_half16) * input_3.data.size());
	cl::Buffer inp_buf_4(context, CL_MEM_READ_ONLY, sizeof(cl_half16) * input_4.data.size());
	cl::Buffer inp_buf_5(context, CL_MEM_READ_WRITE, sizeof(cl_half16) * input_5.data.size());
	cl::Buffer inp_buf_6(context, CL_MEM_READ_WRITE, sizeof(float));
	cl::Buffer inp_buf_7(context, CL_MEM_READ_WRITE, sizeof(float));
	cl::Buffer inp_buf_8(context, CL_MEM_READ_WRITE, sizeof(cl_half16));

	std::cout << "== Create command queue" << std::endl;
	// create queue to which we will push commands for the device.
	cl::CommandQueue queue(context, default_device);

	// write arrays A and B to the device
	auto buf_0 = input_0.to_int();
	auto buf_1 = input_1.to_half();
	auto buf_2 = input_2.to_half();
	auto buf_3 = input_3.to_half();
	
	std::cout << "== enqueueWriteBuffer" << std::endl;
	queue.enqueueWriteBuffer(inp_buf_0, CL_TRUE, 0, sizeof(int) * input_0.data.size(), buf_0);
	queue.enqueueWriteBuffer(inp_buf_1, CL_TRUE, 0, sizeof(ushort) * input_1.data.size(), buf_1);
	queue.enqueueWriteBuffer(inp_buf_2, CL_TRUE, 0, sizeof(ushort) * input_2.data.size(), buf_2);
	queue.enqueueWriteBuffer(inp_buf_3, CL_TRUE, 0, sizeof(ushort) * input_3.data.size(), buf_3);
	queue.enqueueWriteBuffer(inp_buf_4, CL_TRUE, 0, sizeof(ushort) * input_4.data.size(), buf_3);
	queue.enqueueWriteBuffer(inp_buf_5, CL_TRUE, 0, sizeof(ushort) * input_5.data.size(), buf_3);

	float f32_zero = 0;
	ushort f16_zero = 0;
	queue.enqueueWriteBuffer(inp_buf_6, CL_TRUE, 0, sizeof(float), &f32_zero);
	queue.enqueueWriteBuffer(inp_buf_7, CL_TRUE, 0, sizeof(float), &f32_zero);
	queue.enqueueWriteBuffer(inp_buf_8, CL_TRUE, 0, sizeof(ushort), &f16_zero);
	
	std::cout << "== Create Kernel with program and run." << std::endl;
	// alternative way to run the kernel
	cl::Kernel kernel_sdpa = cl::Kernel(program, "sdpa_opt_multi_tokens_8660372428234100028_0_0__sa");    // Construct kernel 2
	kernel_sdpa.setArg(0, inp_buf_0);
	kernel_sdpa.setArg(1, inp_buf_1);
	kernel_sdpa.setArg(2, inp_buf_2);
	kernel_sdpa.setArg(3, inp_buf_3);
	kernel_sdpa.setArg(4, inp_buf_4);
	kernel_sdpa.setArg(5, inp_buf_5);
	kernel_sdpa.setArg(6, inp_buf_6);
	kernel_sdpa.setArg(7, inp_buf_7);
	kernel_sdpa.setArg(8, inp_buf_8);

	auto kernel_name = kernel_sdpa.getInfo<CL_KERNEL_FUNCTION_NAME>();
	std::cout << "== Test get kernel name from cl::Kernel, kernel_name = " << kernel_name << std::endl;
	size_t loop_num = test_performance ? 150 : 1;
	for (size_t i = 0; i < loop_num; i++)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		queue.enqueueNDRangeKernel(kernel_sdpa, cl::NullRange, cl::NDRange(28, 64, 256), cl::NDRange(1, 1, 256));
		queue.finish();
		auto t2 = std::chrono::high_resolution_clock::now();
		if (test_performance)
			std::cout << "== Infer " << i << ", time = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " micro sec." << std::endl;
	}

	std::cout << "== Read result." << std::endl;
	// read result C from the device to array C
	// ushort* output = (ushort*)malloc(sizeof(ushort) * output_expected.data.size());
	// queue.enqueueReadBuffer(out_buf_0, CL_TRUE, 0, sizeof(ushort) * output_expected.data.size(), output);

	// if (!test_performance)
	// {
	// 	std::cout << "== Start compare result with expected.\n";
	// 	bool is_expected = true;
	// 	for (int i = 0; i < output_expected.data.size(); i++)
	// 	{
	// 		if (std::to_string(output_expected.data[i]) != std::to_string(half_to_float(output[i])))
	// 		{
	// 			std::cout << "Index: " << i << " Result " << half_to_float(output[i]) << "!=" << " Expected " << output_expected.data[i] << std::endl;
	// 			is_expected = false;
	// 		}
	// 	}
	// }

	// Dump result:
	// std::string dump_fn = argc == 2 ? argv[1] : "out_opencl_sdpa_kernel.log";
	// std::cout << "== Start dump result: " << dump_fn << std::endl;
	// std::filebuf fb;
	// fb.open(dump_fn, std::ios::out);
	// std::ostream os(&fb);
	// os.precision(6);
	// for (int i = 0; i < output_expected.data.size(); i++)
	// {
	// 	os << std::fixed<< half_to_float(output[i]) << std::endl;
	// }
	// fb.close();

	// std::cout << "== Result is_expected = " << is_expected << std::endl;
	std::cout << "== Done." << std::endl;
	return 0;
}
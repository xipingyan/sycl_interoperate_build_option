#include "private.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>

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

// Reference:
// https://github.com/intel/compute-runtime/blob/master/level_zero/core/test/black_box_tests/zello_world_jitc_ocloc.cpp
// https://github.com/intel/compute-runtime/blob/master/level_zero/core/test/black_box_tests/common/zello_compile.cpp

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

#include <sstream> 
#include "ze_api_wrap.hpp"

std::string get_device_id(ze_device_handle_t device)
{
	ze_device_properties_t deviceProperties = {ZE_STRUCTURE_TYPE_DEVICE_PROPERTIES};
	SUCCESS_OR_TERMINATE(zeDeviceGetProperties(device, &deviceProperties));
	// LevelZeroBlackBoxTests::printDeviceProperties(deviceProperties);

	std::stringstream ss;
	ss.setf(std::ios::hex, std::ios::basefield);
	ss << "0x" << deviceProperties.deviceId;

	std::cout << "  == Get device id: " << ss.str() << std::endl;
	return ss.str();
}

std::vector<uint8_t> compile_ocl_to_bin(sycl::queue &q, const std::string& source, std::string func_name) {
	auto zeDevice = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(
		q.get_device());
	auto zeContext = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(
		q.get_context());

	auto device = get_device_id(zeDevice);
	std::string revisionId = "0x2f";

	std::string mainFileName = "src.cl";
	FILE *pf = fopen(mainFileName.c_str(), "wb");
	fwrite(source.c_str(), source.length(), 1, pf);
	fclose(pf);
	std::string output = "output.bin";

	const char *argv[] = {"ocloc", "-v", "-device", device.c_str(), "-revision_id", revisionId.c_str(), "-file", mainFileName.c_str(), "-o", "output.bin", "", "", "", "", "", ""};
	int param_num = sizeof(argv) / sizeof(char *);
	std::string cmd;
	for (int i = 0; i < param_num; i++) {
		cmd += argv[i] + std::string(" ");
	}
	system(cmd.c_str());

	auto binary = load_kernel(output);
	std::vector<uint8_t> ret;
	ret.assign(binary.begin(), binary.begin() + binary.length());
	return ret;
}

using LoadedDeviceQueueProperties = std::map<ze_device_handle_t, std::vector<ze_command_queue_group_properties_t>>;
static LoadedDeviceQueueProperties deviceQueueProperties;

std::vector<ze_command_queue_group_properties_t> &getDeviceQueueProperties(ze_device_handle_t device) {
    auto &queueProperties = deviceQueueProperties[device];
    if (queueProperties.size() == 0) {
        uint32_t numQueueGroups = 0;
        SUCCESS_OR_TERMINATE(zeDeviceGetCommandQueueGroupProperties(device, &numQueueGroups, nullptr));
        if (numQueueGroups == 0) {
            std::cerr << "No queue groups found!\n";
            std::terminate();
        }
        queueProperties.resize(numQueueGroups);
        SUCCESS_OR_TERMINATE(zeDeviceGetCommandQueueGroupProperties(device, &numQueueGroups,
                                                                    queueProperties.data()));
    }
    return queueProperties;
}

uint32_t getCommandQueueOrdinal(ze_device_handle_t &device, bool useCooperativeFlag)
{
	std::vector<ze_command_queue_group_properties_t> &queueProperties = getDeviceQueueProperties(device);

	ze_command_queue_group_property_flags_t computeFlags = ZE_COMMAND_QUEUE_GROUP_PROPERTY_FLAG_COMPUTE;
	if (useCooperativeFlag)
	{
		computeFlags |= ZE_COMMAND_QUEUE_GROUP_PROPERTY_FLAG_COOPERATIVE_KERNELS;
	}

	uint32_t computeQueueGroupOrdinal = std::numeric_limits<uint32_t>::max();
	for (uint32_t i = 0; i < queueProperties.size(); i++)
	{
		if (queueProperties[i].flags & computeFlags)
		{
			computeQueueGroupOrdinal = i;
			break;
		}
	}
	return computeQueueGroupOrdinal;
}

sycl::event launchOpenCLKernel_OCLOC(sycl::queue &q, std::string source,
									 std::string func_name, std::vector<std::pair<void *, size_t>> &params,
									 sycl::event &dep_event, bool test_performance)
{
	auto device = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(
		q.get_device());
	auto context = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(
		q.get_context());
	auto qqq = sycl::get_native<sycl::backend::ext_oneapi_level_zero>(q);

	std::cout << "  == Start to kernel_bundle opencl source" << std::endl;
	auto bin = compile_ocl_to_bin(q, source, func_name);

	std::string buildFlags = "-cl-std=CL3.0";

	// Create module
	ze_module_handle_t module = nullptr;
	ze_module_desc_t moduleDesc = {ZE_STRUCTURE_TYPE_MODULE_DESC};
    moduleDesc.format = ZE_MODULE_FORMAT_NATIVE;
	moduleDesc.pInputModule = bin.data();
	moduleDesc.inputSize = bin.size();
	moduleDesc.pBuildFlags = buildFlags.c_str();
	SUCCESS_OR_TERMINATE(zeModuleCreate(context, device, &moduleDesc, &module, nullptr));

	// create kernel
	ze_kernel_handle_t kernel = nullptr;
    ze_kernel_desc_t kernelDesc = {ZE_STRUCTURE_TYPE_KERNEL_DESC};
    kernelDesc.pKernelName = func_name.c_str();
    SUCCESS_OR_TERMINATE(zeKernelCreate(module, &kernelDesc, &kernel));


    ze_result_t result;
    ze_command_queue_desc_t cmdQueueDesc = {ZE_STRUCTURE_TYPE_COMMAND_QUEUE_DESC};
    cmdQueueDesc.ordinal = getCommandQueueOrdinal(device, false);
    cmdQueueDesc.index = 0;

	// if (isImmediate) {
	cmdQueueDesc.mode = ZE_COMMAND_QUEUE_MODE_SYNCHRONOUS;
	ze_command_list_handle_t cmdList;
	SUCCESS_OR_TERMINATE(zeCommandListCreateImmediate(context, device, &cmdQueueDesc, &cmdList));

	uint32_t poolSize = 1;
	ze_event_pool_handle_t eventPool;
 	ze_event_pool_desc_t eventPoolDesc{ZE_STRUCTURE_TYPE_EVENT_POOL_DESC};
    eventPoolDesc.count = poolSize;
    eventPoolDesc.flags = ZE_EVENT_POOL_FLAG_HOST_VISIBLE;
    SUCCESS_OR_TERMINATE(zeEventPoolCreate(context, &eventPoolDesc, 1, &device, &eventPool));

    ze_event_handle_t event;
    ze_event_desc_t eventDesc = {ZE_STRUCTURE_TYPE_EVENT_DESC};
    for (uint32_t i = 0; i < poolSize; i++) {
        // if (counterEvents) {
        //     SUCCESS_OR_TERMINATE(zexCounterBasedEventCreate2Func(context, device, counterBasedDesc, events + i));
        // } else {
            eventDesc.index = i;
            eventDesc.signal = ZE_EVENT_SCOPE_FLAG_HOST;
            eventDesc.wait = ZE_EVENT_SCOPE_FLAG_HOST;
            SUCCESS_OR_TERMINATE(zeEventCreate(eventPool, &eventDesc, &event));
        // }
    }
	std::cout << "--------------1" << std::endl;
	for (size_t i = 0; i < params.size(); i++)
	{
		try
		{
			// sycl malloc's usm buffer, seems not work here.
			// crash.
			zeKernelSetArgumentValue(kernel, i, sizeof(void *), params[i].first);
		}
		catch (const std::exception &e)
		{
			std::cerr << e.what() << '\n';
		}
	}
std::cout << "--------------2" << std::endl;
	uint32_t groupSizeX = 1;
	uint32_t groupSizeY = 2u;
	uint32_t groupSizeZ = 192u;
	SUCCESS_OR_TERMINATE(zeKernelSetGroupSize(kernel, groupSizeX, groupSizeY, groupSizeZ));

	ze_group_count_t dispatchTraits;
	dispatchTraits.groupCountX = 1;
	dispatchTraits.groupCountY = 14u;
	dispatchTraits.groupCountZ = 192u;

	ze_command_queue_handle_t cmdQueue;
	ze_command_queue_desc_t descriptor = {};
	descriptor.stype = ZE_STRUCTURE_TYPE_COMMAND_QUEUE_DESC;

	descriptor.pNext = nullptr;
	descriptor.flags = 0;
	descriptor.mode = ZE_COMMAND_QUEUE_MODE_DEFAULT;
	descriptor.priority = ZE_COMMAND_QUEUE_PRIORITY_NORMAL;

	descriptor.ordinal = getCommandQueueOrdinal(device, false);
	descriptor.index = 0;
	SUCCESS_OR_TERMINATE(zeCommandQueueCreate(context, device, &descriptor, &cmdQueue));
std::cout << "--------------3" << std::endl;
	for (auto i = 0; i < 150; i++)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		auto hTSEvent = create_event_timestamp(eventPool);
		// std::cout << "Create even timestamp: hTSEvent = " << hTSEvent << std::endl;

		// Append a signal of a timestamp event into the command list after the kernel executes
		auto r = zeCommandListAppendLaunchKernel(cmdList, kernel, &dispatchTraits, hTSEvent, 0, nullptr);
		CHECK_RET(r)

		SUCCESS_OR_TERMINATE(zeCommandListClose(cmdList));

		// Execute the command list with the signal
		// std::cout << "Command queue start to execute command list." << std::endl;
		r = zeCommandQueueExecuteCommandLists(cmdQueue, 1, &cmdList, nullptr);
		CHECK_RET(r)

		zeCommandListReset(cmdList);

		// r = zeEventHostSynchronize(hEvent, 0);
		// CHECK_RET(r)

		auto t2 = std::chrono::high_resolution_clock::now();
		std::cout << "i = " << i << ", tm = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " micro sec.\n";
	}

	return sycl::event();
}

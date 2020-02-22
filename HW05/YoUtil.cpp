#include "YoUtil.hpp"
#include <iostream>			// for I/O
#include <fstream>			// for file I/O

using namespace std;

namespace YoUtil {

bool GPU::verbose = true;
string GPU::compilationOptions = "-cl-std=CL1.2 -w -cl-kernel-arg-info";

GPU::GPU(cl_device_type deviceType, bool enableProfiling, const string& vendorNameFilter, const string& vendorVersionFilter) {
	vector< cl::Platform > platforms;
	cl::Platform::get(&platforms);
	
	for (cl::Platform &platform : platforms) {
		string vendor = platform.getInfo<CL_PLATFORM_VENDOR>();
		string vendorVersion = platform.getInfo<CL_PLATFORM_VERSION>();
		if( vendorNameFilter.size()>0 && vendor.find( vendorNameFilter ) == vendor.npos ) continue;
		if (vendorVersionFilter.size()>0 && vendorVersion.find(vendorVersionFilter) == vendor.npos) continue;

		platform.getDevices(deviceType, &devices);
		if (devices.size() == 0) continue;
		// create a context that can manage all devices given in devices
		context = new cl::Context(devices);

		// create command queues in the context
		devices = context->getInfo<CL_CONTEXT_DEVICES>();
		cl_command_queue_properties props = 0;
		if( enableProfiling ) props = CL_QUEUE_PROFILING_ENABLE;
		for (auto& device : devices ) {
			cmdQueues.push_back(cl::CommandQueue(*context, device, props));
		}

		if (verbose) {
			cout << "\n=== GPU constructor ===";
			cout << "\n\tPlatform by : " << vendor;
			cout << "\n\tThere is/are " << context->getInfo<CL_CONTEXT_NUM_DEVICES>() << " device(s) in the defined context.";
			cout << "\n\t# of Command queues: " << cmdQueues.size() << "\n";
		}
		return ;
	}

	if (verbose) cerr << "\n\t***Cannot find any platform for given filters or device type: " << deviceType;
	context = nullptr;
}

GPU::~GPU() {
	delete context;
	for(auto& it: programs) {
		delete it.second;
	}
}

cl::Context& GPU::operator() () {
	return *context;
}

size_t GPU::getNoDevices() const {
	return devices.size();
}

cl::CommandQueue& GPU::getCommandQueue(size_t which) {
	return cmdQueues[which];
}

bool GPU::addProgramByFile(const char *programName, const char *filename) {
	string code = readSourceCode( filename );
	auto prg = compile( *context, code, compilationOptions.c_str() );
	programs[ programName ] = prg;

	return true;
}

cl::Program* GPU::getProgram(const char *kernelName) {
	auto it = programs.find( kernelName );
	if( it == programs.end() ) {
		cerr << "\nError, cannot find specified program: " << kernelName << endl;
		return nullptr;
	}
	else {
		return it->second;
	}
}

string GPU::readSourceCode(const char *filename) {
	ifstream inp(filename);
	if (!inp) {
		cerr << "\nError opening file: " << filename << endl;
		return "";
	}

	string kernel((istreambuf_iterator<char>(inp)), istreambuf_iterator<char>());

	inp.close();
	return kernel;
}

cl::Program* GPU::compile(cl::Context &context, const string &source, const char *options) {
	cl::Program *prog = new cl::Program(context, source);

	try {
		// prog->build("-cl-std=CL1.2 -w -cl-kernel-arg-info");
		prog->build(compilationOptions.c_str());
	}
	catch (cl::Error &e) {
		cerr << "\nFile: " << __FILE__ << ", line: " << __LINE__ << e.what();
		cerr << "\nError no: " << e.err() << endl;
		for (auto& device : context.getInfo<CL_CONTEXT_DEVICES>()) {
			cout << "\n=== " << device.getInfo<CL_DEVICE_NAME>() << " ===";
			cout << "\nBuild log: " << prog->getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
			cout << "\nBuild options used:" << prog->getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device);
		}
		return nullptr;
	}

	if( verbose ) {
		// See Table 5.13 of OpenCL 1.2 specification for information that can be queried to program objects
		cout << "\n\t# devices associated with the program: " << prog->getInfo<CL_PROGRAM_NUM_DEVICES>();
		cout << "\n\t# Kernels defined: " << prog->getInfo<CL_PROGRAM_NUM_KERNELS>();
		cout << "\n\tProgram kernel names: " << prog->getInfo<CL_PROGRAM_KERNEL_NAMES>();
		cout << "\n\tProg sizes: ";  for (auto s : prog->getInfo<CL_PROGRAM_BINARY_SIZES>()) cout << s << ";";
	}

	return prog;
}

//
// List all OpenCL platforms found on this host
//
void GPU::showPlatforms(ostream& outp) {
	try {
		std::vector< cl::Platform > platforms;
		cl::Platform::get(&platforms);
		outp << "\nnumPlatforms: " << platforms.size();

		// go through each platform
		for (auto &platform : platforms) {
			// See table 4.1 of OpenCL specification 1.2 for more query keys
			outp << "\n=== Platform Info ===";
			outp << "\nPlatform name: " << platform.getInfo<CL_PLATFORM_NAME>();
			outp << "\nPlatform vendor: " << platform.getInfo<CL_PLATFORM_VENDOR>();
			outp << "\nPlatform version: " << platform.getInfo<CL_PLATFORM_VERSION>();
			outp << "\nPlatform profile: " << platform.getInfo<CL_PLATFORM_PROFILE>();
			outp << "\nPlatform extensions: ";

			string exts = platform.getInfo<CL_PLATFORM_EXTENSIONS>();
			size_t pos = 0;
			while ((pos = exts.find(' ')) != string::npos) {
				outp << "\n\t\t" << exts.substr(0, pos);
				exts = exts.substr(pos + 1);
			}
		}
	}
	catch (cl::Error &e) {
		cerr << "\nMain: " << e.what();
		cerr << "\nError no: " << e.err() << endl;
	}
}

//
// Show information about found devices 
//
void GPU::showDevices(std::ostream& cout) {
	for (auto& device : devices) {
		// See table 4.3 of OpenCL speficiation 1.2 for device query keys
		cout << "\n=== Device Info ===";
		cout << "\n\tDevice name: " << device.getInfo<CL_DEVICE_NAME>();
		cout << "\n\tDevice vendor: " << device.getInfo<CL_DEVICE_VENDOR>();
		cout << "\n\tDevice profile: " << device.getInfo<CL_DEVICE_PROFILE>();
		cout << "\n\tDevice version: " << device.getInfo<CL_DEVICE_VERSION>();
		cout << "\n\tDriver version: " << device.getInfo<CL_DRIVER_VERSION>();
		cout << "\n\tDevice OpenCL_C version: " << device.getInfo<CL_DEVICE_OPENCL_C_VERSION>();

		cout << "\n\tDevice global memory size (MB): " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() / 1024 / 1024;
		cout << "\n\tDevice max. allocable memory size (MB): " << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() / 1024 / 1024;
		cout << "\n\tDevice address bits: " << device.getInfo<CL_DEVICE_ADDRESS_BITS>();
		cout << "\n\tDevice mem base address align: " << device.getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>();

		cout << "\n\tMax. clock: " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>();
		cout << "\n\tMax. compute units: " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();

		cout << "\n\tMax. Work Group size: " << device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
		cout << "\n\tMax Work Item Dimensions: " << device.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>();

		auto sizes = device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
		cout << "\n\tMax. Work Item Size: "; for (auto i : sizes) cout << i << " ";


		cout << "\n\n\tDevice extensions: ";
		// show extensions
		string exts = device.getInfo<CL_DEVICE_EXTENSIONS>();
		size_t pos = 0;
		while ((pos = exts.find(' ')) != string::npos) {
			cout << "\n\t\t" << exts.substr(0, pos);
			exts = exts.substr(pos + 1);
		}
		cout << "\n\t\t" << exts;
	}
}



} // End namespace YoUtil

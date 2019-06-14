#ifndef YO_OPENCL_UTIL_HEADER
#define YO_OPENCL_UTIL_HEADER

#define CL_HPP_ENABLE_EXCEPTIONS
#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#include <CL/cl.hpp>

#include <string>			// for C++ string
#include <vector>			// for vector container
#include <chrono>			// for Timing
#include <map>				// for storing kernel programs

namespace YoUtil {
	class GPU {
		cl::Context *context;
		std::vector<cl::CommandQueue> cmdQueues;
		std::map<std::string, cl::Program*> programs;
		std::vector<cl::Device> devices;

		// utility function to load OpenCL source code and store it into a string
		std::string readSourceCode(const char *filename);

		// utility function compile the CL source code into Programs
		cl::Program *compile(cl::Context &, const std::string &source, const char *option = "-cl-std=CL1.2 -w -cl-kernel-arg-info");
	public:
		GPU(cl_device_type type=CL_DEVICE_TYPE_GPU, bool enableProfiling=false, const std::string& vendorNameFilter="", const std::string& vendorVersionFilter="");
		~GPU();

		// get OpenCl context (needed for many OpenCL functions)
		cl::Context& operator() ();

		// show selected device information
		void showDevices(std::ostream& outp);

		// get the number of devices available to work
		size_t getNoDevices() const; 

		// get the command queue for the device specified by which
		cl::CommandQueue& getCommandQueue(size_t which); 

		// load source code of kernel files from given filename, and kernel name is simply an identifier
		bool addProgramByFile(const char* kernelName, const char *filename);

		// get the compiled program specified by kernelName
		cl::Program* getProgram(const char *kernelName);

		// List OpenCL platforms found on the current host
		static void showPlatforms(std::ostream& outp);

		static bool verbose;
		static std::string compilationOptions;
	};

}	// End of YoUtil namespace

#endif

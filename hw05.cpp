#include "hw05.hpp"
#include <fstream>
#include <iostream>
#include <string>
using namespace std;
using namespace YoUtil;

// Prepare the needed GPU stuff (context, device, commandQueue)
void prepareGPU(mygpu& gpu) {
	//cout << "\nDo not forget to implement prepareGPU! Line " << __LINE__ << "@ " << __FILE__; 
    try {
        GPU::verbose = false;
        gpu.gpu = new GPU(CL_DEVICE_TYPE_GPU, false, "");
        gpu.gpu->addProgramByFile("hw05", "hw05.cl");
        gpu.cmdQueue = gpu.gpu->getCommandQueue(0);
        gpu.device = gpu.cmdQueue.getInfo<CL_QUEUE_DEVICE>();
        gpu.context = gpu.cmdQueue.getInfo<CL_QUEUE_CONTEXT>();
        gpu.program = gpu.gpu->getProgram("hw05");
    } catch (cl::Error &e) {
        cerr << "HW05: " << e.what() << endl;
        cerr << "Error no: " << e.err() << endl;
    }
}

void releaseGPU(mygpu& gpu) {
	//cout << "\nDo not forget to implement releaseGPU! Line " << __LINE__ << "@ " << __FILE__; 
    delete gpu.gpu;
}

// Read the original point cloud data and store them into GPU buffer
bool readPointCloud(const char *filename, mydata& data, mygpu& gpu) {
	fstream inp(filename);
	if (!inp) {
		cerr << "\nError opening point cloud data file: " << filename;
		return false;
	}
	size_t nData;

	// your likely want to store the number of points somewhere in your struct
	inp >> nData;
    data.size = nData;

	// some pts file has a 'D' character after the number of points...:(
	char pp = inp.peek(); 
	if (pp == 'D' || pp == 'd') inp >> pp; 

	// you should now allocate enough memory to store all the data...
	#ifdef aos
		// prepare AoS data struct
		//cout << "\nAllocate for AoS data struct"; 
        data.point = new cl_double [data.size * 3];
	#endif
	#ifdef soa
		// prepare SoA data struct
		//cout << "\nAllocate for SoA data struct"; 
        data.point = new cl_double [data.size * 3];
	#endif

	for (size_t i = 0; i < nData; ++i) {
		double x, y, z;
		string rest;		// some pts files have only intensity, and some have RGB values, so we treat them all as a long string...
		inp >> x >> y >> z;
		std::getline(inp, rest); 
		// if(i==0) cout << "\n" << x << ", " << y << ", " << z << ": " << rest; 
		// now you have coordinates in x, y, z, and additional information in the string rest.  You need to store them into your data struct...
		// we are going to discard the additional information in this assignment.
        #ifdef aos
            data.point[i * 3    ] = x;
            data.point[i * 3 + 1] = y;
            data.point[i * 3 + 2] = z;
        #endif
        #ifdef soa
            data.point[i            ] = x;
            data.point[nData + i    ] = y;
            data.point[nData * 2 + i] = z;
        #endif
	}
	inp.close();

	return true;
}

// This function moves all points in the point cloud so that all points that was in the AABB of [xmin, ymin, zmin] - [xmax, ymax, zmax] 
// will be in the region of [-Lx, -Ly, -Lz] - [Lx, Ly, Lz], where Lx = 0.5*(xmax+xmin), Ly=0.5*(ymax+ymin), Lz=0.5*(zmax+zmin)
// Also, this function returns sqrt( (0.5*(xmax-xmin))**2 + (0.5*(ymax-ymin))**2 + (0.5*(zmax-zmin))**2) ).  This is the maximum possible rho.
// Do this with OpenCL GPU code
double centerPointCloudToOrigin(mydata &data, mygpu& gpu) {
	//cout << "\nDo not forget to implement centerPointCloudToOrigin! Line " << __LINE__ << "@ " << __FILE__; 
    cl_mem_flags flag = CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR;

    #ifdef aos
        cl_uint bytes = sizeof(data.point[0]) * data.size * 3;
        gpu.point_buf = cl::Buffer(gpu.context, flag, bytes, data.point);

        cl_double max[3] = { data.point[0], data.point[1], data.point[2] };
        cl_double min[3] = { data.point[0], data.point[1], data.point[2] };
        cl::Buffer max_buf(gpu.context, flag, sizeof(max[0]) * 3, max);
        cl::Buffer min_buf(gpu.context, flag, sizeof(min[0]) * 3, min);
    #endif
    #ifdef soa
        cl_uint bytes = sizeof(data.point[0]) * data.size * 3;
        gpu.point_buf = cl::Buffer(gpu.context, flag, bytes, data.point);

        cl_double max[3] = { data.point[0], data.point[data.size], data.point[data.size * 2] };
        cl_double min[3] = { data.point[0], data.point[data.size], data.point[data.size * 2] };
        cl::Buffer max_buf(gpu.context, flag, sizeof(max[0]) * 3, max);
        cl::Buffer min_buf(gpu.context, flag, sizeof(min[0]) * 3, min);
    #endif

    #ifdef aos
        cl::Kernel maxmin(*gpu.program, "maxmin_aos");
        maxmin.setArg(0, (cl_uint) data.size);
        maxmin.setArg(1, gpu.point_buf);
        maxmin.setArg(2, max_buf);
        maxmin.setArg(3, min_buf);

        cl_uint localSize = 256;
        cl::NDRange local(localSize);
        cl::NDRange global((data.size + localSize - 1) / localSize * localSize);
        // Start running kernel
        gpu.cmdQueue.enqueueNDRangeKernel(maxmin, cl::NullRange, global, local);
        // Read result
        gpu.cmdQueue.enqueueReadBuffer(gpu.point_buf, CL_TRUE, 0, bytes, data.point);
        gpu.cmdQueue.enqueueReadBuffer(max_buf, CL_TRUE, 0, sizeof(max), max);
        gpu.cmdQueue.enqueueReadBuffer(min_buf, CL_TRUE, 0, sizeof(min), min);
    #endif
    #ifdef soa
        cl::Kernel maxmin(*gpu.program, "maxmin_soa");
        maxmin.setArg(0, (cl_uint) data.size);
        maxmin.setArg(1, gpu.point_buf);
        maxmin.setArg(2, max_buf);
        maxmin.setArg(3, min_buf);

        cl_uint localSize = 256;
        cl::NDRange local(localSize);
        cl::NDRange global((data.size + localSize - 1) / localSize * localSize);
        // Start running kernel
        gpu.cmdQueue.enqueueNDRangeKernel(maxmin, cl::NullRange, global, local);
        // Read result
        gpu.cmdQueue.enqueueReadBuffer(gpu.point_buf, CL_TRUE, 0, bytes, data.point);
        gpu.cmdQueue.enqueueReadBuffer(max_buf, CL_TRUE, 0, sizeof(max), max);
        gpu.cmdQueue.enqueueReadBuffer(min_buf, CL_TRUE, 0, sizeof(min), min);
    #endif

    cl_double AABB[3];
    AABB[0] = 0.5 * (min[0] + max[0]);
    AABB[1] = 0.5 * (min[1] + max[1]);
    AABB[2] = 0.5 * (min[2] + max[2]);

    cl::Buffer AABB_buf(gpu.context, flag, sizeof(AABB), AABB);

    #ifdef aos
        cl::Kernel center(*gpu.program, "center_aos");
        center.setArg(0, (cl_uint) data.size);
        center.setArg(1, gpu.point_buf);
        center.setArg(2, AABB_buf);

        // Start running kernel
        gpu.cmdQueue.enqueueNDRangeKernel(center, cl::NullRange, global, local);
        // Read result
        gpu.cmdQueue.enqueueReadBuffer(gpu.point_buf, CL_TRUE, 0, bytes, data.point);
    #endif
    #ifdef soa
        cl::Kernel center(*gpu.program, "center_soa");
        center.setArg(0, (cl_uint) data.size);
        center.setArg(1, gpu.point_buf);
        center.setArg(3, AABB_buf);

        // Start running kernel
        gpu.cmdQueue.enqueueNDRangeKernel(center, cl::NullRange, global, local);
        // Read result
        gpu.cmdQueue.enqueueReadBuffer(gpu.point_buf, CL_TRUE, 0, bytes, data.point);
    #endif

    // Get maximum possible rho
    cl_double diff[3];
    diff[0] = 0.5 * (max[0] - min[0]);
    diff[1] = 0.5 * (max[1] - min[1]);
    diff[2] = 0.5 * (max[2] - min[2]);

    double max_rho = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);

	return max_rho; 
}

// This function prepare the accumulator struct votes so that it will have sufficient memory for storing all votes.
// Allocate a buffer for the voting process.  You may need additional buffers for storing parameters (theta, rho, phi, etc.)  that are accessible from GPU
void prepareAccumulator(accumulator &votes, const double rho_max, const size_t n_theta, const size_t n_phi, const size_t n_rho, mygpu& gpu) {
	//cout << "\nDo not forget to implement prepareAccumulator! Line " << __LINE__ << "@ " << __FILE__ ; 
    votes.n_theta = n_theta;
    votes.n_phi = n_phi;
    votes.n_rho = n_rho;

    votes.d_theta = M_PI / (double) n_theta;
    votes.d_phi = M_PI / (2.0 * (double) (n_phi - 1));
    votes.d_rho = (2.0 * rho_max) / ((double) n_rho - 1);

    votes.rho_max = rho_max;

    votes.accum = new size_t [n_rho * n_theta * n_phi];
}

// This function release the allocated buffers for the accumulator votes.
void releaseAccumulator(accumulator &votes) {
	//cout << "\nDo not forget to implement releaseAccumulator! Line " << __LINE__ << "@ " << __FILE__; 
    delete [] votes.accum;
}

// This function conducts the Hough Transform to cast votes in the rho, theta, phi parametric space using OpenCL kernel functions.
void houghTransform(const mydata &data, accumulator &votes, mygpu& gpu) {
	//cout << "\nDo not forget to implement houghTransform! Line " << __LINE__ << "@ " << __FILE__; 
    cl_mem_flags flag = CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR;

    cl_uint bytes = sizeof(votes.accum[0]) * votes.n_rho * votes.n_theta * votes.n_phi;
    gpu.accum_buf = cl::Buffer(gpu.context, flag, bytes, votes.accum);

    #ifdef aos
        cl::Kernel hough(*gpu.program, "hough_aos");
        hough.setArg(0, (cl_uint) data.size);
        hough.setArg(1, gpu.point_buf);
        hough.setArg(2, (cl_uint) votes.n_theta);
        hough.setArg(3, (cl_uint) votes.n_phi);
        hough.setArg(4, (cl_uint) votes.n_rho);
        hough.setArg(5, (cl_double) votes.rho_max);
        hough.setArg(6, (cl_double) votes.d_theta);
        hough.setArg(7, (cl_double) votes.d_phi);
        hough.setArg(8, (cl_double) votes.d_rho);
        hough.setArg(9, gpu.accum_buf);
    
        cl_uint localSize = 4;
        cl_uint thetaSize = (votes.n_theta + localSize - 1) / localSize * localSize;
        cl_uint phiSize = (votes.n_phi + localSize - 1) / localSize * localSize;
        cl::NDRange local(localSize, localSize);
        cl::NDRange global(thetaSize, phiSize);
        // Start running kernel
        gpu.cmdQueue.enqueueNDRangeKernel(hough, cl::NullRange, global, local);
        // Read result
        gpu.cmdQueue.enqueueReadBuffer(gpu.accum_buf, CL_TRUE, 0, bytes, votes.accum);
    #endif 
    #ifdef soa
        cl::Kernel hough(*gpu.program, "hough_aos");
        hough.setArg(0, (cl_uint) data.size);
        hough.setArg(1, gpu.point_buf);
        hough.setArg(2, (cl_uint) votes.n_theta);
        hough.setArg(3, (cl_uint) votes.n_phi);
        hough.setArg(4, (cl_uint) votes.n_rho);
        hough.setArg(5, (cl_double) votes.rho_max);
        hough.setArg(6, (cl_double) votes.d_theta);
        hough.setArg(7, (cl_double) votes.d_phi);
        hough.setArg(8, (cl_double) votes.d_rho);
        hough.setArg(9, gpu.accum_buf);
    
        cl_uint localSize = 4;
        cl_uint thetaSize = (votes.n_theta + localSize - 1) / localSize * localSize;
        cl_uint phiSize = (votes.n_phi + localSize - 1) / localSize * localSize;
        cl::NDRange local(localSize, localSize);
        cl::NDRange global(thetaSize, phiSize);
        // Start running kernel
        gpu.cmdQueue.enqueueNDRangeKernel(hough, cl::NullRange, global, local);
        // Read result
        gpu.cmdQueue.enqueueReadBuffer(gpu.accum_buf, CL_TRUE, 0, bytes, votes.accum);
    #endif 
}

// find the largest vote and its corresponding Hough parameter
void identifyPlaneParameters(const accumulator& votes, houghPlanes &results, mygpu& gpu) {
	//cout << "\nDo not forget to implement identifyPlaneParameters! Line " << __LINE__ << "@ " << __FILE__; 
    results.votes = 0;

    for (size_t i = 0; i < votes.n_theta; i++) {
        for (size_t j = 0; j < votes.n_phi; j++) {
            for (size_t k = 0; k < votes.n_rho; k++) {
                size_t pos = i * votes.n_phi * votes.n_rho + j * votes.n_rho + k;

                if (results.votes < votes.accum[pos]) {
                    results.votes = votes.accum[pos];
                    results.theta = (double) i * votes.d_theta;
                    results.phi   = (double) j * votes.d_phi;
                    results.rho   = (double) k * votes.d_rho - votes.rho_max;
                }
            }
        }
    }
}

// This function de-allocate memory allocated for planes
void releaseHoughPlanes(houghPlanes &planes) {
	cout << "\nDo not forget to implement releaseHoughPlanes! Line " << __LINE__ << "@ " << __FILE__; 
}

// task 10 - release allocated memory for point-cloud data
void release(mydata& data) {
	cout << "\nDo not forget to implement release! Line " << __LINE__ << "@ " << __FILE__; 
}

// Task 8a - Output transformed point cloud data
bool outputPtxFile(const mydata& data, const houghPlanes &results, const accumulator& votes, const char *outputCloudData, mygpu& gpu) {
	cout << "\nRemeber to complete outputPtxFile! Line " << __LINE__ << "@ " << __FILE__; 
	ofstream outp(outputCloudData);
	if (!outp) return false; 

	// Output header of PLY file format
	outp << "ply\nformat ascii 1.0\nelement vertex " << 8 /* <-- Please replace 8 with the number of points from your data struct ... */
		<< "\nproperty float x\nproperty float y\nproperty float z"
		<< "\nproperty uchar red\nproperty uchar green\nproperty uchar blue"
		<< "\nend_header";

	// go through every point in your data struct and output x, y, z, R, G, B
    for (size_t i = 0; i < data.size; ++i) {
        #ifdef aos
            double point[3] = { data.point[i * 3], data.point[i * 3 + 1], data.point[i * 3 + 2] };
        #endif
        #ifdef soa
            double point[3] = { data.point[i], data.point[data.size + i], data.point[data.size * 2 + i] };
        #endif

        double cos_theta = cos(results.theta);
        double sin_theta = sin(results.theta);

        double cos_phi = cos(results.phi);
        double sin_phi = sin(results.phi);

        double rho = point[0] * cos_theta * sin_phi +
                     point[1] * sin_theta * sin_phi +
                     point[2] * cos_phi;

        outp << "\n" << point[0] << " " << point[1] << " "<< point[2];

        if (abs(rho - results.rho) < votes.d_rho) {
            outp << "64 64 255"; // Dark washed blue
        } else {
            outp << "128 128 128"; // Gray
        }
    }

	outp.close(); 

	return true;
}


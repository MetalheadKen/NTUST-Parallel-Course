#include <iostream>
#include <fstream>
#include "stopwatch.hpp"
#include "hw05.hpp"
#include <fstream>

using namespace std;

int main(int argc, char **argv) {
	if (argc != 6) {
		cerr << argc << endl; 
		cerr << "\n" << argv[0] << " [point cloud data] [n_theta] [n_phi] [n_rho] [output ptx filename]"; 
		return 255;
	}

	// task a: get command line arguments
	const char *inputCloudData = argv[1];
	const size_t n_theta = atoi(argv[2]); 
	const size_t n_phi = atoi(argv[3]);
	const size_t n_rho = atoi(argv[4]);
	const char *outputCloudData = argv[5];

	stopwatch t[4];

	// task b: prepare OpenCL device, command queue, and context.
	t[0].start(); 
	mygpu gpu; 
	prepareGPU(gpu); 

	// task c: read point-cloud data into OpenCL buffers
	mydata data;
	if (!readPointCloud(inputCloudData, data, gpu)) {
		cerr << "\nError reading the point-cloud data, file: " << inputCloudData; 
		return 254;
	}

	// Task d: Re-center the point-cloud data and get the maximum possible rho
	double rho_max = centerPointCloudToOrigin(data, gpu);
	t[0].stop(); 

	// task e: Prepare OpenCL buffer (filled with zero) for the accumulator and likely another buffer for storing voting parameters (rho, theta, phi, etc.)
	t[1].start(); 
	accumulator votes; 
	prepareAccumulator(votes, rho_max, n_theta, n_phi, n_rho, gpu); 

	// task f: Conduct Hough transform voting process.
	houghTransform(data, votes, gpu); 
	t[1].stop(); 

	// task g: Find the Hough plane with the most votes
	t[2].start(); 
	houghPlanes results; 
	identifyPlaneParameters(votes, results, gpu); 
	t[2].stop(); 

	// task h: Output PTX file with points on the plane with most votes in color of RGB(64, 64, 255) (dark washed blue) and 
	// points not on that plane with color of RGB(128, 128, 128) (gray).
	t[3].start();
	outputPtxFile(data, results, votes, outputCloudData, gpu); 
	t[3].stop(); 
	   	 
	// task 9 - output timing
	cout << "\n"; 
	for (auto &i : t) {
		cout << i.elapsedTime() << ", "; 
	}

	// task 10 - release all allocated memory
	releaseHoughPlanes(results); 
	releaseAccumulator(votes); 
	release(data); 
	releaseGPU(gpu);

	return 0;
}

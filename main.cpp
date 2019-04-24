#include <iostream>
#include <fstream>
#include "stopwatch.hpp"
#include "hw03.hpp"
#include <fstream>

using namespace std;

int main(int argc, char **argv) {
	if (argc != 7) {
		cerr << argc << endl; 
		cerr << "\n" << argv[0] << " [point cloud data] [n_theta] [n_phi] [n_rho] [threshold] [output ptx filename]"; 
		return 255;
	}

	// task 1 - get command line arguments
	const char *inputCloudData = argv[1];
	const size_t n_theta = atoi(argv[2]); 
	const size_t n_phi = atoi(argv[3]);
	const size_t n_rho = atoi(argv[4]);
	const double threshold = atoi(argv[5]); 
	const char *outputCloudData = argv[6];

	stopwatch t[4];

	// task 2 & 3 - allocate memory & read point-cloud data
	t[0].start(); 
	mydata data;
	if (!readPointCloud(inputCloudData, data)) {
		cerr << "\nError reading the point-cloud data, file: " << inputCloudData; 
		return 254;
	}

	// Task 4 - Center the point-cloud data and get the maximum possible rho
	double rho_max = centerPointCloudToOrigin(data);
	t[0].stop(); 

	// task 5 - prepare accumulator
	t[1].start(); 
	accumulator votes; 
	prepareAccumulator(votes, rho_max, n_theta, n_phi, n_rho); 

	// task 6 - do 3D standard Hough Transform
	houghTransform(data, votes); 
	t[1].stop(); 

	// task 7 - Find resultant plane parameters
	t[2].start(); 
	houghPlanes results; 
	identifyPlaneParameters(votes, threshold, results); 
	t[2].stop(); 

	// task 8 - Output ptx file
	t[3].start();
	outputPtxFile(data, results, votes, outputCloudData); 
	t[3].stop(); 
	   	 
	// task 9 - output timing
	cout << "\n"; 
	for (auto &i : t) {
		cout << i.elapsedTime() << ", "; 
	}
    cout << endl;

	// task 10 - release all allocated memory
	releaseHoughPlanes(results); 
	releaseAccumulator(votes); 
	release(data); 

	return 0;
}

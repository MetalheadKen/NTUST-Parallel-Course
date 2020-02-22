#include <iostream>
#include <fstream>
#include "stopwatch.hpp"
#include "hw01.hpp"
#include <fstream>

using namespace std;

// This is the main program that you should not change.

// task 4 - read the transformation matrix
// This function reads a 4x4 transformation matrix from the file specified by the filename
bool readTransformationMatrix(const char *filename, double (&M)[4][4]) {
	std::ifstream inp(filename); 
	if (!inp) return false; 

	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			inp >> M[i][j];
		}
	}

	inp.close();
	return true;
}

int main(int argc, char **argv) {
	if (argc != 6) {
		cerr << argc << endl; 
		cerr << "\n" << argv[0] << " [point cloud data] [transformation matrix] [AABB filename] [centroid filename] [output point cloud data]";
		return 255;
	}
	// task 1 - get command line arguments
	const char *inputCloudData = argv[1];
	const char *transformationMatrix = argv[2];
	const char *AABB = argv[3];
	const char *centroid = argv[4];
	const char *outputCloudData = argv[5];

	// task 2 & 3 - allocate memory & read point-cloud data
	mydata data;
	double M[4][4]; 
	if (!readPointCloud(inputCloudData, data)) {
		cerr << "\nError reading the point-cloud data, file: " << inputCloudData; 
		return 254;
	}
	// task 4 - read the transformation matrix
	if (!readTransformationMatrix(transformationMatrix, M)) {
		cerr << "\nError reading transformation matrix from file: " << transformationMatrix;
		return 253;
	}
	
	stopwatch t1, t2, t3, t4; 
	// Task 5 - Apply coordinate transformation on the point-cloud data
	t1.start();
	applyTransformation(M, data);
	t1.stop();

	// task 6 - FindAABB
	t2.start(); 
	findAABB(data);
	t2.stop(); 

	// task 7 - Find Centroid
	t3.start();
	findCentroid(data);
	t3.stop();

	// Task 8 - Output results to files.
	t4.start(); 
	if (!writePointCloud(outputCloudData, data)) {
		cerr << "\nFailed to write the output data: " << outputCloudData; 
	}
	if (!writeCentroid(centroid, data)) {
		cerr << "\nFailed to write the centroid data: " << centroid;
	}
	if (!writeAABB(AABB, data)) {
		cerr << "\nFailed to write the AABB data: " << AABB; 
	}
	t4.stop(); 

	// task 9 - output timing
	cout
		<< "\n" << t1.elapsedTime()
		<< "\t" << t2.elapsedTime()
		<< "\t" << t3.elapsedTime()
		<< "\t" << t4.elapsedTime() << endl; 

	// task 10 - release allocated memory
	release(data); 

	return 0;
}
#include "hw05.hpp"
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

// Prepare the needed GPU stuff (context, device, commandQueue)
void prepareGPU(mygpu& gpu) {
	cout << "\nDo not forget to implement prepareGPU! Line " << __LINE__ << "@ " << __FILE__; 
}

void releaseGPU(mygpu& gpu) {
	cout << "\nDo not forget to implement releaseGPU! Line " << __LINE__ << "@ " << __FILE__; 
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

	// some pts file has a 'D' character after the number of points...:(
	char pp = inp.peek(); 
	if (pp == 'D' || pp == 'd') inp >> pp; 

	// you should now allocate enough memory to store all the data...
	#ifdef aos
		// prepare AoS data struct
		cout << "\nAllocate for AoS data struct"; 
	#endif
	#ifdef soa
		// prepare SoA data struct
		cout << "\nAllocate for SoA data struct"; 
	#endif

	for (size_t i = 0; i < nData; ++i) {
		double x, y, z;
		string rest;		// some pts files have only intensity, and some have RGB values, so we treat them all as a long string...
		inp >> x >> y >> z;
		std::getline(inp, rest); 
		// if(i==0) cout << "\n" << x << ", " << y << ", " << z << ": " << rest; 
		// now you have coordinates in x, y, z, and additional information in the string rest.  You need to store them into your data struct...
		// we are going to discard the additional information in this assignment.
	}
	inp.close();

	return true;
}

// This function moves all points in the point cloud so that all points that was in the AABB of [xmin, ymin, zmin] - [xmax, ymax, zmax] 
// will be in the region of [-Lx, -Ly, -Lz] - [Lx, Ly, Lz], where Lx = 0.5*(xmax+xmin), Ly=0.5*(ymax+ymin), Lz=0.5*(zmax+zmin)
// Also, this function returns sqrt( (0.5*(xmax-xmin))**2 + (0.5*(ymax-ymin))**2 + (0.5*(zmax-zmin))**2) ).  This is the maximum possible rho.
// Do this with OpenCL GPU code
double centerPointCloudToOrigin(mydata &data, mygpu& gpu) {
	cout << "\nDo not forget to implement centerPointCloudToOrigin! Line " << __LINE__ << "@ " << __FILE__; 
	return 0.0; 
}

// This function prepare the accumulator struct votes so that it will have sufficient memory for storing all votes.
// Allocate a buffer for the voting process.  You may need additional buffers for storing parameters (theta, rho, phi, etc.)  that are accessible from GPU
void prepareAccumulator(accumulator &votes, const double rho_max, const size_t n_theta, const size_t n_phi, const size_t n_rho, mygpu& gpu) {
	cout << "\nDo not forget to implement prepareAccumulator! Line " << __LINE__ << "@ " << __FILE__ ; 
}

// This function release the allocated buffers for the accumulator votes.
void releaseAccumulator(accumulator &votes) {
	cout << "\nDo not forget to implement releaseAccumulator! Line " << __LINE__ << "@ " << __FILE__; 
}

// This function conducts the Hough Transform to cast votes in the rho, theta, phi parametric space using OpenCL kernel functions.
void houghTransform(const mydata &data, accumulator &votes, mygpu& gpu) {
	cout << "\nDo not forget to implement houghTransform! Line " << __LINE__ << "@ " << __FILE__; 
}

// find the largest vote and its corresponding Hough parameter
void identifyPlaneParameters(const accumulator& votes, houghPlanes &results, mygpu& gpu) {
	cout << "\nDo not forget to implement identifyPlaneParameters! Line " << __LINE__ << "@ " << __FILE__; 
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

	// The following outputs 8 points with different colors.  Please delete them after you have completed your outputPtxFile....
	outp << "\n0 0 0 255 0 0"; 	// (0, 0, 0) with color RGB(255, 0, 0)
	outp << "\n1 0 0 0 255 0";  // (1, 0, 0) with color RGB(0, 255, 0)
	outp << "\n1 1 0 0 0 255";  // (1, 1, 0) with color RGB(0, 0, 255)
	outp << "\n0 1 0 255 255 0"; // (0, 1, 0) with color RGB(255, 255, 0); 
	outp << "\n0 0 1 255 0 255"; 	// (0, 0, 1) with color RGB(255, 0, 255)
	outp << "\n1 0 1 0 255 255";  // (1, 0, 1) with color RGB(0, 255, 255)
	outp << "\n1 1 1 255 255 255";  // (1, 1, 1) with color RGB(255,255,255)
	outp << "\n0 1 1 0 0 0"; // (0, 1, 1) with color RGB(0, 0, 0); 

	// go through every point in your data struct and output x, y, z, R, G, B

	outp.close(); 

	return true;
}


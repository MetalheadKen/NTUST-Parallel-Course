#include "hw01.hpp"
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

// task 2 & 3 goes in here.
bool readPointCloud(const char *filename, mydata& data) {
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

	for (size_t i = 0; i < nData; ++i) {
		double x, y, z;
		string rest;		// some pts files have only intensity, and some have RGB values, so we treat them all as a long string...
		inp >> x >> y >> z;
		std::getline(inp, rest); 

		// now all the data for a single point is read in x, y, z and rest.  You should store it somewhere...
	}
	inp.close();

	return true;
}

// task 10 - release allocated memory
void release(mydata& data) {
	cout << "\nAlso need to de-allocate allocated memory...";
}

// task 5 - FindAABB
void findAABB(mydata& data) {
	cout << "\nRemeber to implement findAABB..."; 
}

// task 6 - Find Centroid
void findCentroid(mydata& data) {
	cout << "\nRemeber to implement findCentroid...";
}

// Task 7 - Apply coordinate transformation on the point-cloud data
void applyTransformation(const double(&M)[4][4], mydata& data) {
	cout << "\nRemeber to implement applyTransformation...";
}

// Task 8a - Output transformed point cloud data, returns true if success and false if failure.
bool writePointCloud(const char *outputFilename, const mydata& data) {
	cout << "\nRemeber to implement writePointCloud...";

	return true;
}

// Task 8b - Output centroid of the point-cloud, returns true if success and false if failure.
bool writeCentroid(const char *filename, const mydata &data) {
	cout << "\nRemeber to implement writeCentroid...";

	return true;
}

// Task 8c - Output AABB of the point-cloud
bool writeAABB(const char *filename, const mydata &data) {
	cout << "\nRemeber to implement writeAABB...";

	return true;
}

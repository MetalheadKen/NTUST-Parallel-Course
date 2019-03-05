#pragma once
#include <string>

struct mydata {
	// you need to determine how you want to store all your data in this data struct
};

// read the pts file and allocate all needed memory space here
bool readPointCloud(const char *filename, mydata& data);

// release all allocated memory, often at the end of the program
void release(mydata& data); 

// Applicated specified transformation in M matrix, and results are stored in the struct data.
void applyTransformation(const double (&M)[4][4], mydata& data);

// find AABB of the given point cloud data.  The results are also stored in the same struct.
void findAABB(mydata& data);

// find Centroid of the the given point cloud data.  The results are also stored in the same struct.
void findCentroid(mydata& data);

// write point cloud data in pts file format in the file specified by the filename.
bool writePointCloud(const char *filename, const mydata& data);

// write the centroid of the point cloud data in the file specified by the filename.
bool writeCentroid(const char *filename, const mydata &data); 

// write AABB of the point cloud in the file specified by the filename.
bool writeAABB(const char *filename, const mydata &data);


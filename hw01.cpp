#include "hw01.hpp"
#include <fstream>
#include <sstream>
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
    data.size = nData;

	// some pts file has a 'D' character after the number of points...:(
	char pp = inp.peek(); 
	if (pp == 'D' || pp == 'd') inp >> pp; 

	// you should now allocate enough memory to store all the data...
    data.data = new struct data [data.size];

	for (size_t i = 0; i < nData; ++i) {
		double x, y, z;
		string rest;		// some pts files have only intensity, and some have RGB values, so we treat them all as a long string...
		inp >> x >> y >> z;
		std::getline(inp, rest); 

		// now all the data for a single point is read in x, y, z and rest.  You should store it somewhere...
        data.data[i].point.x = x;
        data.data[i].point.y = y;
        data.data[i].point.z = z;

        std::istringstream in(rest);
        double intensity = -1.0;
        int r = -1, g = -1, b = -1;

        in >> intensity >> r >> g >> b;

        data.data[i].intensity = intensity;
        data.data[i].color.r = r;
        data.data[i].color.g = g;
        data.data[i].color.b = b;
	}
	inp.close();

	return true;
}

// task 10 - release allocated memory
void release(mydata& data) {
    delete [] data.data;
}

// task 5 - FindAABB
void findAABB(mydata& data) {
    data.max.x = data.max.y = data.max.z = DBL_MIN;
    data.min.x = data.min.y = data.min.z = DBL_MAX;

    for (size_t i = 0; i < data.size; i++) {
        if (data.data[i].point.x > data.max.x) data.max.x = data.data[i].point.x;
        if (data.data[i].point.y > data.max.y) data.max.y = data.data[i].point.y;
        if (data.data[i].point.z > data.max.z) data.max.z = data.data[i].point.z;

        if (data.data[i].point.x < data.min.x) data.min.x = data.data[i].point.x;
        if (data.data[i].point.y < data.min.y) data.min.y = data.data[i].point.y;
        if (data.data[i].point.z < data.min.z) data.min.z = data.data[i].point.z;
    }
}

// task 6 - Find Centroid
void findCentroid(mydata& data) {
    data.centroid.x = data.centroid.y = data.centroid.z = 0;

    for (size_t i = 0; i < data.size; i++) {
        data.centroid.x += data.data[i].point.x;
        data.centroid.y += data.data[i].point.y;
        data.centroid.z += data.data[i].point.z;
    }

    data.centroid.x /= data.size;
    data.centroid.y /= data.size;
    data.centroid.z /= data.size;
}

// Task 7 - Apply coordinate transformation on the point-cloud data
void applyTransformation(const double(&M)[4][4], mydata& data) {
    for (size_t i = 0; i < data.size; i++) {
        struct point trans;

        trans.x = data.data[i].point.x * M[0][0] + 
                  data.data[i].point.y * M[1][0] +
                  data.data[i].point.z * M[2][0] +
                  1.0 * M[3][0];

        trans.y = data.data[i].point.x * M[0][1] +
                  data.data[i].point.y * M[1][1] +
                  data.data[i].point.z * M[2][1] +
                  1.0 * M[3][1];

        trans.z = data.data[i].point.x * M[0][2] +
                  data.data[i].point.y * M[1][2] +
                  data.data[i].point.z * M[2][2] +
                  1.0 * M[3][2];

        data.data[i].point = trans;
    }
}

// Task 8a - Output transformed point cloud data, returns true if success and false if failure.
bool writePointCloud(const char *outputFilename, const mydata& data) {
    ofstream outp(outputFilename);
    if (!outp) {
       cerr << "\nError opening point cloud data file: " << outputFilename;
       return false;
    }

    outp << data.size << endl;
    for (size_t i = 0; i < data.size; i++) {
        outp << data.data[i].point.x << "\t" << data.data[i].point.y << "\t" << data.data[i].point.z;
        
        if (data.data[i].intensity != -1.0) outp << "\t" << data.data[i].intensity;
        if (data.data[i].color.r != -1) outp << "\t" << data.data[i].color.r;
        if (data.data[i].color.g != -1) outp << "\t" << data.data[i].color.g;
        if (data.data[i].color.b != -1) outp << "\t" << data.data[i].color.b;

        outp << endl;
    }
    outp.close();

	return true;
}

// Task 8b - Output centroid of the point-cloud, returns true if success and false if failure.
bool writeCentroid(const char *filename, const mydata &data) {
    ofstream outp(filename);
    if (!outp) {
        cerr << "\nError opening centroid data file: " << filename;
        return false;
    }

    outp << data.centroid.x << "\t" << data.centroid.y << "\t" << data.centroid.z;
    outp.close();

	return true;
}

// Task 8c - Output AABB of the point-cloud
bool writeAABB(const char *filename, const mydata &data) {
    ofstream outp(filename);
    if (!outp) {
        cerr << "\nError opening AABB data file: " << filename;
        return false;
    }

    outp << data.min.x << "\t" << data.min.y << "\t" << data.min.z << endl;
    outp << data.max.x << "\t" << data.max.y << "\t" << data.max.z;
    outp.close();

	return true;
}

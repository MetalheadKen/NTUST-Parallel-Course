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
    data.point = new struct point [data.size];
    data.color = new struct color [data.size];
    data.intensity = new double [data.size];

	for (size_t i = 0; i < nData; ++i) {
		double x, y, z;
		string rest;		// some pts files have only intensity, and some have RGB values, so we treat them all as a long string...
		inp >> x >> y >> z;
		std::getline(inp, rest); 

		// now all the data for a single point is read in x, y, z and rest.  You should store it somewhere...
        data.point[i].x = x;
        data.point[i].y = y;
        data.point[i].z = z;

        std::istringstream in(rest);
        double intensity = -1.0;
        int r = -1, g = -1, b = -1;

        in >> intensity >> r >> g >> b;

        data.intensity[i] = intensity;
        data.color[i].r = r;
        data.color[i].g = g;
        data.color[i].b = b;
	}
	inp.close();

	return true;
}

// task 10 - release allocated memory
void release(mydata& data) {
    delete [] data.point;
    delete [] data.color;
    delete [] data.intensity;
}

// task 5 - FindAABB
void findAABB(mydata& data) {
    struct point max, min;

    max.x = max.y = max.z = DBL_MIN;
    min.x = min.y = min.z = DBL_MAX;

    for (size_t i = 0; i < data.size; i++) {
        struct point *ppoint = data.point + i;

        if (ppoint->x > max.x) max.x = ppoint->x;
        if (ppoint->y > max.y) max.y = ppoint->y;
        if (ppoint->z > max.z) max.z = ppoint->z;

        if (ppoint->x < min.x) min.x = ppoint->x;
        if (ppoint->y < min.y) min.y = ppoint->y;
        if (ppoint->z < min.z) min.z = ppoint->z;
    }

    data.max = max;
    data.min = min;
}

// task 6 - Find Centroid
void findCentroid(mydata& data) {
    struct point centroid;
    
    centroid.x = centroid.y = centroid.z = 0;

    for (size_t i = 0; i < data.size; i++) {
        struct point *ppoint = data.point + i;
        
        centroid.x += ppoint->x;
        centroid.y += ppoint->y;
        centroid.z += ppoint->z;
    }

    centroid.x /= data.size;
    centroid.y /= data.size;
    centroid.z /= data.size;

    data.centroid = centroid;
}

// Task 7 - Apply coordinate transformation on the point-cloud data
void applyTransformation(const double(&M)[4][4], mydata& data) {
    for (size_t i = 0; i < data.size; i++) {
        struct point trans;
        struct point *ppoint = data.point + i;

        trans.x = ppoint->x * M[0][0] +
                  ppoint->y * M[1][0] +
                  ppoint->z * M[2][0] +
                  1.0 * M[3][0];

        trans.y = ppoint->x * M[0][1] +
                  ppoint->y * M[1][1] +
                  ppoint->z * M[2][1] +
                  1.0 * M[3][1];

        trans.z = ppoint->x * M[0][2] +
                  ppoint->y * M[1][2] +
                  ppoint->z * M[2][2] +
                  1.0 * M[3][2];

        data.point[i] = trans;
    }
}

// Task 8a - Output transformed point cloud data, returns true if success and false if failure.
bool writePointCloud(const char *outputFilename, const mydata& data) {
    FILE *outp = fopen(outputFilename, "w");

    if (!outp) {
       cerr << "\nError opening point cloud data file: " << outputFilename;
       return false;
    }

    fprintf(outp, "%zu\n", data.size);
    for (size_t i = 0; i < data.size; i++) {
        fprintf(outp, "%g\t%g\t%g", data.point[i].x, data.point[i].y, data.point[i].z);
        
        if (data.intensity[i] != -1.0) fprintf(outp, "\t%g", data.intensity[i]);
        if (data.color[i].r != -1) fprintf(outp, "\t%hu", data.color[i].r);
        if (data.color[i].g != -1) fprintf(outp, "\t%hu", data.color[i].g);
        if (data.color[i].b != -1) fprintf(outp, "\t%hu", data.color[i].b);

        fprintf(outp, "\n");
    }
    fclose(outp);

	return true;
}

// Task 8b - Output centroid of the point-cloud, returns true if success and false if failure.
bool writeCentroid(const char *filename, const mydata &data) {
    FILE *outp = fopen(filename, "w");

    if (!outp) {
        cerr << "\nError opening centroid data file: " << filename;
        return false;
    }

    fprintf(outp, "%g\t%g\t%g", data.centroid.x, data.centroid.y, data.centroid.z);
    fclose(outp);

	return true;
}

// Task 8c - Output AABB of the point-cloud
bool writeAABB(const char *filename, const mydata &data) {
    FILE *outp = fopen(filename, "w");

    if (!outp) {
        cerr << "\nError opening AABB data file: " << filename;
        return false;
    }

    fprintf(outp, "%.0f\t%.0f\t%.0f\n", data.min.x, data.min.y, data.min.z);
    fprintf(outp, "%.0f\t%.0f\t%.0f", data.max.x, data.max.y, data.max.z);
    fclose(outp);

	return true;
}

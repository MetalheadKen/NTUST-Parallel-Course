#pragma once
#include <cstring>

struct mygpu {
}; 

void prepareGPU(mygpu& gpu); 
void releaseGPU(mygpu& gpu); 

struct mydata {
};
bool readPointCloud(const char *filename, mydata& data, mygpu& gpu);
double centerPointCloudToOrigin(mydata &data, mygpu &gpu);

struct accumulator {
};
void prepareAccumulator(accumulator& votes, const double rho_max, const size_t n_theta, const size_t n_phi, const size_t n_rho, mygpu& gpu);
void releaseAccumulator(accumulator& votes);

void houghTransform(const mydata &data, accumulator &votes, mygpu& gpu);

// This struct should store resultant plane parameters in descending order
struct houghPlanes {
};

void identifyPlaneParameters(const accumulator &votes, houghPlanes &planes, mygpu& gpu);
void releaseHoughPlanes(houghPlanes &planes);

bool outputPtxFile(const mydata& data, const houghPlanes &results, const accumulator &votes, const char *outputCloudData, mygpu &gpu);

void release(mydata& data); 

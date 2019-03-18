#pragma once
#include <cstring>

struct mydata {
};
bool readPointCloud(const char *filename, mydata& data);
double centerPointCloudToOrigin(mydata &data);

struct accumulator {
};
void prepareAccumulator(accumulator& votes, const double rho_max, const size_t n_theta, const size_t n_phi, const size_t n_rho);
void releaseAccumulator(accumulator& votes);

void houghTransform(const mydata &data, accumulator &votes);

// This struct should store resultant plane parameters in descending order
struct houghPlanes {
};

void identifyPlaneParameters(const accumulator &votes, const size_t threshold, houghPlanes &planes);
void releaseHoughPlanes(houghPlanes &planes);

bool outputPtxFile(const mydata& data, const houghPlanes &results, const accumulator &votes, const char *outputCloudData);

void release(mydata& data); 

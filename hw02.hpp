#pragma once
#include <cstring>
#include <cmath>
#include <vector>

struct point {
    double x;
    double y;
    double z;
};

struct color {
    int R;
    int G;
    int B;
};

struct mydata {
    size_t size;
    struct point *point;
};

bool readPointCloud(const char *filename, mydata& data);
double centerPointCloudToOrigin(mydata &data);

struct accumulator {
    size_t ***accum;

    size_t n_theta;
    size_t n_phi;
    size_t n_rho;

    double rho_max;

    double d_theta;
    double d_phi;
    double d_rho;
};

void prepareAccumulator(accumulator& votes, const double rho_max, const size_t n_theta, const size_t n_phi, const size_t n_rho);
void releaseAccumulator(accumulator& votes);

void houghTransform(const mydata &data, accumulator &votes);

// This struct should store resultant plane parameters in descending order
struct houghPlane {
    size_t votes;

    double theta;
    double phi;
    double rho;
};

struct houghPlanes {
    std::vector<houghPlane> planes;
};

void identifyPlaneParameters(const accumulator &votes, const size_t threshold, houghPlanes &planes);
void releaseHoughPlanes(houghPlanes &planes);

bool outputPtxFile(const mydata& data, const houghPlanes &results, const accumulator &votes, const char *outputCloudData);

void release(mydata& data); 

#pragma once
#include <cstring>
#include <cmath>
#include "YoUtil.hpp"

struct mygpu {
    YoUtil::GPU *gpu;
    cl::CommandQueue cmdQueue;
    cl::Device device;
    cl::Context context;
    cl::Program *program;

    cl::Buffer point_buf;
    cl::Buffer accum_buf;
};

void prepareGPU(mygpu& gpu); 
void releaseGPU(mygpu& gpu); 

struct mydata {
    size_t size;
    cl_double *point;
};

bool readPointCloud(const char *filename, mydata& data, mygpu& gpu);
double centerPointCloudToOrigin(mydata &data, mygpu &gpu);

struct accumulator {
    size_t *accum;

    size_t n_theta;
    size_t n_phi;
    size_t n_rho;

    cl_double rho_max;

    cl_double d_theta;
    cl_double d_phi;
    cl_double d_rho;
};
void prepareAccumulator(accumulator& votes, const double rho_max, const size_t n_theta, const size_t n_phi, const size_t n_rho, mygpu& gpu);
void releaseAccumulator(accumulator& votes);

void houghTransform(const mydata &data, accumulator &votes, mygpu& gpu);

// This struct should store resultant plane parameters in descending order

struct houghPlanes {
    size_t votes;

    double theta;
    double phi;
    double rho;
};

void identifyPlaneParameters(const accumulator &votes, houghPlanes &planes, mygpu& gpu);
void releaseHoughPlanes(houghPlanes &planes);

bool outputPtxFile(const mydata& data, const houghPlanes &results, const accumulator &votes, const char *outputCloudData, mygpu &gpu);

void release(mydata& data); 

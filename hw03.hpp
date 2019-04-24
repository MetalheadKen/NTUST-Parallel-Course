#pragma once
#include <cstring>
#include <cmath>
#include <vector>

#if defined(aos) || defined(soa) 
    #define p1
#elif defined(p1) || defined(p2) || defined(p3) || defined(p4) || defined(p5) || defined(p6)
    #define aos
#elif defined(best)
    #define aos
    #define p6
#endif

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

#ifdef aos
    struct mydata {
        size_t size;
        struct point *point;
    };
#elif defined(soa)
    struct mydata {
        size_t size;

        struct {
            double *x;
            double *y;
            double *z;
        } point;
    };
#else
    #error "There are no struct mydata here"
#endif

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

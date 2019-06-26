#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

void Atomic_Max(volatile __global cl_double *src, cl_double val) {
    union {
        unsigned int intVal;
        double doubleVal;
    } newVal, prevVal;

    do {
        prevVal.doubleVal = *src;
        newVal.doubleVal = max(prevVal.doubleVal, val);
    } while (atomic_cmpxchg((volatile __global unsigned int *) src, prevVal.intVal, newVal.intVal != prevVal.intVal);
}

void Atomic_Min(volatile __global cl_double *src, cl_double val) {
    union {
        unsigned int intVal;
        double doubleVal;
    } newVal, prevVal;

    do {
        prevVal.doubleVal = *src;
        newVal.doubleVal = min(prevVal.doubleVal, val);
    } while (atomic_cmpxchg((volatile __global unsigned int *) src, prevVal.intVal, newVal.intVal != prevVal.intVal);
}

__kernel
double maxmin_aos(cl_uint nData, __global cl_double *point, __global cl_double *max, __global cl_double *min) {
    cl_uint i = get_global_id(0);

    if (i < nData) {
        if (point[i * 3    ] > max[0])      { Atomic_Max(&max[0], point[i * 3]); }
        else if (point[i * 3    ] < min[0]) { Atomic_Min(&min[0], point[i * 3]); }

        if (point[i * 3 + 1] > max[1])      { Atomic_Max(&max[1], point[i * 3 + 1]); }
        else if (point[i * 3 + 1] < min[1]) { Atomic_Min(&min[1], point[i * 3 + 1]); }
        
        if (point[i * 3 + 2] > max[2])      { Atomic_Max(&max[2], point[i * 3 + 2]); }
        else if (point[i * 3 + 2] < min[2]) { Atomic_Min(&min[2], point[i * 3 + 2]); }
    }   
}

__kernel
void maxmin_soa(cl_uint nData, __global cl_double *point, __global cl_double *max, __global cl_double *min) {
    cl_uint i = get_global_id(0);

    if (i < nData) {
        if (point[i] > max[0])      { Atomic_Max(&max[0], point[i]); }
        else if (point[i] < min[0]) { Atomic_Min(&min[0], point[i]); }

        if (point[nData + i] > max[1])      { Atomic_Max(&max[1], point[nData + i]); }
        else if (point[nData + i] < min[1]) { Atomic_Min(&min[1], point[nData + i]); }
        
        if (point[nData * 2 + i] > max[2])      { Atomic_Max(&max[2], point[nData * 2 + i]); }
        else if (point[nData * 2 + i] < min[2]) { Atomic_Min(&min[2], point[nData * 2 + i]); }
    }
}

__kernel
void center_aos(cl_uint nData, __global cl_double *point, __global cl_double *AABB) {
    cl_uint i = get_global_id(0);

    if (i < nData) {
        point[i * 3] -= AABB[0];
        point[i * 3 + 1] -= AABB[1];
        point[i * 3 + 2] -= AABB[2];
    }
}

__kernel
void center_soa(cl_uint nData, __global cl_double *point, __global cl_double *AABB) {
    cl_uint i = get_global_id(0);

    if (i < nData) {
        point[i] -= AABB[0];
        point[nData + i] -= AABB[1];
        point[nData * 2 + i] -= AABB[2];
    }
}

void hough_aos(cl_uint nData, __global cl_double *point, cl_uint n_theta, cl_uint n_phi, cl_uint n_rho, 
        cl_double rho_max, cl_double d_theta, cl_double d_phi, cl_double d_rho, __global cl_uint *accum) {
    cl_double theta, cos_theta, sin_theta, cos_phi, sin_phi, rho;
    cl_uint i = get_global_id(0);
    cl_uint j = get_global_id(1);

    if (i < n_theta) {
        if (j < n_phi) {
            for (cl_uint pos = 0; pos < nData; ++pos) {
                theta = (cl_double) i * votes.d_theta;
                cos_theta = cos(theta);
                sin_theta = sin(theta);

                phi = (cl_double) j * votes.d_phi;
                cos_phi = cos(phi);
                sin_phi = sin(phi);

                rho = data.point[pos * 3    ] * cos_theta * sin_phi + 
                      data.point[pos * 3 + 1] * sin_theta * sin_phi + 
                      data.point[pos * 3 + 2] * cos_phi;

                cl_uint k = (cl_uint) ((rho + votes.rho_max) / votes.d_rho);
                votes.accum[i * n_phi * n_rho + j * n_rho + k]++;
            }
        }
    }
}

void hough_soa(cl_uint nData, __global cl_double *point, cl_uint n_theta, cl_uint n_phi, cl_uint n_rho, 
        cl_double rho_max, cl_double d_theta, cl_double d_phi, cl_double d_rho, __global cl_uint *accum) {
    cl_double theta, cos_theta, sin_theta, cos_phi, sin_phi, rho;
    cl_uint i = get_global_id(0);
    cl_uint j = get_global_id(1);

    if (i < n_theta) {
        if (j < n_phi) {
            for (cl_uint pos = 0; pos < nData; ++pos) {
                theta = (cl_double) i * votes.d_theta;
                cos_theta = cos(theta);
                sin_theta = sin(theta);

                phi = (cl_double) j * votes.d_phi;
                cos_phi = cos(phi);
                sin_phi = sin(phi);

                rho = data.point[pos            ] * cos_theta * sin_phi + 
                      data.point[nData + pos    ] * sin_theta * sin_phi + 
                      data.point[nData * 2 + pos] * cos_phi;

                cl_uint k = (cl_uint) ((rho + votes.rho_max) / votes.d_rho);
                votes.accum[i * n_phi * n_rho + j * n_rho + k]++;
            }
        }
    }
}


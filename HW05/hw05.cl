#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

void Atomic_Max(volatile __global double *src, double val) {
    union {
        unsigned int intVal;
        double doubleVal;
    } newVal, prevVal;

    do {
        prevVal.doubleVal = *src;
        newVal.doubleVal = max(prevVal.doubleVal, val);
    } while (atomic_cmpxchg((volatile __global unsigned int *) src, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

void Atomic_Min(volatile __global double *src, double val) {
    union {
        unsigned int intVal;
        double doubleVal;
    } newVal, prevVal;

    do {
        prevVal.doubleVal = *src;
        newVal.doubleVal = min(prevVal.doubleVal, val);
    } while (atomic_cmpxchg((volatile __global unsigned int *) src, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

__kernel
void maxmin_aos(uint nData, __global double *point, __global double *max, __global double *min) {
    uint i = get_global_id(0);

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
void maxmin_soa(uint nData, __global double *point, __global double *max, __global double *min) {
    uint i = get_global_id(0);

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
void center_aos(uint nData, __global double *point, __global double *AABB) {
    uint i = get_global_id(0);

    if (i < nData) {
        point[i * 3] -= AABB[0];
        point[i * 3 + 1] -= AABB[1];
        point[i * 3 + 2] -= AABB[2];
    }
}

__kernel
void center_soa(uint nData, __global double *point, __global double *AABB) {
    uint i = get_global_id(0);

    if (i < nData) {
        point[i] -= AABB[0];
        point[nData + i] -= AABB[1];
        point[nData * 2 + i] -= AABB[2];
    }
}

__kernel
void hough_aos(uint nData, __global double *point, uint n_theta, uint n_phi, uint n_rho, 
        double rho_max, double d_theta, double d_phi, double d_rho, __global uint *accum) {
    double theta, phi, rho, cos_theta, sin_theta, cos_phi, sin_phi;
    uint i = get_global_id(0);
    uint j = get_global_id(1);

    if (i < n_theta) {
        if (j < n_phi) {
            for (uint pos = 0; pos < nData; ++pos) {
                theta = (double) i * d_theta;
                cos_theta = cos(theta);
                sin_theta = sin(theta);

                phi = (double) j * d_phi;
                cos_phi = cos(phi);
                sin_phi = sin(phi);

                rho = point[pos * 3    ] * cos_theta * sin_phi + 
                      point[pos * 3 + 1] * sin_theta * sin_phi + 
                      point[pos * 3 + 2] * cos_phi;

                uint k = (uint) ((rho + rho_max) / d_rho);
                uint index = i * n_phi * n_rho + j * n_rho + k;

                if (index < n_theta * n_phi * n_rho) {
                    accum[i * n_phi * n_rho + j * n_rho + k]++;
                }
            }
        }
    }
}

__kernel
void hough_soa(uint nData, __global double *point, uint n_theta, uint n_phi, uint n_rho, 
        double rho_max, double d_theta, double d_phi, double d_rho, __global uint *accum) {
    double theta, phi, rho, cos_theta, sin_theta, cos_phi, sin_phi;
    uint i = get_global_id(0);
    uint j = get_global_id(1);

    if (i < n_theta) {
        if (j < n_phi) {
            for (uint pos = 0; pos < nData; ++pos) {
                theta = (double) i * d_theta;
                cos_theta = cos(theta);
                sin_theta = sin(theta);

                phi = (double) j * d_phi;
                cos_phi = cos(phi);
                sin_phi = sin(phi);

                rho = point[pos            ] * cos_theta * sin_phi + 
                      point[nData + pos    ] * sin_theta * sin_phi + 
                      point[nData * 2 + pos] * cos_phi;

                uint k = (uint) ((rho + rho_max) / d_rho);
                uint index = i * n_phi * n_rho + j * n_rho + k;

                if (index < n_theta * n_phi * n_rho) {
                    accum[i * n_phi * n_rho + j * n_rho + k]++;
                }
            }
        }
    }
}


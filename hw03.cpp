#include "hw03.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

// Read the original point cloud data
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
#ifdef aos
    // prepare AoS data struct
    data.point = new struct point [data.size];
#endif
#ifdef soa
    // prepare SoA data struct
    data.point.x = new double [data.size];
    data.point.y = new double [data.size];
    data.point.z = new double [data.size];
#endif

	for (size_t i = 0; i < nData; ++i) {
		double x, y, z;
		string rest;		// some pts files have only intensity, and some have RGB values, so we treat them all as a long string...
		inp >> x >> y >> z;
		std::getline(inp, rest); 
		// if(i==0) cout << "\n" << x << ", " << y << ", " << z << ": " << rest; 
		// now you have coordinates in x, y, z, and additional information in the string rest.  You need to store them into your data struct...
		// we are going to discard the additional information in this assignment.
    #ifdef aos
        data.point[i].x = x;
        data.point[i].y = y;
        data.point[i].z = z;
    #endif
    #ifdef soa
        data.point.x[i] = x;
        data.point.y[i] = y;
        data.point.z[i] = z;
    #endif
	}
	inp.close();

	return true;
}

// This function moves all points in the point cloud so that all points that was in the AABB of [xmin, ymin, zmin] - [xmax, ymax, zmax] 
// will be in the region of [-Lx, -Ly, -Lz] - [Lx, Ly, Lz], where Lx = 0.5*(xmax+xmin), Ly=0.5*(ymax+ymin), Lz=0.5*(zmax+zmin)
// Also, this function returns sqrt( (0.5*(xmax-xmin))**2 + (0.5*(ymax-ymin))**2 + (0.5*(zmax-zmin))**2) ).  This is the maximum possible rho.
double centerPointCloudToOrigin(mydata &data) {
    // Get AABB
#ifdef aos
    struct point max = { data.point[0].x, data.point[0].y, data.point[0].z };
    struct point min = { data.point[0].x, data.point[0].y, data.point[0].z };
#endif
#ifdef soa
    struct point max = { data.point.x[0], data.point.y[0], data.point.z[0] };
    struct point min = { data.point.x[0], data.point.y[0], data.point.z[0] };
#endif

    for (size_t i = 1; i < data.size; ++i) {
    #ifdef aos
        if (data.point[i].x > max.x)      { max.x = data.point[i].x; }
        else if (data.point[i].x < min.x) { min.x = data.point[i].x; }

        if (data.point[i].y > max.y)      { max.y = data.point[i].y; }
        else if (data.point[i].y < min.y) { min.y = data.point[i].y; }
        
        if (data.point[i].z > max.z)      { max.z = data.point[i].z; }
        else if (data.point[i].z < min.z) { min.z = data.point[i].z; }
    #endif
    #ifdef soa
        if (data.point.x[i] > max.x)      { max.x = data.point.x[i]; }
        else if (data.point.x[i] < min.x) { min.x = data.point.x[i]; }

        if (data.point.y[i] > max.y)      { max.y = data.point.y[i]; }
        else if (data.point.y[i] < min.y) { min.y = data.point.y[i]; }
        
        if (data.point.z[i] > max.z)      { max.z = data.point.z[i]; }
        else if (data.point.z[i] < min.z) { min.z = data.point.z[i]; }
    #endif
    }

    // Center the point-cloud data so that the center of the AABB is at the origin
    struct point AABB;
    AABB.x = 0.5 * (min.x + max.x);
    AABB.y = 0.5 * (min.y + max.y);
    AABB.z = 0.5 * (min.z + max.z);

    for (size_t i = 0; i < data.size; ++i) {
    #ifdef aos
        data.point[i].x -= AABB.x;
        data.point[i].y -= AABB.y;
        data.point[i].z -= AABB.z;
    #endif
    #ifdef soa 
        data.point.x[i] -= AABB.x;
        data.point.y[i] -= AABB.y;
        data.point.z[i] -= AABB.z;
    #endif
    }

    // Get maximum possible rho
    struct point diff;
    diff.x = 0.5 * (max.x - min.x);
    diff.y = 0.5 * (max.y - min.y);
    diff.z = 0.5 * (max.z - min.z);

    double max_rho = sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);

	return max_rho; 
}

// This function prepare the accumulator struct votes so that it will have sufficient memory for storing all votes.
void prepareAccumulator(accumulator &votes, const double rho_max, const size_t n_theta, const size_t n_phi, const size_t n_rho) {
    votes.n_theta = n_theta;
    votes.n_phi = n_phi;
    votes.n_rho = n_rho;

    votes.d_theta = M_PI / (double) n_theta;
    votes.d_phi = M_PI / (2.0 * ((double) n_phi - 1));
    votes.d_rho = (2.0 * rho_max) / ((double) n_rho - 1);

    votes.rho_max = rho_max;

#ifdef p1
    size_t ***accum = new size_t **[n_rho]();

    for (size_t i = 0; i < n_rho; ++i) {
        accum[i] = new size_t *[n_theta]();
        for (size_t j = 0; j < n_theta; ++j) {
            accum[i][j] = new size_t [n_phi]();
        }
    }
    votes.accum = accum;
#endif
#ifdef p2
    size_t ***accum = new size_t **[n_rho]();

    for (size_t i = 0; i < n_rho; ++i) {
        accum[i] = new size_t *[n_phi]();
        for (size_t j = 0; j < n_phi; ++j) {
            accum[i][j] = new size_t [n_theta]();
        }
    }
    votes.accum = accum;
#endif
#ifdef p3
    size_t ***accum = new size_t **[n_theta]();

    for (size_t i = 0; i < n_theta; ++i) {
        accum[i] = new size_t *[n_rho]();
        for (size_t j = 0; j < n_rho; ++j) {
            accum[i][j] = new size_t [n_phi]();
        }
    }
    votes.accum = accum;
#endif
#ifdef p4
    size_t ***accum = new size_t **[n_theta]();

    for (size_t i = 0; i < n_theta; ++i) {
        accum[i] = new size_t *[n_phi]();
        for (size_t j = 0; j < n_phi; ++j) {
            accum[i][j] = new size_t [n_rho]();
        }
    }
    votes.accum = accum;
#endif
#ifdef p5
    size_t ***accum = new size_t **[n_phi]();

    for (size_t i = 0; i < n_phi; ++i) {
        accum[i] = new size_t *[n_rho]();
        for (size_t j = 0; j < n_rho; ++j) {
            accum[i][j] = new size_t [n_theta]();
        }
    }
    votes.accum = accum;
#endif
#ifdef p6
    size_t ***accum = new size_t **[n_phi]();

    for (size_t i = 0; i < n_phi; ++i) {
        accum[i] = new size_t *[n_theta]();
        for (size_t j = 0; j < n_theta; ++j) {
            accum[i][j] = new size_t [n_rho]();
        }
    }
    votes.accum = accum;
#endif
}

// This function release the allocated memory for the accumulator votes.
void releaseAccumulator(accumulator &votes) {
#ifdef p1
    for (size_t i = 0; i < votes.n_rho; ++i) {
        for (size_t j = 0; j < votes.n_theta; ++j) {
            delete [] votes.accum[i][j];
        }
        delete [] votes.accum[i];
    }
    delete [] votes.accum;
#endif
#ifdef p2
    for (size_t i = 0; i < votes.n_rho; ++i) {
        for (size_t j = 0; j < votes.n_phi; ++j) {
            delete [] votes.accum[i][j];
        }
        delete [] votes.accum[i];
    }
    delete [] votes.accum;
#endif
#ifdef p3
    for (size_t i = 0; i < votes.n_theta; ++i) {
        for (size_t j = 0; j < votes.n_rho; ++j) {
            delete [] votes.accum[i][j];
        }
        delete [] votes.accum[i];
    }
    delete [] votes.accum;
#endif
#ifdef p4
    for (size_t i = 0; i < votes.n_theta; ++i) {
        for (size_t j = 0; j < votes.n_phi; ++j) {
            delete [] votes.accum[i][j];
        }
        delete [] votes.accum[i];
    }
    delete [] votes.accum;
#endif
#ifdef p5
    for (size_t i = 0; i < votes.n_phi; ++i) {
        for (size_t j = 0; j < votes.n_rho; ++j) {
            delete [] votes.accum[i][j];
        }
        delete [] votes.accum[i];
    }
    delete [] votes.accum;
#endif
#ifdef p6
    for (size_t i = 0; i < votes.n_phi; ++i) {
        for (size_t j = 0; j < votes.n_theta; ++j) {
            delete [] votes.accum[i][j];
        }
        delete [] votes.accum[i];
    }
    delete [] votes.accum;
#endif
}

// This function conducts the Hough Transform to cast votes in the rho, theta, phi parametric space.
void houghTransform(const mydata &data, accumulator &votes) {
    double theta, phi, rho;
    double cos_theta, sin_theta;
    double cos_phi, sin_phi;

    for (size_t i = 0; i < votes.n_theta; ++i) {
        theta = (double) i * votes.d_theta;
        cos_theta = cos(theta);
        sin_theta = sin(theta);

        for (size_t j = 0; j < votes.n_phi; ++j) {
            phi = (double) j * votes.d_phi;
            cos_phi = cos(phi);
            sin_phi = sin(phi);

            for (size_t pos = 0; pos < data.size; ++pos) {
            #ifdef aos
                rho = data.point[pos].x * cos_theta * sin_phi + 
                      data.point[pos].y * sin_theta * sin_phi + 
                      data.point[pos].z * cos_phi;
            #endif
            #ifdef soa
                rho = data.point.x[pos] * cos_theta * sin_phi + 
                      data.point.y[pos] * sin_theta * sin_phi + 
                      data.point.z[pos] * cos_phi;
            #endif

                int k = (int) ((rho + votes.rho_max) / votes.d_rho);
            #ifdef p1
                votes.accum[k][i][j]++;
            #endif
            #ifdef p2
                votes.accum[k][j][i]++;
            #endif
            #ifdef p3
                votes.accum[i][k][j]++;
            #endif
            #ifdef p4
                votes.accum[i][j][k]++;
            #endif
            #ifdef p5
                votes.accum[j][k][i]++;
            #endif
            #ifdef p6
                votes.accum[j][i][k]++;
            #endif
            }
        }
    }
}

static inline bool compare(houghPlane planeA, houghPlane planeB) {
    return planeA.votes > planeB.votes;
}

// find votes that are larger than threshold and store its parameters into the results data struct
void identifyPlaneParameters(const accumulator& votes, const size_t threshold, houghPlanes &results) {
    for (size_t i = 0; i < votes.n_theta; ++i) {
        for (size_t j = 0; j < votes.n_phi; ++j) {
            for (size_t k = 0; k < votes.n_rho; ++k) {
            #ifdef p1
                if (votes.accum[k][i][j] > threshold) {
                    houghPlane plane;

                    plane.theta = (double) i * votes.d_theta;
                    plane.phi = (double) j * votes.d_phi;
                    plane.rho = (double) k * votes.d_rho - votes.rho_max;
                    plane.votes = votes.accum[k][i][j];

                    results.planes.push_back(plane);
                }
            #endif
            #ifdef p2
                if (votes.accum[k][j][i] > threshold) {
                    houghPlane plane;

                    plane.theta = (double) i * votes.d_theta;
                    plane.phi = (double) j * votes.d_phi;
                    plane.rho = (double) k * votes.d_rho - votes.rho_max;
                    plane.votes = votes.accum[k][j][i];

                    results.planes.push_back(plane);
                }
            #endif
            #ifdef p3
                if (votes.accum[i][k][j] > threshold) {
                    houghPlane plane;

                    plane.theta = (double) i * votes.d_theta;
                    plane.phi = (double) j * votes.d_phi;
                    plane.rho = (double) k * votes.d_rho - votes.rho_max;
                    plane.votes = votes.accum[i][k][j];

                    results.planes.push_back(plane);
                }
            #endif
            #ifdef p4
                if (votes.accum[i][j][k] > threshold) {
                    houghPlane plane;

                    plane.theta = (double) i * votes.d_theta;
                    plane.phi = (double) j * votes.d_phi;
                    plane.rho = (double) k * votes.d_rho - votes.rho_max;
                    plane.votes = votes.accum[i][j][k];

                    results.planes.push_back(plane);
                }
            #endif
            #ifdef p5
                if (votes.accum[j][k][i] > threshold) {
                    houghPlane plane;

                    plane.theta = (double) i * votes.d_theta;
                    plane.phi = (double) j * votes.d_phi;
                    plane.rho = (double) k * votes.d_rho - votes.rho_max;
                    plane.votes = votes.accum[j][k][i];

                    results.planes.push_back(plane);
                }
            #endif
            #ifdef p6
                if (votes.accum[j][i][k] > threshold) {
                    houghPlane plane;

                    plane.theta = (double) i * votes.d_theta;
                    plane.phi = (double) j * votes.d_phi;
                    plane.rho = (double) k * votes.d_rho - votes.rho_max;
                    plane.votes = votes.accum[j][i][k];

                    results.planes.push_back(plane);
                }
            #endif
            }
        }
    }

    // sort planes vector ordered by votes in descending order
    sort(results.planes.begin(), results.planes.end(), compare);
}

// This function de-allocate memory allocated for planes
void releaseHoughPlanes(houghPlanes &planes) {
    // vector will automatically free allocated memory
}

// task 10 - release allocated memory for point-cloud data
void release(mydata& data) {
#ifdef aos
    delete [] data.point;
#endif
#ifdef soa
    delete [] data.point.x;
    delete [] data.point.y;
    delete [] data.point.z;
#endif
}

// Task 8a - Output transformed point cloud data
bool outputPtxFile(const mydata& data, const houghPlanes &results, const accumulator& votes, const char *outputCloudData) {
	ofstream outp(outputCloudData);
	if (!outp) return false; 

	// Output header of PLY file format
	outp << "ply\nformat ascii 1.0\nelement vertex " << data.size /* <-- Please replace 8 with the number of points from your data struct ... */
		<< "\nproperty float x\nproperty float y\nproperty float z"
		<< "\nproperty uchar red\nproperty uchar green\nproperty uchar blue"
		<< "\nend_header";

    // go through every point in your data struct and output x, y, z, R, G, B
    for (size_t i = 0; i < data.size; ++i) {
    #ifdef aos
        struct point &ppoint = data.point[i];
    #endif
    #ifdef soa
        struct point ppoint = { data.point.x[i], data.point.y[i], data.point.z[i] };
    #endif
    struct color color = { 32, 32, 32 }; // Gray color

        for (houghPlane plane : results.planes) {
            double cos_theta = cos(plane.theta);
            double sin_theta = sin(plane.theta);

            double cos_phi = cos(plane.phi);
            double sin_phi = sin(plane.phi);

            double rho = ppoint.x * cos_theta * sin_phi +
                         ppoint.y * sin_theta * sin_phi +
                         ppoint.z * cos_phi;

            if (abs(rho - plane.rho) < votes.d_rho) {
                color.R = 255 * cos_theta * sin_phi;
                color.G = 255 * sin_theta * sin_phi;
                color.B = 255 * cos_phi;
                break;
            }
        }
        outp << "\n" << ppoint.x << " " << ppoint.y << " "<< ppoint.z;
        outp << " " << color.R << " " << color.G << " " << color.B;
    }

	outp.close(); 

	return true;
}

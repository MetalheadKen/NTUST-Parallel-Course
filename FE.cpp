#include "FE.hpp"
#include "PNGio.hpp"
#include <cmath>
#include <iostream>
using std::cout; 
// 
int readImage(const char *filename, inputImage &png) {
    cout << "\nCheck readImage in " << __FILE__ << " @ " << __LINE__; 
    // If you are okay with how image data is stored (defined in PNGio.cpp), you can simply use the following line of code for this function
    return pngRead(filename, png); 
}

void toGrayScale(const inputImage& input, outputImage &output) {
    output.width = input.width; 
    output.height = input.height; 
    output.pixels = new png_byte[output.width * output.height]; 
    for(int y=0;y<input.height;++y) {
        for(int x=0;x<input.width;++x) {
            unsigned char r = input.row_pointers[y][4 * x];
            unsigned char g = input.row_pointers[y][4 * x + 1];
            unsigned char b = input.row_pointers[y][4 * x + 2];

            // Turn RGBA into gray-scale image
            output.pixels[y*input.width+x] = 0.299 * r + 0.587 * g + 0.114 * b; 
        }
    }
}

// https://en.wikipedia.org/wiki/Sobel_operator
void toEdge(const outputImage &input, outputImage &output) {
    const short GX[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
    const short GY[3][3] = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };
    
    output.width = input.width; 
    output.height = input.height; 
    output.pixels = new png_byte[output.width * output.height]; 

    for(png_uint_32 y = 1; y < input.height - 1; ++y) {
        for(png_uint_32 x = 1; x < input.width - 1; ++x) {
            short mask[3][3], gx = 0, gy = 0;

            mask[0][0] = input.pixels[(y - 1) * input.width + (x - 1)];
            mask[0][1] = input.pixels[(y - 1) * input.width + (x    )];
            mask[0][2] = input.pixels[(y - 1) * input.width + (x + 1)];
            mask[1][0] = input.pixels[(y    ) * input.width + (x - 1)];
            mask[1][1] = input.pixels[(y    ) * input.width + (x    )];
            mask[1][2] = input.pixels[(y    ) * input.width + (x + 1)];
            mask[2][0] = input.pixels[(y + 1) * input.width + (x - 1)];
            mask[2][1] = input.pixels[(y + 1) * input.width + (x    )];
            mask[2][2] = input.pixels[(y + 1) * input.width + (x + 1)];

            for (png_uint_16 i = 0; i < 3; i++) {
                for (png_uint_16 j = 0; j < 3; j++) {
                    gx += mask[i][j] * GX[i][j];
                    gy += mask[i][j] * GY[i][j];
                }
            }

            double g = sqrt(gx * gx + gy * gy);
            g = (g < 0) ? 0 : (g > 255) ? 255 : g;
            
            output.pixels[y * input.width + x] = (unsigned char) g; 
        }
    }
}

void prepareAccumulator (const outputImage &img, uint32_t r_min, uint32_t r_max, uint32_t r_step, uint8_t pixel_threshold, double threshold, accumulator& votes) {
    votes.r_min = r_min;
    votes.r_max = r_max;
    votes.r_step = r_step;

    votes.pixel_threshold = pixel_threshold;
    votes.threshold = threshold;

    votes.width = img.width;
    votes.height = img.height;
    votes.nRadii = (r_max - r_min) / r_step + 1;

    votes.accum = new size_t [votes.width * votes.height * votes.nRadii]();
}

uint32_t CHT(const outputImage &input, accumulator& votes) {
    cout << "\nImplement CHT function in " << __FILE__ << " @ " << __LINE__; 
    // for each pixel that are above pixel_threshold, cast votes 

    // After voting, filter significant votes based on the criterion specified in the instruction of the exam
    // each significant vote becomes an identified circle
    
    // return the number of identified circles at this stage.
    return 0; 
}

uint32_t extractCircles(accumulator& votes, circles& circles) {
    cout << "\nImplement extractCircles function in " << __FILE__ << " @ " << __LINE__; 
    
    // the following line of code demonstrate how to put detected circles into the circles object
    // you shoudl remove it after you have extract circles from votes...
    circles.data.push_back({100, 150, 200, 1000}); // store a circle @ (100, 150) with radius 200 and 1000 votes.

    return 0; 
}

void freeAccumulator(accumulator& votes) {
    delete [] votes.accum;
}

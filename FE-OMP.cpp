#include "FE.hpp"
#include "PNGio.hpp"
#include <cmath>
#include <iostream>
using std::cout; 
// 
int readImage(const char *filename, inputImage &png) {
    cout << "\n[OMP]Check readImage in " << __FILE__ << " @ " << __LINE__; 
    // If you are okay with how image data is stored (defined in PNGio.cpp), you can simply use the following line of code for this function
    return pngRead(filename, png); 
}

void toGrayScale(const inputImage& input, outputImage &output) {
    cout << "\n[OMP]Modify toGrayScale function in " << __FILE__ << " @ " << __LINE__;
    // The following code only use the blue intensity as the output
    output.width = input.width; 
    output.height = input.height; 
    output.pixels = new png_byte[output.width * output.height]; 
    for(int y=0;y<input.height;++y) {
        for(int x=0;x<input.width;++x) {
            output.pixels[y*input.width+x] = input.row_pointers[y][4*x+2]; 
        }
    }
}

// https://en.wikipedia.org/wiki/Sobel_operator
void toEdge(const outputImage &input, outputImage &output) {
    cout << "\n[OMP]Modify toEdge function in " << __FILE__ << " @ " << __LINE__; 
    // The following code only copies input image to the output image
    output.width = input.width; 
    output.height = input.height; 
    output.pixels = new png_byte[output.width * output.height]; 
    for(png_uint_32 y=0;y<input.height;++y) {
        for(png_uint_32 x=0;x<input.width;++x) {
            output.pixels[y*input.width+x] = input.pixels[y*input.width+x]; 
        }
    }    
}

void prepareAccumulator (const outputImage &img, uint32_t r_min, uint32_t r_max, uint32_t r_step, uint8_t pixel_threshold, double threshold, accumulator& votes) {
    cout << "\n[OMP]Implement prepareAccumulator function in " << __FILE__ << " @ " << __LINE__; 
}

uint32_t CHT(const outputImage &input, accumulator& votes) {
    cout << "\n[OMP]Implement CHT function in " << __FILE__ << " @ " << __LINE__; 
    // for each pixel that are above pixel_threshold, cast votes 

    // After voting, filter significant votes based on the criterion specified in the instruction of the exam
    // each significant vote becomes an identified circle
    
    // return the number of identified circles at this stage.
    return 0; 
}

uint32_t extractCircles(accumulator& votes, circles& circles) {
    cout << "\n[OMP]Implement extractCircles function in " << __FILE__ << " @ " << __LINE__; 
    
    // the following line of code demonstrate how to put detected circles into the circles object
    // you shoudl remove it after you have extract circles from votes...
    circles.data.push_back({100, 150, 200, 1000}); // store a circle @ (100, 150) with radius 200 and 1000 votes.

    return 0; 
}

void freeAccumulator(accumulator& votes) {
    cout << "\n[OMP]Implement freeAccumulator function in " << __FILE__ << " @ " << __LINE__;
}

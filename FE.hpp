#pragma once 
#include "PNGio.hpp"

// Read PNG image file specified by the filename and store the image into inputImage struct called png.
int readImage(const char *filename, inputImage &png);

// Convert input color image with R, G, B, A components/channels into gray-scale image in the outputImage called output
void toGrayScale(const inputImage& input, outputImage &output); 

// Apply Sobel operator to turn gray-scale image in the input into edge-plot in the output
void toEdge(const outputImage &input, outputImage &output); 

// design your own accumulator struct to store parameters, voting results, and identified circles
struct accumulator {

}; 

// prepare the accumulator struct
void prepareAccumulator (const outputImage&, uint32_t r_min, uint32_t r_max, uint32_t r_step, uint8_t pixel_threshold, double threshold, accumulator& tobeprepared); 

// Do Circle Hough Transform to cast votes on detected pixels.  
// After voting, store detected circle parameters (based on significant votes) in the accumulator.
// This function returns the number of detected circles.
uint32_t CHT(const outputImage &input, accumulator& votes); 

// This function combines/merges detected circles based on certain criterion and stores the final result into circles struct
// This function returns the number of final/merged circles.
uint32_t extractCircles(accumulator& votes, circles& circles); 

// This function frees memory allocated in prepareAccumulator
void freeAccumulator(accumulator& votes);


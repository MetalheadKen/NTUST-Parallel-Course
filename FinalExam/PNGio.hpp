#pragma once
#include <png.h>

#include <cstdint>
#include <list>
#include <array>

using std::list; 
using std::array; 

// struct for PNG image
struct inputImage {
    int width, height;          // width & height of the image
    png_byte depth;             // color depth of the image
    png_byte color_type;        // type of the color data (e.g. PNG_COLOR_TYPE_GRAY, PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB_ALPHA ...)
	png_uint_32 stride; 		// bytes per row
    png_bytep *row_pointers; 	
}; 

// struct for output gray-scale images
struct outputImage {
    png_uint_32 width, height; 
    png_byte *pixels; // width * height
}; 

// struct for storing circles identified by HCT
struct circles {
    list<array<uint32_t, 4>> data;   // a list of arrays of size 4, each array contains center-x, center-y, radius, and # of votes
}; 

// This function reads a PNG image file specified by the filename and stores the image data into data
int pngRead(const char *filename, inputImage &data); 

// This function de-allocates memory allocated by pngRead or inputImage struct
void pngFree(inputImage &data); 

// This function saves gray-scale image data into PNG file specified by the filename.
int pngWrite(const char *filename, outputImage& data); 

// This function de-allocates memory for output image struct
void pngFree(outputImage &data); 

// This function allocates memory for grya-scale images stored in outputImage struct with specified width and height.  This function also blanks the canvas, i.e. fills the canvas with black color.
void blank(outputImage &canvas, png_uint_32 width, png_uint_32 height); 

// This function draws circles based on the data recorded in the circleData struct
void drawCircles(outputImage& canvas, circles& circleData, const png_byte v=0xff); 

#include <iostream>
using namespace std; 

#include "stopwatch.hpp"

#ifdef FEOMP
	#include "FE-OMP.hpp"
#else 
	#ifdef FEOCL
		#include "FE-OCL.hpp"
	#else 
		#include "FE.hpp"
	#endif
#endif


int main(int argc, char **argv) {
	if(argc<8) {
		cerr << "\n" << argv[0] << " [PNG filename] [OUTPUT filename] [r_min] [r_max] [r_step] [pixel threshold] [threshold]"; 
		return 255;
	}
	const uint32_t r_min = atoi(argv[3]); 
	const uint32_t r_max = atoi(argv[4]); 
	const uint32_t r_step = atoi(argv[5]); 
	const uint8_t pixel_threshold = atoi(argv[6]); 
	const double threshold = atof(argv[7]); 

	stopwatch timers[6]; 
	uint32_t counts[2]; 
	timers[0].start(); timers[1].start(); 

	//
	// 1. read PNG image file
	//
	inputImage sourceImage; 
	int status = readImage(argv[1], sourceImage); 
	if(status) {
		return status;
	}

	//
	// 2. Convert to gray-scale image 
	// 
	outputImage step1; 
	toGrayScale(sourceImage, step1); 
	timers[1].stop(); timers[2].start(); 

	//
	// 3. Apply Sobel filter
	//
	outputImage step2; 
	toEdge(step1, step2); 
	timers[2].stop(); timers[3].start(); 

	//
	// 4. prepare accumulator
	//
    accumulator votes; 
	prepareAccumulator(step2, r_min, r_max, r_step, pixel_threshold, threshold, votes); 
	timers[3].stop(); timers[4].start(); 

	//
	// 5. Do Hough circle Transform (HCT) to find circles
	//
	counts[0] = CHT(step2, votes); 
	timers[4].stop(); timers[5].start(); 

	//
	// 6. merge/combine similar circles and store them into circles
	//
	circles circles;
    counts[1] = extractCircles(votes, circles); 
	timers[5].stop();

	// 
	// 7. Output identified circles in the output image
	// 
	outputImage step3; 
	blank(step3, step2.width, step2.height); 
	drawCircles(step3, circles); 
	if( pngWrite(argv[2], step3) ) {
		cerr << "\nFailed to write PNG file: " << argv[2];
	}

	// pngWrite("gray.png", step1); 
	// pngWrite("sobol.png", step2); 
	
	pngFree(step3); 
	pngFree(step2); 
	pngFree(step1); 
	pngFree(sourceImage);

	freeAccumulator(votes); 

	timers[0].stop(); 

	cout << "\n###, "; 
	for(auto& t: timers) cout << t.elapsedTime() << ", "; 
	cout << counts[0] << ", " << counts[1] << endl; 
	
	return 0; 
}

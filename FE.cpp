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
    output.pixels = new png_byte [output.width * output.height]; 

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

static inline void _vote(accumulator& votes, int x, int y, int r, uint32_t& total) {
    if ((x > 0 && (uint32_t) x < votes.width) && (y > 0 && (uint32_t) y < votes.height)) {
        votes.accum[y * votes.width + x + r * votes.width * votes.height]++;
        total++;
    }
}

uint32_t CHT(const outputImage &input, accumulator& votes) {
    //cout << "\nImplement CHT function in " << __FILE__ << " @ " << __LINE__; 
    uint32_t total_votes = 0;
    
    // for each pixel that are above pixel_threshold, cast votes 
    for (png_uint_32 y = 0; y < input.height; y++) {
        for (png_uint_32 x = 0; x < input.width; x++) {
            uint32_t pos = y * input.width + x;

            if (input.pixels[pos] > votes.pixel_threshold) {
                for (uint32_t r = 0; r < votes.nRadii; r++) {
                    uint32_t radii = votes.r_min + r * votes.r_step;

                    _vote(votes, x, y - radii, r, total_votes);
                    _vote(votes, x, y + radii, r, total_votes);
                    _vote(votes, x - radii, y, r, total_votes);
                    _vote(votes, x + radii, y, r, total_votes);

                    for (uint32_t dy = 1; dy < 0.71 * radii; dy++) {
                        const uint32_t dx = radii * sqrt(1.0 - (double) dy * dy / (radii * radii));

                        _vote(votes, x - dx, y - dy, r, total_votes);
                        _vote(votes, x - dx, y + dy, r, total_votes);
                        _vote(votes, x + dx, y - dy, r, total_votes);
                        _vote(votes, x + dx, y + dy, r, total_votes);

                        _vote(votes, x - dy, y - dx, r, total_votes);
                        _vote(votes, x - dy, y + dx, r, total_votes);
                        _vote(votes, x + dy, y - dx, r, total_votes);
                        _vote(votes, x + dy, y + dx, r, total_votes);

                    }
                }
            }
        }
    }

    // After voting, filter significant votes based on the criterion specified in the instruction of the exam
    // each significant vote becomes an identified circle
    const double vote_threshold = votes.threshold * total_votes / (votes.width * votes.height * votes.nRadii);
    for (png_uint_32 y = 0; y < input.height; y++) {
        for (png_uint_32 x = 0; x < input.width; x++) {
            for (uint32_t r = 0; r < votes.nRadii; r++) {
                uint32_t pos = y * input.width + x + r * input.width * input.height;

                if (votes.accum[pos] > vote_threshold) {
                    struct circle circle;

                    circle.x     = x;
                    circle.y     = y;
                    circle.r     = votes.r_min + r * votes.r_step;
                    circle.votes = votes.accum[pos];

                    votes.circles.push_back(circle);
                }
            }
        }
    }
    
    // return the number of identified circles at this stage.
    return votes.circles.size(); 
}

uint32_t extractCircles(accumulator& votes, circles& circles) {
    //cout << "\nImplement extractCircles function in " << __FILE__ << " @ " << __LINE__; 
    const int diff_dist = 20;
    bool *is_merge = new bool [votes.circles.size()];

    // Find the similar circle and only pick one of them
    for (uint32_t i = 0; i < votes.circles.size(); i++) {
        if (is_merge[i]) continue;

        is_merge[i] = true;
        auto c1 = votes.circles[i];
        std::vector<circle> similar;

        for (uint32_t j = 0; j < votes.circles.size(); j++) {
            if (is_merge[j]) continue;

            auto c2 = votes.circles[j];

            if (abs(c1.x - c2.x) < diff_dist && 
                abs(c1.y - c2.y) < diff_dist &&
                abs(c1.r - c2.r) < diff_dist) {
                is_merge[j] = true;
                similar.push_back(c2);
            }
        }

        int size = similar.size();
        struct circle sum = { 0 };
        for (const auto& c : similar) {
            sum.x     += c.x / size;
            sum.y     += c.y / size;
            sum.r     += c.r / size;
            sum.votes += c.votes / size;
        }

        circles.data.push_back({sum.x, sum.y, sum.r, sum.votes});
    }
    
    // the following line of code demonstrate how to put detected circles into the circles object
    // you shoudl remove it after you have extract circles from votes...
    // circles.data.push_back({100, 150, 200, 1000}); // store a circle @ (100, 150) with radius 200 and 1000 votes.

    return circles.data.size(); 
}

void freeAccumulator(accumulator& votes) {
    delete [] votes.accum;
}

#include "FE-OCL.hpp"
#include "PNGio.hpp"
#include "YoUtil.hpp"
#include <cmath>
#include <iostream>
using std::cout; 
// 
int readImage(const char *filename, inputImage &png) {
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

uint32_t CHT(const outputImage &input, accumulator& votes) {
    uint32_t total_votes = 0;

    cl::CommandQueue cmdQueue;
    cl::Context context;
    cl::Device device;
    cl::Program *program;

    // Prepare GPU
    //try {
        YoUtil::GPU gpu(CL_DEVICE_TYPE_GPU, false, "");
        gpu.addProgramByFile("FE-OCL", "FE-OCL.cl");
        cmdQueue = gpu.getCommandQueue(0);
        context = cmdQueue.getInfo<CL_QUEUE_CONTEXT>();
        device = cmdQueue.getInfo<CL_QUEUE_DEVICE>();
        program = gpu.getProgram("FE-OCL"); 
    //} catch (cl::Error &e) {
        //std::cerr << "CHT: " << e.what() << std::endl;
        //std::cerr << "Error no: " << e.err() << std::endl;
    //}

    cl_uint in_bytes = sizeof(input.pixels[0]) * input.width * input.height;
    cl_uint out_bytes = sizeof(votes.accum[0]) * votes.width * votes.height * votes.nRadii;

    cl_mem_flags in_flag = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
    cl_mem_flags out_flag = CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR;

    cl::Buffer in_buf(context, in_flag, in_bytes, input.pixels);
    cl::Buffer out_buf(context, out_flag, out_bytes, votes.accum);
    cl::Buffer cnt_buf(context, out_flag, sizeof(total_votes), &total_votes);

    // Setting kernel arguments
    cl::Kernel kernel(*program, "CHT_vote");
    kernel.setArg(0, (cl_uint) votes.width);
    kernel.setArg(1, (cl_uint) votes.height);
    kernel.setArg(2, (cl_uint) votes.nRadii);
    kernel.setArg(3, (cl_uint) votes.r_min);
    kernel.setArg(4, (cl_uint) votes.r_step);
    kernel.setArg(5, (cl_uint) votes.pixel_threshold);
    kernel.setArg(6, in_buf);
    kernel.setArg(7, out_buf);
    kernel.setArg(8, cnt_buf);

    // Setting work group
    cl_uint local_size = 32;
    cl_uint global_size_x = (votes.width + local_size - 1) / local_size * local_size;
    cl_uint global_size_y = (votes.height + local_size - 1) / local_size * local_size;
    cl::NDRange local(local_size, local_size);
    cl::NDRange global(global_size_x, global_size_y);

    // Start running kernel
    cmdQueue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
    cmdQueue.enqueueReadBuffer(out_buf, CL_TRUE, 0, out_bytes, votes.accum);
    cmdQueue.enqueueReadBuffer(cnt_buf, CL_TRUE, 0, sizeof(total_votes), &total_votes);

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
    const int diff_dist = 20;
    bool *is_merge = new bool [votes.circles.size()]();

    // Find the similar circle and only pick one of them
    for (uint32_t i = 0; i < votes.circles.size(); i++) {
        if (is_merge[i]) continue;

        is_merge[i] = true;
        auto c1 = votes.circles[i];

        std::vector<circle> similar;
        similar.push_back(c1);

        for (uint32_t j = 0; j < votes.circles.size(); j++) {
            if (is_merge[j]) continue;

            auto c2 = votes.circles[j];

            if (abs((int32_t) c1.x - (int32_t) c2.x) < diff_dist && 
                abs((int32_t) c1.y - (int32_t) c2.y) < diff_dist &&
                abs((int32_t) c1.r - (int32_t) c2.r) < diff_dist) {
                is_merge[j] = true;
                similar.push_back(c2);
            }
        }

        // Average all circles parameters
        struct circle sum = { 0 };
        for (const auto& c : similar) {
            sum.x     += c.x;
            sum.y     += c.y;
            sum.r     += c.r;
            sum.votes += c.votes;
        }

        int size = similar.size();
        sum.x     /= size;
        sum.y     /= size;
        sum.r     /= size;
        sum.votes /= size;

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

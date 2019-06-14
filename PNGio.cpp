#include <cstdio>
#include <cmath>
#include "PNGio.hpp"

using std::nothrow;

// https://gist.github.com/niw/5963798
// http://zarb.org/~gc/html/libpng.html
// http://www.libpng.org/pub/png/libpng-manual.txt
// http://www.libpng.org/pub/png/book/chapter13.html
int pngRead(const char *filename, inputImage &data) {
    FILE *fp = fopen(filename, "rb"); 
    if(!fp) return 1;
    
    // create PNG read struct
    png_structp png_read = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png_read) {
        fclose(fp);
        return 2;
    }

    // create PNG png_info struct
    png_infop png_info = png_create_info_struct(png_read);
    if(!png_info) {
        png_destroy_read_struct(&png_read, NULL, NULL); 
        fclose(fp); 
        return 3;
    }

    // something to do with error handling
    if(setjmp(png_jmpbuf(png_read))) {
        png_destroy_read_struct(&png_read, &png_info, NULL); 
        fclose(fp); 
        return 4;
    }

    png_init_io(png_read, fp);
    png_read_info(png_read, png_info);

    data.width      = png_get_image_width(png_read, png_info);
    data.height     = png_get_image_height(png_read, png_info);
    data.color_type = png_get_color_type(png_read, png_info);
    data.depth  = png_get_bit_depth(png_read, png_info);
    
    if(data.depth == 16) png_set_strip_16(png_read); 
    if(data.color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png_read);
    if(data.color_type == PNG_COLOR_TYPE_GRAY && data.depth < 8) png_set_expand_gray_1_2_4_to_8(png_read);

    if( png_get_valid(png_read, png_info, PNG_INFO_tRNS) ) png_set_tRNS_to_alpha(png_read);

    // These color_type don't have an alpha channel then fill it with 0xff.
    if(data.color_type == PNG_COLOR_TYPE_RGB || data.color_type == PNG_COLOR_TYPE_GRAY || data.color_type == PNG_COLOR_TYPE_PALETTE) png_set_filler(png_read, 0xFF, PNG_FILLER_AFTER);

    // Expand gray-scale images to color images
    if(data.color_type == PNG_COLOR_TYPE_GRAY || data.color_type == PNG_COLOR_TYPE_GRAY_ALPHA) png_set_gray_to_rgb(png_read);
    // 
    // if (data.color_type == PNG_COLOR_TYPE_RGB || data.color_type == PNG_COLOR_TYPE_RGB_ALPHA) png_set_rgb_to_gray(png_read, 1, -1, -1); 

    // if( data.color_type == PNG_COLOR_TYPE_RGB_ALPHA || data.color_type == PNG_COLOR_TYPE_GRAY_ALPHA) 
    png_set_strip_alpha(png_read); 

    png_read_update_info(png_read, png_info);

	// now read the actual PNG file image data
	// data.row_pointers = (png_bytep*) malloc (sizeof(png_bytep) * height);
	data.stride=png_get_rowbytes(png_read, png_info);
	data.row_pointers = new (nothrow) png_bytep[ data.height ]; 
	if(! data.row_pointers ) {
		png_destroy_read_struct(&png_read, &png_info, NULL);
		fclose(fp); 
		return 5;
	}
	data.row_pointers[0] = new (nothrow) png_byte[ data.height * data.stride ]; 
	if(! data.row_pointers[0] ) {
		delete[] data.row_pointers;
		png_destroy_read_struct(&png_read, &png_info, NULL);
		fclose(fp); 
		return 6;
	}
	for(int y = 1; y < data.height; y++) {
		data.row_pointers[y] = data.row_pointers[0] + y * data.stride; 
	}

	png_read_image( png_read, data.row_pointers ); 
	png_read_end( png_read, NULL );

	// Finally, free up used memory
	png_destroy_read_struct(&png_read, &png_info, NULL);

    fclose(fp); 

    return 0;
}

void pngFree(inputImage &data) {
	delete[] data.row_pointers[0]; 
	delete[] data.row_pointers; 
}

void pngFree(outputImage &data) {
    delete[] data.pixels; 
}

int pngWrite(const char *filename, outputImage &data) {
    FILE *outp = fopen(filename, "wb");
    if(!outp) return 1; 

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fclose(outp); 
        return 2;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        png_destroy_write_struct(&png, NULL); 
        fclose(outp); 
        return 3;
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct(&png, &info); 
        fclose(outp); 
        return 4;
    }

    png_init_io(png, outp);

    const int depth=8;
    png_set_IHDR(
        png,
        info,
        data.width, data.height,
        depth,
        PNG_COLOR_TYPE_GRAY,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);

    // To remove the alpha channel for PNG_COLOR_TYPE_RGB format,
    // Use png_set_filler().
    //png_set_filler(png, 0, PNG_FILLER_AFTER);
    png_bytep *row_pointers = new (nothrow) png_bytep[ data.height ]; 
    if(!row_pointers) {
        png_destroy_write_struct(&png, &info); 
        fclose(outp); 
        return 5;
    }
    for(png_uint_32 i=0;i<data.height;++i) {
        row_pointers[i] = data.pixels + i * data.width;
    }
    png_write_image(png, row_pointers);
    png_write_end(png, NULL);

    png_destroy_write_struct(&png, &info); 
    delete[] row_pointers;
    fclose(outp);
    return 0;
}

void blank(outputImage &canvas, png_uint_32 width, png_uint_32 height) {
    canvas.width = width; 
    canvas.height = height; 
    canvas.pixels = new png_byte[width * height]; 
    for(size_t i=0;i<width*height;++i) canvas.pixels[i] = 0; 
}

void drawPixel(outputImage&canvas, uint32_t x, uint32_t y, const png_byte v) {
    if(x >= canvas.width || x < 0 || y >= canvas.height || y<0 ) return;
    auto data = canvas.pixels + y * canvas.width + x;
    *data = v;
}

void drawCircles(outputImage& canvas, circles& circleData, const png_byte v) {
    for(auto &i: circleData.data) {
        uint32_t x = i[0]; 
        uint32_t y = i[1];
        uint32_t r = i[2]; 
        drawPixel(canvas, x+r, y, v); 
        drawPixel(canvas, x-r, y, v); 
        drawPixel(canvas, x, y+r, v); 
        drawPixel(canvas, x, y-r, v); 

        for(int dy = 1; dy < 0.71*r; dy++) {
            auto dx = uint32_t(r * sqrt(1.0 - double(dy)*dy/r/r)); 
            drawPixel(canvas, x+dx, y-dy, v); 
            drawPixel(canvas, x-dx, y-dy, v); 
            drawPixel(canvas, x+dx, y+dy, v); 
            drawPixel(canvas, x-dx, y+dy, v); 
            drawPixel(canvas, x+dy, y-dx, v); 
            drawPixel(canvas, x-dy, y-dx, v); 
            drawPixel(canvas, x+dy, y+dx, v); 
            drawPixel(canvas, x-dy, y+dx, v); 
        }
    }
}

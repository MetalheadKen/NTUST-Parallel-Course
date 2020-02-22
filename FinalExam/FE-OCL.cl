//#pragma OPENCL EXTENSION cl_khr_int64_base_atomics:enable
//#pragma OPENCL EXTENSION cl_khr_fp:enable

// Put your OpenCL kernel code in this file.
void _vote(
    __global uint *accum, 
    uint x, uint y, uint r, uint width, uint height, __global uint *cnt) {
    if ((x > 0 && x < width) && (y > 0 && y < height)) {
        atomic_inc(accum + (y * width + x + r * width * height));
        atomic_inc(cnt);
    }
}

__kernel
void CHT_vote(
    uint width, uint height, uint nRadii,
    uint r_min, uint r_step, uint pixel_threshold,
    __global uchar *pixels, __global uint *accum, __global uint *cnt) {
    uint x = get_global_id(0);
    uint y = get_global_id(1);

    // for each pixel that are above pixel_threshold, cast votes 
    if ((x > 0 && x < width) && (y > 0 && y < height)) {
        uint pos = y * width + x;

        if (pixels[pos] > pixel_threshold) {
            for (uint r = 0; r < nRadii; r++) {
                uint radii = r_min + r * r_step;

                _vote(accum, x, y - radii, r, width, height, cnt);
                _vote(accum, x, y + radii, r, width, height, cnt);
                _vote(accum, x - radii, y, r, width, height, cnt);
                _vote(accum, x + radii, y, r, width, height, cnt);

                for (uint dy = 1; dy < 0.71 * radii; dy += 8) {
                    const uint dx = radii * sqrt(1.0 - (double) dy * dy / (radii * radii));

                    _vote(accum, x - dx, y - dy, r, width, height, cnt);
                    _vote(accum, x - dx, y + dy, r, width, height, cnt);
                    _vote(accum, x + dx, y - dy, r, width, height, cnt);
                    _vote(accum, x + dx, y + dy, r, width, height, cnt);

                    _vote(accum, x - dy, y - dx, r, width, height, cnt);
                    _vote(accum, x - dy, y + dx, r, width, height, cnt);
                    _vote(accum, x + dy, y - dx, r, width, height, cnt);
                    _vote(accum, x + dy, y + dx, r, width, height, cnt);

                }
            }
        }
    }
}


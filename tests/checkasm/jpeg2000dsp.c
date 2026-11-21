/*
 * Copyright (c) 2015 James Almer
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with FFmpeg; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "checkasm.h"
#include "libavcodec/jpeg2000dsp.h"
#include "libavutil/common.h"
#include "libavutil/internal.h"
#include "libavutil/intreadwrite.h"
#include "libavutil/mem_internal.h"

#define BUF_SIZE 512

#define randomize_buffers()                 \
    do {                                    \
        int i;                              \
        for (i = 0; i < BUF_SIZE*3; i++)    \
            src[i] = rnd();                 \
    } while (0)

#define randomize_buffers_float()                    \
    do {                                             \
        int i;                                       \
        for (i = 0; i < BUF_SIZE*3; i++)             \
            src[i] = (float)rnd() / (UINT_MAX >> 5); \
    } while (0)

static const int mct_sizes[] = {
    0, 1, 2, 3, 4, 5, 7, 8, 15, 16, 17, 31, 32, 33,
    63, 64, 65, 127, 128, 129, 255, 256, 257, 511, BUF_SIZE
};

static void check_rct_int(void)
{
    LOCAL_ALIGNED_32(int32_t, src, [BUF_SIZE*3]);
    LOCAL_ALIGNED_32(int32_t, ref, [BUF_SIZE*3]);
    LOCAL_ALIGNED_32(int32_t, new, [BUF_SIZE*3]);
    int32_t *ref0 = &ref[BUF_SIZE*0], *new0 = &new[BUF_SIZE*0];
    int32_t *ref1 = &ref[BUF_SIZE*1], *new1 = &new[BUF_SIZE*1];
    int32_t *ref2 = &ref[BUF_SIZE*2], *new2 = &new[BUF_SIZE*2];
    size_t n;

    declare_func(void, void *src0, void *src1, void *src2, int csize);

    for (n = 0; n < FF_ARRAY_ELEMS(mct_sizes); n++) {
        int csize = mct_sizes[n];

        randomize_buffers();
        memcpy(ref, src, BUF_SIZE * 3 * sizeof(*src));
        memcpy(new, src, BUF_SIZE * 3 * sizeof(*src));
        call_ref(ref0, ref1, ref2, csize);
        call_new(new0, new1, new2, csize);
        /* Compare the csize samples the transform returns. SIMD that
         * rounds up to a full vector may touch the caller's padding,
         * which the decoder allocates. */
        if (memcmp(ref0, new0, csize * sizeof(*src)) ||
            memcmp(ref1, new1, csize * sizeof(*src)) ||
            memcmp(ref2, new2, csize * sizeof(*src))) {
            fprintf(stderr, "csize=%d\n", csize);
            fail();
        }
    }
    randomize_buffers();
    memcpy(new, src, BUF_SIZE * 3 * sizeof(*src));
    bench_new(new0, new1, new2, BUF_SIZE);
}

static void check_ict_float(void)
{
    LOCAL_ALIGNED_32(float, src, [BUF_SIZE*3]);
    LOCAL_ALIGNED_32(float, ref, [BUF_SIZE*3]);
    LOCAL_ALIGNED_32(float, new, [BUF_SIZE*3]);
    float *ref0 = &ref[BUF_SIZE*0], *new0 = &new[BUF_SIZE*0];
    float *ref1 = &ref[BUF_SIZE*1], *new1 = &new[BUF_SIZE*1];
    float *ref2 = &ref[BUF_SIZE*2], *new2 = &new[BUF_SIZE*2];

    size_t n;

    declare_func(void, void *src0, void *src1, void *src2, int csize);

    for (n = 0; n < FF_ARRAY_ELEMS(mct_sizes); n++) {
        int csize = mct_sizes[n];

        randomize_buffers_float();
        memcpy(ref, src, BUF_SIZE * 3 * sizeof(*src));
        memcpy(new, src, BUF_SIZE * 3 * sizeof(*src));
        call_ref(ref0, ref1, ref2, csize);
        call_new(new0, new1, new2, csize);
        /* Compare the csize samples the transform returns. SIMD that
         * rounds up to a full vector may touch the caller's padding,
         * which the decoder allocates. */
        if (!float_near_abs_eps_array(ref0, new0, 1.0e-5, csize) ||
            !float_near_abs_eps_array(ref1, new1, 1.0e-5, csize) ||
            !float_near_abs_eps_array(ref2, new2, 1.0e-5, csize)) {
            fprintf(stderr, "csize=%d\n", csize);
            fail();
        }
    }
    randomize_buffers_float();
    memcpy(new, src, BUF_SIZE * 3 * sizeof(*src));
    bench_new(new0, new1, new2, BUF_SIZE);
}

void checkasm_check_jpeg2000dsp(void)
{
    Jpeg2000DSPContext h;

    ff_jpeg2000dsp_init(&h);

    if (check_func(h.mct_decode[FF_DWT53], "jpeg2000_rct_int"))
        check_rct_int();
    if (check_func(h.mct_decode[FF_DWT97], "jpeg2000_ict_float"))
        check_ict_float();

    report("mct_decode");
}

/*
 * Discrete wavelet transform
 * Copyright (c) 2007 Kamil Nowosad
 * Copyright (c) 2013 Nicolas Bertrand <nicoinattendu@gmail.com>
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include <inttypes.h>
#include <stdio.h>

#include "libavutil/common.h"
#include "libavutil/intfloat.h"

#include "libavcodec/jpeg2000dwt.c"

#include "libavutil/lfg.h"

#define MAX_W 256

static int test_dwt(int *array, int *ref, int border[2][2], int decomp_levels, int type, int max_diff) {
    int ret, j;
    DWTContext s1={{{0}}}, *s= &s1;
    int64_t err2 = 0;

    ret = ff_jpeg2000_dwt_init(s,  border, decomp_levels, type);
    if (ret < 0) {
        fprintf(stderr, "ff_jpeg2000_dwt_init failed\n");
        return 1;
    }
    ret = ff_dwt_encode(s, array);
    if (ret < 0) {
        fprintf(stderr, "ff_dwt_encode failed\n");
        return 1;
    }
    if (type == FF_DWT97_INT) {
        // pre-scaling to simulate dequantization which places the binary point at 1 bit above from LSB
        for (j = 0; j< s->linelen[decomp_levels-1][0] * s->linelen[decomp_levels-1][1]; j++)
            array[j] = (uint32_t)array[j] << I_PRESHIFT;
    }
    ret = ff_dwt_decode(s, array);
    if (ret < 0) {
        fprintf(stderr, "ff_dwt_encode failed\n");
        return 1;
    }
    for (j = 0; j<MAX_W * MAX_W; j++) {
        if (FFABS(array[j] - (int64_t)ref[j]) > max_diff) {
            fprintf(stderr, "mismatch at %d (%d != %d) decomp:%d border %d %d %d %d\n",
                    j, array[j], ref[j],decomp_levels, border[0][0], border[0][1], border[1][0], border[1][1]);
            return 2;
        }
        err2 += (array[j] - ref[j]) * (int64_t)(array[j] - ref[j]);
        array[j] = ref[j];
    }
    ff_dwt_destroy(s);

    printf("%s, decomp:%2d border %3d %3d %3d %3d milli-err2:%9"PRId64"\n",
           type == FF_DWT53 ? "5/3i" : "9/7i",
           decomp_levels, border[0][0], border[0][1], border[1][0], border[1][1],
           1000*err2 / ((border[0][1] - border[0][0])*(border[1][1] - border[1][0])));

    return 0;
}

static int test_dwtf(float *array, float *ref, int border[2][2], int decomp_levels, float max_diff) {
    int ret, j;
    DWTContext s1={{{0}}}, *s= &s1;
    double err2 = 0;

    ret = ff_jpeg2000_dwt_init(s,  border, decomp_levels, FF_DWT97);
    if (ret < 0) {
        fprintf(stderr, "ff_jpeg2000_dwt_init failed\n");
        return 1;
    }
    ret = ff_dwt_encode(s, array);
    if (ret < 0) {
        fprintf(stderr, "ff_dwt_encode failed\n");
        return 1;
    }
    ret = ff_dwt_decode(s, array);
    if (ret < 0) {
        fprintf(stderr, "ff_dwt_encode failed\n");
        return 1;
    }
    for (j = 0; j<MAX_W * MAX_W; j++) {
        if (FFABS(array[j] - ref[j]) > max_diff) {
            fprintf(stderr, "mismatch at %d (%f != %f) decomp:%d border %d %d %d %d\n",
                    j, array[j], ref[j],decomp_levels, border[0][0], border[0][1], border[1][0], border[1][1]);
            return 2;
        }
        err2 += (array[j] - ref[j]) * (array[j] - ref[j]);
        array[j] = ref[j];
    }
    ff_dwt_destroy(s);

    printf("9/7f, decomp:%2d border %3d %3d %3d %3d err2:%20.3f\n",
           decomp_levels, border[0][0], border[0][1], border[1][0], border[1][1],
           err2 / ((border[0][1] - border[0][0])*(border[1][1] - border[1][0])));

    return 0;
}

static int array[MAX_W * MAX_W];
static int ref  [MAX_W * MAX_W];
static float arrayf[MAX_W * MAX_W];
static float reff  [MAX_W * MAX_W];
static float arrayf_dec[MAX_W * MAX_W];
static float arrayf_dec_ref[MAX_W * MAX_W];

/* Pre-fusion boundary extension from base 17de9f8597. Kept local so the
 * reference does not share code with the production path. */
static void extend97_float_ref(float *p, int i0, int i1)
{
    int i;

    for (i = 1; i <= 4; i++) {
        p[i0 - i]     = p[i0 + i];
        p[i1 + i - 1] = p[i1 - i - 1];
    }
}

/* Pre-fusion 1-D 9/7 float lifting loops from base 17de9f8597. */
static void sr_1d97_float_ref(float *p, int i0, int i1)
{
    int i;

    if (i1 <= i0 + 1) {
        if (i0 == 1)
            p[1] *= F_LFTG_K/2;
        else
            p[0] *= F_LFTG_X;
        return;
    }

    extend97_float_ref(p, i0, i1);

    for (i = (i0 >> 1) - 1; i < (i1 >> 1) + 2; i++)
        p[2 * i]     -= F_LFTG_DELTA * (p[2 * i - 1] + p[2 * i + 1]);
    /* step 4 */
    for (i = (i0 >> 1) - 1; i < (i1 >> 1) + 1; i++)
        p[2 * i + 1] -= F_LFTG_GAMMA * (p[2 * i]     + p[2 * i + 2]);
    /*step 5*/
    for (i = (i0 >> 1); i < (i1 >> 1) + 1; i++)
        p[2 * i]     += F_LFTG_BETA  * (p[2 * i - 1] + p[2 * i + 1]);
    /* step 6 */
    for (i = (i0 >> 1); i < (i1 >> 1); i++)
        p[2 * i + 1] += F_LFTG_ALPHA * (p[2 * i]     + p[2 * i + 2]);
}

static int check_sr_1d97_float(int i0, int i1, AVLFG *prng)
{
    float buf_a[MAX_W + 16], buf_b[MAX_W + 16];
    size_t k;

    for (k = 0; k < FF_ARRAY_ELEMS(buf_a); k++) {
        float v = (av_lfg_get(prng) & 0xFFFFFF) / 4096.0f - 2048.0f;
        buf_a[k] = buf_b[k] = v;
    }

    sr_1d97_float(buf_a + 5, i0, i1);
    sr_1d97_float_ref(buf_b + 5, i0, i1);

    for (k = 0; k < FF_ARRAY_ELEMS(buf_a); k++) {
        uint32_t ua = av_float2int(buf_a[k]);
        uint32_t ub = av_float2int(buf_b[k]);
        if (ua != ub) {
            fprintf(stderr,
                    "1D 9/7f mismatch: i0=%d i1=%d at idx %d: "
                    "%g (0x%08"PRIx32") != %g (0x%08"PRIx32")\n",
                    i0, i1, (int)k - 5, buf_a[k], ua, buf_b[k], ub);
            return 1;
        }
    }

    return 0;
}

static int test_sr_1d97_float(void)
{
    static const struct { int i0, i1; } edge_cases[] = {
        { 0, 0 }, { 1, 1 }, { 2, 2 },
        { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 4 },
    };
    AVLFG prng;
    size_t idx;
    int i0, len, trial, ret;

    av_lfg_init(&prng, 3);

    for (idx = 0; idx < FF_ARRAY_ELEMS(edge_cases); idx++) {
        for (trial = 0; trial < 8; trial++) {
            ret = check_sr_1d97_float(edge_cases[idx].i0, edge_cases[idx].i1,
                                      &prng);
            if (ret)
                return ret;
        }
    }

    /* Exhaustive short lengths covering overlapping boundary extension,
     * odd/even lengths, and both parities (i0 = 0, 1). */
    for (len = 2; len <= 35; len++) {
        for (i0 = 0; i0 <= 1; i0++) {
            ret = check_sr_1d97_float(i0, i0 + len, &prng);
            if (ret)
                return ret;
        }
    }

    /* Deterministic randomized sweep over (i0, i1) ranges, including
     * length 1 (i1 <= i0 + 1), odd/even lengths, and full MAX_W spans. */
    for (trial = 0; trial < 512; trial++) {
        i0 = av_lfg_get(&prng) % 4;
        len = 1 + (int)(av_lfg_get(&prng) % MAX_W);
        ret = check_sr_1d97_float(i0, i0 + len, &prng);
        if (ret)
            return ret;
    }

    return 0;
}

/*
 * Pre-patch 9/7 float synthesis schedule using the unfused four-pass lifting
 * loops (sr_1d97_float_ref) and the original per-row and per-column
 * gather/scatter order so both the fused 1-D kernel and the blocked vertical
 * pass are checked against it.
 */
static void dwt_decode97_float_ref(DWTContext *s, float *t)
{
    int lev;
    int w       = s->linelen[s->ndeclevels - 1][0];
    float *line = s->f_linebuf;
    float *data = t;
    /* position at index O of line range [0-5,w+5] cf. extend function */
    line += 5;

    for (lev = 0; lev < s->ndeclevels; lev++) {
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        float *l;
        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++) {
            int i, j = 0;
            // copy with interleaving
            for (i = mh; i < lh; i += 2, j++)
                l[i] = data[w * lp + j];
            for (i = 1 - mh; i < lh; i += 2, j++)
                l[i] = data[w * lp + j];

            sr_1d97_float_ref(line, mh, mh + lh);

            for (i = 0; i < lh; i++)
                data[w * lp + i] = l[i];
        }

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;
            // copy with interleaving
            for (i = mv; i < lv; i += 2, j++)
                l[i] = data[w * j + lp];
            for (i = 1 - mv; i < lv; i += 2, j++)
                l[i] = data[w * j + lp];

            sr_1d97_float_ref(line, mv, mv + lv);

            for (i = 0; i < lv; i++)
                data[w * i + lp] = l[i];
        }
    }
}

/* The production decoder and the unfused scalar reference perform identical
 * arithmetic in the same evaluation order, so outputs are compared bitwise
 * with exact equality (no tolerance). */
static int test_dwt97_decode(int border[2][2], int decomp_levels, AVLFG *prng)
{
    DWTContext s = { 0 }, ref = { 0 };
    int w, h, ret;
    int first = -1;
    float first_a = 0, first_b = 0, max_diff = 0;
    uint32_t first_ua = 0, first_ub = 0;
    size_t i, n;

    ret = ff_jpeg2000_dwt_init(&s, border, decomp_levels, FF_DWT97);
    if (ret < 0) {
        fprintf(stderr, "ff_jpeg2000_dwt_init failed\n");
        return 1;
    }
    ret = ff_jpeg2000_dwt_init(&ref, border, decomp_levels, FF_DWT97);
    if (ret < 0) {
        fprintf(stderr, "ff_jpeg2000_dwt_init (reference) failed\n");
        ff_dwt_destroy(&s);
        return 1;
    }

    w = s.linelen[decomp_levels - 1][0];
    h = s.linelen[decomp_levels - 1][1];
    n = (size_t)w * h;
    if (n > MAX_W * MAX_W) {
        fprintf(stderr, "test buffer too small: %dx%d\n", w, h);
        ff_dwt_destroy(&s);
        ff_dwt_destroy(&ref);
        return 1;
    }

    for (i = 0; i < n; i++) {
        float v = (av_lfg_get(prng) % 2000001) / 1000000.0f - 1.0f;
        arrayf_dec[i] = arrayf_dec_ref[i] = v;
    }

    ret = ff_dwt_decode(&s, arrayf_dec);
    if (ret < 0) {
        fprintf(stderr, "ff_dwt_decode failed\n");
        ff_dwt_destroy(&s);
        ff_dwt_destroy(&ref);
        return 1;
    }
    dwt_decode97_float_ref(&ref, arrayf_dec_ref);

    for (i = 0; i < n; i++) {
        uint32_t ua = av_float2int(arrayf_dec[i]);
        uint32_t ub = av_float2int(arrayf_dec_ref[i]);
        float diff = FFABS(arrayf_dec[i] - arrayf_dec_ref[i]);
        max_diff = FFMAX(max_diff, diff);
        if (ua != ub && first < 0) {
            first = (int)i;
            first_a = arrayf_dec[i];
            first_b = arrayf_dec_ref[i];
            first_ua = ua;
            first_ub = ub;
        }
    }

    if (first >= 0) {
        fprintf(stderr,
                "9/7f decode mismatch: %dx%d border %d %d %d %d levels %d "
                "at %d: %g (0x%08"PRIx32") != %g (0x%08"PRIx32") (max diff %g)\n",
                w, h, border[0][0], border[0][1], border[1][0], border[1][1],
                decomp_levels, first, first_a, first_ua, first_b, first_ub,
                max_diff);
        ret = 1;
    } else {
        ret = 0;
    }

    ff_dwt_destroy(&s);
    ff_dwt_destroy(&ref);
    return ret;
}

static int test_dwt97_decode_ref(void)
{
    /* Widths below F_DWT_BLOCK use the scalar fallback, 16..31 cover every
     * block tail once, 32..47 add two-block tails, and the larger widths
     * exercise multi-block tails. Heights below 8 exercise the scalar
     * vertical path. The fixed matrix is 4,288 cases, followed by a
     * 512-case seeded randomized sweep over widths/heights 1..MAX_W. */
    static const int ws[] = { 8, 9, 12, 15,
                              16, 17, 18, 19, 20, 21, 22, 23,
                              24, 25, 26, 27, 28, 29, 30, 31,
                              32, 33, 34, 40, 47, 48, 63, 64, 100, 128 };
    static const int hs[] = { 7, 8, 9, 12, 16, 17, 54, 55, 90, 100 };
    AVLFG prng;
    size_t wi, hi;
    int lev, x0, y0, i, ret;
    int border[2][2];

    av_lfg_init(&prng, 2);

    for (wi = 0; wi < FF_ARRAY_ELEMS(ws); wi++) {
        for (hi = 0; hi < FF_ARRAY_ELEMS(hs); hi++) {
            int w = ws[wi], h = hs[hi];
            int maxlev = 0, t = FFMIN(w, h);

            while ((1 << (maxlev + 1)) <= t && maxlev < 6)
                maxlev++;

            for (lev = 1; lev <= maxlev; lev++) {
                for (x0 = 0; x0 <= 1; x0++) {
                    for (y0 = 0; y0 <= 1; y0++) {
                        border[0][0] = x0;
                        border[0][1] = x0 + w;
                        border[1][0] = y0;
                        border[1][1] = y0 + h;

                        ret = test_dwt97_decode(border, lev, &prng);
                        if (ret)
                            return ret;
                    }
                }
            }
        }
    }

    for (i = 0; i < 512; i++) {
        int w = 1 + (int)(av_lfg_get(&prng) % MAX_W);
        int h = 1 + (int)(av_lfg_get(&prng) % MAX_W);
        int maxlev = 0, t = FFMIN(w, h);

        while ((1 << (maxlev + 1)) <= t && maxlev < 6)
            maxlev++;
        if (maxlev < 1)
            maxlev = 1;

        lev = 1 + (int)(av_lfg_get(&prng) % maxlev);
        x0  = av_lfg_get(&prng) & 1;
        y0  = av_lfg_get(&prng) & 1;

        border[0][0] = x0;
        border[0][1] = x0 + w;
        border[1][0] = y0;
        border[1][1] = y0 + h;

        ret = test_dwt97_decode(border, lev, &prng);
        if (ret)
            return ret;
    }

    return 0;
}

int main(void) {
    AVLFG prng;
    int i,j;
    int border[2][2];
    int ret, decomp_levels;

    av_lfg_init(&prng, 1);

    ret = test_sr_1d97_float();
    if (ret)
        return ret;

    ret = test_dwt97_decode_ref();
    if (ret)
        return ret;

    for (i = 0; i<MAX_W * MAX_W; i++)
        arrayf[i] = reff[i] = array[i] = ref[i] =  av_lfg_get(&prng) % 2048;

    for (i = 0; i < 100; i++) {
        for (j=0; j<4; j++)
            border[j>>1][j&1] = av_lfg_get(&prng) % MAX_W;
        if (border[0][0] >= border[0][1] || border[1][0] >= border[1][1])
            continue;
        decomp_levels = av_lfg_get(&prng) % FF_DWT_MAX_DECLVLS;

        ret = test_dwt(array, ref, border, decomp_levels, FF_DWT53, 0);
        if (ret)
            return ret;
        ret = test_dwt(array, ref, border, decomp_levels, FF_DWT97_INT, FFMIN(7+5*decomp_levels, 15+3*decomp_levels));
        if (ret)
            return ret;
        ret = test_dwtf(arrayf, reff, border, decomp_levels, 0.05);
        if (ret)
            return ret;
    }

    return 0;
}

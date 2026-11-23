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

/**
 * @file
 * Discrete wavelet transform
 */

#include <string.h>

#include "libavutil/error.h"
#include "libavutil/macros.h"
#include "libavutil/mem.h"
#include "jpeg2000dwt.h"

/* Defines for 9/7 DWT lifting parameters.
 * Parameters are in float. */
#define F_LFTG_ALPHA  1.586134342059924f
#define F_LFTG_BETA   0.052980118572961f
#define F_LFTG_GAMMA  0.882911075530934f
#define F_LFTG_DELTA  0.443506852043971f

/* Lifting parameters in integer format.
 * Computed as param = (float param) * (1 << 16) */
#define I_LFTG_ALPHA_PRIME   38413ll // = 103949 - 65536, (= alpha - 1.0)
#define I_LFTG_BETA           3472ll
#define I_LFTG_GAMMA         57862ll
#define I_LFTG_DELTA         29066ll
#define I_LFTG_K             80621ll
#define I_LFTG_X             53274ll

static inline void extend53(int *p, int i0, int i1)
{
    p[i0 - 1] = p[i0 + 1];
    p[i1]     = p[i1 - 2];
    p[i0 - 2] = p[i0 + 2];
    p[i1 + 1] = p[i1 - 3];
}

static inline void extend97_float(float *p, int i0, int i1)
{
    int i;

    for (i = 1; i <= 4; i++) {
        p[i0 - i]     = p[i0 + i];
        p[i1 + i - 1] = p[i1 - i - 1];
    }
}

static inline void extend97_int(int32_t *p, int i0, int i1)
{
    int i;

    for (i = 1; i <= 4; i++) {
        p[i0 - i]     = p[i0 + i];
        p[i1 + i - 1] = p[i1 - i - 1];
    }
}

static void sd_1d53(int *p, int i0, int i1)
{
    int i;

    if (i1 <= i0 + 1) {
        if (i0 == 1)
            p[1] *= 2;
        return;
    }

    extend53(p, i0, i1);

    for (i = ((i0+1)>>1) - 1; i < (i1+1)>>1; i++)
        p[2*i+1] -= (p[2*i] + p[2*i+2]) >> 1;
    for (i = ((i0+1)>>1); i < (i1+1)>>1; i++)
        p[2*i] += (p[2*i-1] + p[2*i+1] + 2) >> 2;
}

static void dwt_encode53(DWTContext *s, int *t)
{
    int lev,
        w = s->linelen[s->ndeclevels-1][0];
    int *line = s->i_linebuf;
    line += 3;

    for (lev = s->ndeclevels-1; lev >= 0; lev--){
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        int *l;

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;

            for (i = 0; i < lv; i++)
                l[i] = t[w*i + lp];

            sd_1d53(line, mv, mv + lv);

            // copy back and deinterleave
            for (i =   mv; i < lv; i+=2, j++)
                t[w*j + lp] = l[i];
            for (i = 1-mv; i < lv; i+=2, j++)
                t[w*j + lp] = l[i];
        }

        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++){
            int i, j = 0;

            for (i = 0; i < lh; i++)
                l[i] = t[w*lp + i];

            sd_1d53(line, mh, mh + lh);

            // copy back and deinterleave
            for (i =   mh; i < lh; i+=2, j++)
                t[w*lp + j] = l[i];
            for (i = 1-mh; i < lh; i+=2, j++)
                t[w*lp + j] = l[i];
        }
    }
}
static void sd_1d97_float(float *p, int i0, int i1)
{
    int i;

    if (i1 <= i0 + 1) {
        if (i0 == 1)
            p[1] *= F_LFTG_X * 2;
        else
            p[0] *= F_LFTG_K;
        return;
    }

    extend97_float(p, i0, i1);
    i0++; i1++;

    for (i = (i0>>1) - 2; i < (i1>>1) + 1; i++)
        p[2*i+1] -= 1.586134 * (p[2*i] + p[2*i+2]);
    for (i = (i0>>1) - 1; i < (i1>>1) + 1; i++)
        p[2*i] -= 0.052980 * (p[2*i-1] + p[2*i+1]);
    for (i = (i0>>1) - 1; i < (i1>>1); i++)
        p[2*i+1] += 0.882911 * (p[2*i] + p[2*i+2]);
    for (i = (i0>>1); i < (i1>>1); i++)
        p[2*i] += 0.443506 * (p[2*i-1] + p[2*i+1]);
}

static void dwt_encode97_float(DWTContext *s, float *t)
{
    int lev,
        w = s->linelen[s->ndeclevels-1][0];
    float *line = s->f_linebuf;
    line += 5;

    for (lev = s->ndeclevels-1; lev >= 0; lev--){
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        float *l;

        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++){
            int i, j = 0;

            for (i = 0; i < lh; i++)
                l[i] = t[w*lp + i];

            sd_1d97_float(line, mh, mh + lh);

            // copy back and deinterleave
            for (i =   mh; i < lh; i+=2, j++)
                t[w*lp + j] = l[i];
            for (i = 1-mh; i < lh; i+=2, j++)
                t[w*lp + j] = l[i];
        }

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;

            for (i = 0; i < lv; i++)
                l[i] = t[w*i + lp];

            sd_1d97_float(line, mv, mv + lv);

            // copy back and deinterleave
            for (i =   mv; i < lv; i+=2, j++)
                t[w*j + lp] = l[i];
            for (i = 1-mv; i < lv; i+=2, j++)
                t[w*j + lp] = l[i];
        }
    }
}

static void sd_1d97_int(int *p, int i0, int i1)
{
    int i;

    if (i1 <= i0 + 1) {
        if (i0 == 1)
            p[1] = (p[1] * I_LFTG_X + (1<<14)) >> 15;
        else
            p[0] = (p[0] * I_LFTG_K + (1<<15)) >> 16;
        return;
    }

    extend97_int(p, i0, i1);
    i0++; i1++;

    for (i = (i0>>1) - 2; i < (i1>>1) + 1; i++) {
        const int64_t sum = p[2 * i] + p[2 * i + 2];
        p[2 * i + 1] -= sum;
        p[2 * i + 1] -= (I_LFTG_ALPHA_PRIME * sum + (1 << 15)) >> 16;
    }
    for (i = (i0>>1) - 1; i < (i1>>1) + 1; i++)
        p[2 * i]     -= (I_LFTG_BETA  * (p[2 * i - 1] + p[2 * i + 1]) + (1 << 15)) >> 16;
    for (i = (i0>>1) - 1; i < (i1>>1); i++)
        p[2 * i + 1] += (I_LFTG_GAMMA * (p[2 * i]     + p[2 * i + 2]) + (1 << 15)) >> 16;
    for (i = (i0>>1); i < (i1>>1); i++)
        p[2 * i]     += (I_LFTG_DELTA * (p[2 * i - 1] + p[2 * i + 1]) + (1 << 15)) >> 16;
}

static void dwt_encode97_int(DWTContext *s, int *t)
{
    int lev;
    int w = s->linelen[s->ndeclevels-1][0];
    int h = s->linelen[s->ndeclevels-1][1];
    int i;
    int *line = s->i_linebuf;
    line += 5;

    for (i = 0; i < w * h; i++)
        t[i] *= 1 << I_PRESHIFT;

    for (lev = s->ndeclevels-1; lev >= 0; lev--){
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        int *l;

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;

            for (i = 0; i < lv; i++)
                l[i] = t[w*i + lp];

            sd_1d97_int(line, mv, mv + lv);

            // copy back and deinterleave
            for (i =   mv; i < lv; i+=2, j++)
                t[w*j + lp] = l[i];
            for (i = 1-mv; i < lv; i+=2, j++)
                t[w*j + lp] = l[i];
        }

        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++){
            int i, j = 0;

            for (i = 0; i < lh; i++)
                l[i] = t[w*lp + i];

            sd_1d97_int(line, mh, mh + lh);

            // copy back and deinterleave
            for (i =   mh; i < lh; i+=2, j++)
                t[w*lp + j] = l[i];
            for (i = 1-mh; i < lh; i+=2, j++)
                t[w*lp + j] = l[i];
        }

    }

    for (i = 0; i < w * h; i++)
        t[i] = (t[i] + ((1<<(I_PRESHIFT))>>1)) >> (I_PRESHIFT);
}

static void sr_1d53(unsigned *p, int i0, int i1)
{
    int i;

    if (i1 <= i0 + 1) {
        if (i0 == 1)
            p[1] = (int)p[1] >> 1;
        return;
    }

    extend53(p, i0, i1);

    for (i = (i0 >> 1); i < (i1 >> 1) + 1; i++)
        p[2 * i] -= (int)(p[2 * i - 1] + p[2 * i + 1] + 2) >> 2;
    for (i = (i0 >> 1); i < (i1 >> 1); i++)
        p[2 * i + 1] += (int)(p[2 * i] + p[2 * i + 2]) >> 1;
}

static void dwt_decode53(DWTContext *s, int *t)
{
    int lev;
    int w     = s->linelen[s->ndeclevels - 1][0];
    int32_t *line = s->i_linebuf;
    line += 3;

    for (lev = 0; lev < s->ndeclevels; lev++) {
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        int *l;

        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++) {
            int i, j = 0;
            // copy with interleaving
            for (i = mh; i < lh; i += 2, j++)
                l[i] = t[w * lp + j];
            for (i = 1 - mh; i < lh; i += 2, j++)
                l[i] = t[w * lp + j];

            sr_1d53(line, mh, mh + lh);

            for (i = 0; i < lh; i++)
                t[w * lp + i] = l[i];
        }

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;
            // copy with interleaving
            for (i = mv; i < lv; i += 2, j++)
                l[i] = t[w * j + lp];
            for (i = 1 - mv; i < lv; i += 2, j++)
                l[i] = t[w * j + lp];

            sr_1d53(line, mv, mv + lv);

            for (i = 0; i < lv; i++)
                t[w * i + lp] = l[i];
        }
    }
}

static void sr_1d97_float(float *p, int i0, int i1)
{
    if (i1 <= i0 + 1) {
        if (i0 == 1)
            p[1] *= F_LFTG_K/2;
        else
            p[0] *= F_LFTG_X;
        return;
    }

    extend97_float(p, i0, i1);

    /* Steps 1+2 fused in one ascending pass. Step 1 writes E1[i] and step 2
     * reads E1[i] and E1[i+1], so O2[i-1] can be emitted as soon as E1[i] is
     * known; the original O[i-1] is carried from the previous iteration
     * because its slot has already been overwritten by then. Each pass now
     * touches every element once instead of twice, and the operands and
     * evaluation order match the unfused four-pass scalar loops bit for bit
     * under matching FP contraction. */
    const int lo = (i0 >> 1) - 1, hi = i1 >> 1;
    float op = p[2 * lo - 1], e1p = 0.0f;
    int i;

    for (i = lo; i < hi + 2; i++) {
        float o1 = p[2 * i + 1];
        float e1 = p[2 * i] - F_LFTG_DELTA * (op + o1);

        p[2 * i] = e1;
        if (i > lo)
            p[2 * i - 1] = op - F_LFTG_GAMMA * (e1p + e1);

        op   = o1;
        e1p  = e1;
    }

    /* Steps 3+4 fused the same way, reading back E1 and O2. */
    const int lo3 = i0 >> 1, hi3 = i1 >> 1;
    float o2p = p[2 * lo3 - 1], e3p = 0.0f;

    for (i = lo3; i < hi3 + 1; i++) {
        float o2c = p[2 * i + 1];
        float e3  = p[2 * i] + F_LFTG_BETA * (o2p + o2c);

        p[2 * i] = e3;
        if (i > lo3)
            p[2 * i - 1] = o2p + F_LFTG_ALPHA * (e3p + e3);

        o2p = o2c;
        e3p = e3;
    }
}

static void dwt_decode97_float(DWTContext *s, float *t)
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
            /* copy with interleaving, fused into one loop. mh is 0 or 1: the
             * first original loop (l[mh + 2j], count c1) reads the row's first
             * c1 samples, the second (l[1 - mh + 2j], count c2) continues at
             * c1. When mh == 1 the two streams swap which line parity they
             * land on; swapping the source offsets instead keeps the line
             * positions a statically adjacent pair (2j, 2j+1) so the compiler
             * can emit interleaved stores instead of per-lane scatters. The
             * values are moved unmodified, so this is bit-exact by
             * construction. */
            int j, c1 = (lh - mh + 1) >> 1, c2 = (lh + mh) >> 1;
            int eoff = mh ? c1 : 0, ooff = mh ? 0 : c1;
            int cE = mh ? c2 : c1, cO = mh ? c1 : c2;
            const int cmin = FFMIN(cE, cO);
            const float *restrict src = data + (size_t)w * lp;
            const float *restrict srcE = src + eoff;
            const float *restrict srcO = src + ooff;
            float *restrict dst = l;

            /* Four index pairs per iteration: the two loads become one
             * 128-bit load per source, and the two interleaved stores one
             * 256-bit store, so the permutation runs at memory speed
             * instead of one scalar load/store pair per coefficient. */
            for (j = 0; j + 4 <= cmin; j += 4) {
                float e0 = srcE[j], e1 = srcE[j + 1];
                float e2 = srcE[j + 2], e3 = srcE[j + 3];
                float o0 = srcO[j], o1 = srcO[j + 1];
                float o2 = srcO[j + 2], o3 = srcO[j + 3];

                dst[2 * j]      = e0;
                dst[2 * j + 1]  = o0;
                dst[2 * j + 2]  = e1;
                dst[2 * j + 3]  = o1;
                dst[2 * j + 4]  = e2;
                dst[2 * j + 5]  = o2;
                dst[2 * j + 6]  = e3;
                dst[2 * j + 7]  = o3;
            }
            for (; j < cmin; j++) {
                dst[2 * j]     = srcE[j];
                dst[2 * j + 1] = srcO[j];
            }
            for (; j < cE; j++)
                dst[2 * j] = srcE[j];
            for (; j < cO; j++)
                dst[2 * j + 1] = srcO[j];

            sr_1d97_float(line, mh, mh + lh);

            /* The lifted line is moved back unchanged; memcpy lets the
             * compiler emit a vector copy instead of an aliasing-limited
             * scalar loop. */
            memcpy(data + (size_t)w * lp, l, (size_t)lh * sizeof(float));
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

            sr_1d97_float(line, mv, mv + lv);

            for (i = 0; i < lv; i++)
                data[w * i + lp] = l[i];
        }
    }
}

static void sr_1d97_int(int32_t *p, int i0, int i1)
{
    int i;

    if (i1 <= i0 + 1) {
        if (i0 == 1)
            p[1] = (p[1] * I_LFTG_K + (1<<16)) >> 17;
        else
            p[0] = (p[0] * I_LFTG_X + (1<<15)) >> 16;
        return;
    }

    extend97_int(p, i0, i1);

    for (i = (i0 >> 1) - 1; i < (i1 >> 1) + 2; i++)
        p[2 * i]     -= (I_LFTG_DELTA * (p[2 * i - 1] + (int64_t)p[2 * i + 1]) + (1 << 15)) >> 16;
    /* step 4 */
    for (i = (i0 >> 1) - 1; i < (i1 >> 1) + 1; i++)
        p[2 * i + 1] -= (I_LFTG_GAMMA * (p[2 * i]     + (int64_t)p[2 * i + 2]) + (1 << 15)) >> 16;
    /*step 5*/
    for (i = (i0 >> 1); i < (i1 >> 1) + 1; i++)
        p[2 * i]     += (I_LFTG_BETA  * (p[2 * i - 1] + (int64_t)p[2 * i + 1]) + (1 << 15)) >> 16;
    /* step 6 */
    for (i = (i0 >> 1); i < (i1 >> 1); i++) {
        const int64_t sum = p[2 * i] + (int64_t) p[2 * i + 2];
        p[2 * i + 1] += sum;
        p[2 * i + 1] += (I_LFTG_ALPHA_PRIME * sum + (1 << 15)) >> 16;
    }
}

static void dwt_decode97_int(DWTContext *s, int32_t *t)
{
    int lev;
    int w       = s->linelen[s->ndeclevels - 1][0];
    int h       = s->linelen[s->ndeclevels - 1][1];
    int i;
    int32_t *line = s->i_linebuf;
    int32_t *data = t;
    /* position at index O of line range [0-5,w+5] cf. extend function */
    line += 5;

    for (lev = 0; lev < s->ndeclevels; lev++) {
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        int32_t *l;
        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++) {
            int i, j = 0;
            // interleaving
            for (i = mh; i < lh; i += 2, j++)
                l[i] = data[w * lp + j];
            for (i = 1 - mh; i < lh; i += 2, j++)
                l[i] = data[w * lp + j];

            sr_1d97_int(line, mh, mh + lh);

            for (i = 0; i < lh; i++)
                data[w * lp + i] = l[i];
        }

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;
            // interleaving
            for (i = mv; i < lv; i += 2, j++)
                l[i] = data[w * j + lp];
            for (i = 1 - mv; i < lv; i += 2, j++)
                l[i] = data[w * j + lp];

            sr_1d97_int(line, mv, mv + lv);

            for (i = 0; i < lv; i++)
                data[w * i + lp] = l[i];
        }
    }

    for (i = 0; i < w * h; i++)
        // We shift down by `I_PRESHIFT` because the input coefficients `datap[]` were shifted up by `I_PRESHIFT` to secure the precision
        data[i] = (int32_t)(data[i] + ((1LL<<(I_PRESHIFT))>>1)) >> (I_PRESHIFT);
}

int ff_jpeg2000_dwt_init(DWTContext *s, int border[2][2],
                         int decomp_levels, int type)
{
    int i, j, lev = decomp_levels, maxlen,
        b[2][2];

    s->ndeclevels = decomp_levels;
    s->type       = type;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            b[i][j] = border[i][j];

    maxlen = FFMAX(b[0][1] - b[0][0],
                   b[1][1] - b[1][0]);
    while (--lev >= 0)
        for (i = 0; i < 2; i++) {
            s->linelen[lev][i] = b[i][1] - b[i][0];
            s->mod[lev][i]     = b[i][0] & 1;
            for (j = 0; j < 2; j++)
                b[i][j] = (b[i][j] + 1) >> 1;
        }
    switch (type) {
    case FF_DWT97:
        s->f_linebuf = av_malloc_array((maxlen + 12), sizeof(*s->f_linebuf));
        if (!s->f_linebuf)
            return AVERROR(ENOMEM);
        break;
     case FF_DWT97_INT:
        s->i_linebuf = av_malloc_array((maxlen + 12), sizeof(*s->i_linebuf));
        if (!s->i_linebuf)
            return AVERROR(ENOMEM);
        break;
    case FF_DWT53:
        s->i_linebuf = av_malloc_array((maxlen +  6), sizeof(*s->i_linebuf));
        if (!s->i_linebuf)
            return AVERROR(ENOMEM);
        break;
    default:
        return -1;
    }
    return 0;
}

int ff_dwt_encode(DWTContext *s, void *t)
{
    if (s->ndeclevels == 0)
        return 0;

    switch(s->type){
        case FF_DWT97:
            dwt_encode97_float(s, t); break;
        case FF_DWT97_INT:
            dwt_encode97_int(s, t); break;
        case FF_DWT53:
            dwt_encode53(s, t); break;
        default:
            return -1;
    }
    return 0;
}

int ff_dwt_decode(DWTContext *s, void *t)
{
    if (s->ndeclevels == 0)
        return 0;

    switch (s->type) {
    case FF_DWT97:
        dwt_decode97_float(s, t);
        break;
    case FF_DWT97_INT:
        dwt_decode97_int(s, t);
        break;
    case FF_DWT53:
        dwt_decode53(s, t);
        break;
    default:
        return -1;
    }
    return 0;
}

void ff_dwt_destroy(DWTContext *s)
{
    av_freep(&s->f_linebuf);
    av_freep(&s->i_linebuf);
}

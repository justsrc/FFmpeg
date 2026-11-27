/*
 * MQ-coder differential test
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
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "libavutil/common.h"
#include "libavutil/lfg.h"

#include "libavcodec/mqc.h"

#define DATA_SIZE   4096
#define BUF_SIZE    8192
#define NUM_TRIALS  1024
#define NUM_SYMBOLS 4096

/* The packed table is hand-built; check every entry against the scalar
 * tables, including states the decoder may never reach. */
static int check_packed_table(void)
{
    int i;

    for (i = 0; i < 2 * 47; i++) {
        uint32_t expected = ff_mqc_qe[i] |
                            ((uint32_t)(i & 1) << 16) |
                            ((uint32_t)ff_mqc_nmps[i] << 17) |
                            ((uint32_t)ff_mqc_nlps[i] << 24);

        if (ff_mqc_packed[i] != expected) {
            fprintf(stderr,
                    "packed[%d] mismatch: 0x%08"PRIx32" != 0x%08"PRIx32"\n",
                    i, ff_mqc_packed[i], expected);
            return 1;
        }
    }

    return 0;
}

/* Scalar reference decoder copied from the pre-rewrite mqcdec.c. */
static void bytein_ref(MqcState *mqc)
{
    if (*mqc->bp == 0xff) {
        if (*(mqc->bp + 1) > 0x8f)
            mqc->c++;
        else {
            mqc->bp++;
            mqc->c += 2 + 0xfe00 - (*mqc->bp << 9);
        }
    } else {
        mqc->bp++;
        mqc->c += 1 + 0xff00 - (*mqc->bp << 8);
    }
}

static int exchange_ref(MqcState *mqc, uint8_t *cxstate, int lps)
{
    int d;
    if ((mqc->a < ff_mqc_qe[*cxstate]) ^ (!lps)) {
        if (lps)
            mqc->a = ff_mqc_qe[*cxstate];
        d = *cxstate & 1;
        *cxstate = ff_mqc_nmps[*cxstate];
    } else {
        if (lps)
            mqc->a = ff_mqc_qe[*cxstate];
        d = 1 - (*cxstate & 1);
        *cxstate = ff_mqc_nlps[*cxstate];
    }
    // do RENORMD: see ISO/IEC 15444-1:2002 §C.3.3
    do {
        if (!(mqc->c & 0xff)) {
            mqc->c -= 0x100;
            bytein_ref(mqc);
        }
        mqc->a += mqc->a;
        mqc->c += mqc->c;
    } while (!(mqc->a & 0x8000));
    return d;
}

static void ff_mqc_initdec_ref(MqcState *mqc, uint8_t *bp, int raw, int reset)
{
    mqc->raw = raw;
    if (reset)
        ff_mqc_init_contexts(mqc);
    mqc->bp = bp;
    mqc->c  = (*mqc->bp ^ 0xff) << 16;
    bytein_ref(mqc);
    mqc->c = mqc->c << 7;
    mqc->a = 0x8000;
}

static int mqc_decode_bypass_ref(MqcState *mqc)
{
    int bit = !(mqc->c & 0x40000000);
    if (!(mqc->c & 0xff)) {
        mqc->c -= 0x100;
        bytein_ref(mqc);
    }
    mqc->c += mqc->c;
    return bit;
}

static int ff_mqc_decode_ref(MqcState *mqc, uint8_t *cxstate)
{
    if (mqc->raw)
        return mqc_decode_bypass_ref(mqc);
    mqc->a -= ff_mqc_qe[*cxstate];
    if ((mqc->c >> 16) < mqc->a) {
        if (mqc->a & 0x8000)
            return *cxstate & 1;
        else
            return exchange_ref(mqc, cxstate, 0);
    } else {
        mqc->c -= mqc->a << 16;
        return exchange_ref(mqc, cxstate, 1);
    }
}

static void fill_buffer(uint8_t *buf, AVLFG *prng, int trial)
{
    static const uint8_t stuffed_tails[] = { 0x00, 0x7f, 0x80, 0x8f };
    static const uint8_t marker_tails[]  = { 0x90, 0xd9, 0xfe, 0xff };
    int i;

    memset(buf, 0, BUF_SIZE);

    if ((trial & 3) == 0) {
        /* Sparse/repetitive patterns exercise long MPS and LPS runs. */
        uint8_t fill = (trial & 4) ? 0xff : 0x00;
        memset(buf, fill, DATA_SIZE);
        for (i = 0; i < 32; i++) {
            int pos = av_lfg_get(prng) % (DATA_SIZE - 4);
            buf[pos]     = 0xff;
            buf[pos + 1] = (trial & 8) ? marker_tails[i & 3] : stuffed_tails[i & 3];
        }
    } else {
        for (i = 0; i < DATA_SIZE; i++)
            buf[i] = av_lfg_get(prng) & 0xff;

        /* On half of the random trials, keep 0xFF bytes stuffed so bp crosses
         * the entire 4 KiB payload into the zero tail. */
        if ((trial & 3) == 1) {
            for (i = 0; i < DATA_SIZE - 1; i++) {
                if (buf[i] == 0xff && buf[i + 1] > 0x8f)
                    buf[i + 1] &= 0x7f;
            }
        }

        for (i = 0; i < 16; i++) {
            int pos = av_lfg_get(prng) % (DATA_SIZE - 8);
            int kind = (trial & 3) == 1 ? (i & 3) : (i & 7);

            if (kind < 4) {
                buf[pos]     = 0xff;
                buf[pos + 1] = stuffed_tails[kind];
            } else if (kind < 7) {
                buf[pos]     = 0xff;
                buf[pos + 1] = marker_tails[kind - 4];
            } else {
                buf[pos]     = 0xff;
                buf[pos + 1] = 0xff;
                buf[pos + 2] = 0xff;
                buf[pos + 3] = stuffed_tails[i & 3];
            }
        }
    }

    /* Exercise deliberate 0xFF pairs right at the stream start and tail boundary. */
    switch (trial) {
    case 0:
        buf[0] = 0xff;
        buf[1] = 0x8f;
        break;
    case 1:
        buf[0] = 0xff;
        buf[1] = 0x00;
        break;
    case 2:
        buf[0] = 0xff;
        buf[1] = 0x90;
        break;
    case 3:
        buf[0] = 0xff;
        buf[1] = 0xff;
        break;
    case 4:
        buf[DATA_SIZE - 1] = 0xff;
        break;
    }
}

static int run_trials(int raw, AVLFG *prng)
{
    uint8_t buf[BUF_SIZE];
    int trial, sym, k;

    for (trial = 0; trial < NUM_TRIALS; trial++) {
        MqcState dec, ref;
        int reset = trial & 1;

        fill_buffer(buf, prng, trial);

        memset(&dec, 0, sizeof(dec));
        memset(&ref, 0, sizeof(ref));
        ff_mqc_init_contexts(&dec);
        ff_mqc_init_contexts(&ref);

        ff_mqc_initdec(&dec, buf, raw, reset);
        ff_mqc_initdec_ref(&ref, buf, raw, reset);

        /* Seed after initdec() so the injected states survive a reset. The
         * eight-trial groups walk 19 consecutive entries each; over the 128
         * groups every packed state is used at least once. */
        if ((trial & 7) >= 6) {
            for (k = 0; k < 19; k++) {
                uint8_t st = ((trial >> 3) * 19 + k) % 94;
                dec.cx_states[k] = st;
                ref.cx_states[k] = st;
            }
        }

        if (memcmp(dec.cx_states, ref.cx_states, sizeof(dec.cx_states))) {
            fprintf(stderr, "init context mismatch (raw=%d trial=%d reset=%d)\n",
                    raw, trial, reset);
            return 1;
        }

        for (sym = 0; sym < NUM_SYMBOLS; sym++) {
            int cx = ((trial & 3) == 0) ? ((trial >> 2) % 19)
                                        : (av_lfg_get(prng) % 19);
            int bit = (trial & 2)
                ? ff_mqc_decode(&dec, dec.cx_states + cx)
                : ff_mqc_decode_inline(&dec, dec.cx_states + cx);
            int bit_ref = ff_mqc_decode_ref(&ref, ref.cx_states + cx);

            if (bit != bit_ref || dec.a != ref.a) {
                fprintf(stderr,
                        "mismatch (raw=%d trial=%d sym=%d cx=%d): "
                        "bit=%d ref_bit=%d a=0x%x ref_a=0x%x\n",
                        raw, trial, sym, cx, bit, bit_ref, dec.a, ref.a);
                return 1;
            }

            if (dec.bp < ref.bp) {
                fprintf(stderr,
                        "bp behind reference (raw=%d trial=%d sym=%d): "
                        "bp=%d ref_bp=%d\n",
                        raw, trial, sym, (int)(dec.bp - buf),
                        (int)(ref.bp - buf));
                return 1;
            }

            if (raw) {
                if (dec.c != ref.c) {
                    fprintf(stderr,
                            "raw c mismatch (trial=%d sym=%d): c=0x%x ref_c=0x%x\n",
                            trial, sym, dec.c, ref.c);
                    return 1;
                }
            } else {
                if ((unsigned)(dec.c64 >> 48) != (unsigned)(ref.c >> 16)) {
                    fprintf(stderr,
                            "code mismatch (trial=%d sym=%d): c=0x%x ref_c=0x%x\n",
                            trial, sym, (unsigned)(dec.c64 >> 48),
                            (unsigned)(ref.c >> 16));
                    return 1;
                }
                for (k = 0; k < 19; k++) {
                    if (dec.cx_states[k] != ref.cx_states[k]) {
                        fprintf(stderr,
                                "cx_states[%d] mismatch (trial=%d sym=%d cx=%d): "
                                "%u != %u\n",
                                k, trial, sym, cx,
                                dec.cx_states[k], ref.cx_states[k]);
                        return 1;
                    }
                }
            }
        }
    }

    return 0;
}

static void fill_sentinel_buffer(uint8_t *buf, int len, AVLFG *prng, int trial)
{
    int i;

    for (i = 0; i < len; i++) {
        if ((trial & 3) == 0)
            buf[i] = 0x00;
        else if ((trial & 3) == 1)
            buf[i] = 0xff;
        else
            buf[i] = av_lfg_get(prng) & 0xff;
    }

    if (trial & 4) {
        for (i = 0; i + 1 < len; i += 5) {
            buf[i]     = 0xff;
            buf[i + 1] = 0xff;
        }
    }

    /* jpeg2000dec sets data[length] and data[length+1] to 0xff and zeroes
     * the four bytes of padding after the block payload. */
    buf[len]     = 0xff;
    buf[len + 1] = 0xff;
    buf[len + 2] = 0x00;
    buf[len + 3] = 0x00;
}

/* Exact-length buffers ending in the cblk 0xff 0xff sentinel, decoded until
 * the read-ahead reaches the tail. */
static int run_sentinel_trials(AVLFG *prng)
{
    static const int lens[] = { 1, 2, 3, 7, 8, 15, 16, 31, 32, 63, 64,
                                127, 128, 1024, 2048, 4096 };
    uint8_t buf[4096 + 4];
    size_t li;
    int trial, sym, k;

    for (li = 0; li < FF_ARRAY_ELEMS(lens); li++) {
        int len = lens[li];
        int trials = len <= 128 ? 8 : 2;

        for (trial = 0; trial < trials; trial++) {
            MqcState dec, ref;

            fill_sentinel_buffer(buf, len, prng, trial);

            memset(&dec, 0, sizeof(dec));
            memset(&ref, 0, sizeof(ref));
            ff_mqc_init_contexts(&dec);
            ff_mqc_init_contexts(&ref);
            ff_mqc_initdec(&dec, buf, 0, 1);
            ff_mqc_initdec_ref(&ref, buf, 0, 1);

            for (sym = 0; sym < 8 * len + 64; sym++) {
                int cx = av_lfg_get(prng) % 19;
                int bit = (sym & 1)
                    ? ff_mqc_decode_inline(&dec, dec.cx_states + cx)
                    : ff_mqc_decode(&dec, dec.cx_states + cx);
                int bit_ref = ff_mqc_decode_ref(&ref, ref.cx_states + cx);

                if (bit != bit_ref || dec.a != ref.a ||
                    (unsigned)(dec.c64 >> 48) != (unsigned)(ref.c >> 16)) {
                    fprintf(stderr,
                            "sentinel mismatch (len=%d trial=%d sym=%d cx=%d): "
                            "bit=%d ref_bit=%d a=0x%x ref_a=0x%x "
                            "c=0x%x ref_c=0x%x\n",
                            len, trial, sym, cx, bit, bit_ref, dec.a, ref.a,
                            (unsigned)(dec.c64 >> 48),
                            (unsigned)(ref.c >> 16));
                    return 1;
                }

                for (k = 0; k < 19; k++) {
                    if (dec.cx_states[k] != ref.cx_states[k]) {
                        fprintf(stderr,
                                "sentinel cx_states[%d] mismatch "
                                "(len=%d trial=%d sym=%d cx=%d): %u != %u\n",
                                k, len, trial, sym, cx,
                                dec.cx_states[k], ref.cx_states[k]);
                        return 1;
                    }
                }
            }

            if (dec.bp > buf + len || ref.bp > buf + len) {
                fprintf(stderr,
                        "sentinel overrun (len=%d trial=%d): "
                        "bp=%d ref_bp=%d\n",
                        len, trial, (int)(dec.bp - buf),
                        (int)(ref.bp - buf));
                return 1;
            }
        }
    }

    return 0;
}

int main(void)
{
    AVLFG prng;

    if (check_packed_table())
        return 1;

    av_lfg_init(&prng, 1);

    if (run_trials(0, &prng))
        return 1;
    if (run_trials(1, &prng))
        return 1;
    if (run_sentinel_trials(&prng))
        return 1;

    return 0;
}

/*
 * MQ-coder: structures, common and decoder functions
 * Copyright (c) 2007 Kamil Nowosad
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

#ifndef AVCODEC_MQC_H
#define AVCODEC_MQC_H

/**
 * MQ-coder
 * @file
 * @author Kamil Nowosad
 */

#include <stdint.h>

#include "libavutil/attributes.h"
#include "libavutil/intmath.h"

#define MQC_CX_UNI 17
#define MQC_CX_RL  18

extern const uint16_t ff_mqc_qe[2 * 47];
extern const uint8_t  ff_mqc_nlps[2 * 47];
extern const uint8_t  ff_mqc_nmps[2 * 47];

/*
 * Probability states packed for the decoder so that one load yields the
 * value, the MPS and both transitions: Qe in bits 0-15, MPS in bit 16, the
 * next MPS state in bits 17-23 and the next LPS state in bits 24-30.
 */
extern const uint32_t ff_mqc_packed[2 * 47];

typedef struct MqcState {
    uint8_t *bp, *bpstart;
    unsigned int a;
    unsigned int c;
    unsigned int ct;
    uint8_t cx_states[19];
    int raw;
    uint64_t c64;         ///< code value in bits 48-63, read-ahead below
    unsigned shift_count; ///< valid read-ahead bits below the code value
} MqcState;

/*
 * Decoder. In MQ mode the C-register holds the code value in bits 48-63 and
 * up to 48 bits of read-ahead below it, so renormalization shifts once for
 * all the bits it needs instead of once per bit, and the symbol decision is
 * branchless. This is the form used by T.88 Annex G.
 *
 * The tier-1 loops call this millions of times per frame; with the state in
 * a caller-local MqcState the compiler keeps C, A and the byte pointer in
 * registers instead of reloading them around a cross-TU call.
 */
static av_always_inline void ff_mqc_bytein_inline(MqcState *mqc)
{
    /* BYTEIN, repeated while another byte fits below the code value. A byte
     * after 0xFF may carry into the bits above it, so it is only read once no
     * read-ahead bits remain; any other byte lands on zero bits and can be
     * read early without changing what the cell sees. */
    while (mqc->shift_count <= 40) {
        if (*mqc->bp == 0xff) {
            if (mqc->bp[1] > 0x8f) {
                mqc->shift_count += 8;
            } else if (!mqc->shift_count) {
                mqc->bp++;
                mqc->c64 += (0x7fULL << 41) - ((uint64_t)*mqc->bp << 41);
                mqc->shift_count = 7;
            } else
                break;
        } else {
            mqc->bp++;
            mqc->c64 += (uint64_t)(*mqc->bp ^ 0xff) << (40 - mqc->shift_count);
            mqc->shift_count += 8;
        }
    }
}

/* RENORMD: shift A back to at least 0x8000, batching the shifts the
 * read-ahead bits cover. */
static av_always_inline void ff_mqc_renorm_inline(MqcState *mqc)
{
    while (!(mqc->a & 0x8000)) {
        unsigned shift;

        ff_mqc_bytein_inline(mqc);
        shift = ff_clz(mqc->a) - 16;
        if (shift > mqc->shift_count)
            shift = mqc->shift_count;
        mqc->a       <<= shift;
        mqc->c64     <<= shift;
        mqc->shift_count -= shift;
    }
}

/* Byte input of the 32-bit C-register, kept for the raw (bypass) segments. */
static av_always_inline void ff_mqc_bytein32_inline(MqcState *mqc)
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

static av_always_inline int ff_mqc_decode_bypass_inline(MqcState *mqc)
{
    int bit = !(mqc->c & 0x40000000);

    if (!(mqc->c & 0xff)) {
        mqc->c -= 0x100;
        ff_mqc_bytein32_inline(mqc);
    }
    mqc->c += mqc->c;

    return bit;
}

static av_always_inline int ff_mqc_decode_inline(MqcState *mqc,
                                                 uint8_t *cxstate)
{
    uint32_t state, qe, a, is_lps, mask, cond, pick_nlps;
    uint64_t mask64;
    unsigned shift;
    uint8_t *cx = cxstate;

    if (mqc->raw)
        return ff_mqc_decode_bypass_inline(mqc);

    state = ff_mqc_packed[*cx];
    qe    = state & 0xffff;
    a     = mqc->a - qe;

    /* The MPS/LPS decision, MPS_EXCHANGE and LPS_EXCHANGE of C.3.2, without
     * branches: both exchanges differ only in whether the condition flips
     * the role of the exchange, so one XOR handles both. */
    is_lps = (mqc->c64 >> 48) >= a;
    mask   = 0 - is_lps;
    mask64 = 0 - (uint64_t)is_lps;
    mqc->c64 -= ((uint64_t)a << 48) & mask64;
    mqc->a    = (a & ~mask) | (qe & mask);

    cond      = a < qe;
    pick_nlps = cond ^ is_lps;
    /* A probability-state transition happens exactly when A renormalizes:
     * on every LPS, and on an MPS that leaves A below 0x8000. */
    *cx = (is_lps | (a < 0x8000))
        ? (state >> (17 + 7 * pick_nlps)) & 0x7f
        : *cx;

    /* Zero shifts in the common case, where A >= 0x8000 already. */
    shift = ff_clz(mqc->a) - 16;
    if (shift <= mqc->shift_count) {
        mqc->a       <<= shift;
        mqc->c64     <<= shift;
        mqc->shift_count -= shift;
    } else {
        ff_mqc_renorm_inline(mqc);
    }

    return ((state >> 16) & 1) ^ pick_nlps;
}

/* encoder */

/** initialize the encoder */
void ff_mqc_initenc(MqcState *mqc, uint8_t *bp);

/** code bit d with context cx */
void ff_mqc_encode(MqcState *mqc, uint8_t *cxstate, int d);

/** flush the encoder [returns number of bytes encoded] */
int ff_mqc_flush_to(MqcState *mqc, uint8_t *dst, int *dst_len);

/* decoder */

/**
 * Initialize MQ-decoder.
 * @param mqc   MQ decoder state
 * @param bp    byte pointer
 * @param raw   raw mode
 * @param reset reset states
 */
void ff_mqc_initdec(MqcState *mqc, uint8_t *bp, int raw, int reset);

/**
 * MQ decoder.
 * @param mqc       MQ decoder state
 * @param cxstate   Context
 * @return          Decision (0 to 1)
 */
int ff_mqc_decode(MqcState *mqc, uint8_t *cxstate);

/* common */

/**
 * MQ-coder context initialisations.
 * @param mqc       MQ-coder context
 */
void ff_mqc_init_contexts(MqcState *mqc);

#endif /* AVCODEC_MQC_H */

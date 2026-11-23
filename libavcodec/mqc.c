/*
 * MQ-coder encoder and decoder common functions
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

/**
 * MQ-coder common (decoder/encoder) functions
 * @file
 * @author Kamil Nowosad
 */

#include <string.h>
#include <stdint.h>

#include "mqc.h"

const uint16_t ff_mqc_qe[2 * 47] = {
    0x5601, 0x5601, 0x3401, 0x3401, 0x1801, 0x1801, 0x0ac1, 0x0ac1,
    0x0521, 0x0521, 0x0221, 0x0221, 0x5601, 0x5601, 0x5401, 0x5401,
    0x4801, 0x4801, 0x3801, 0x3801, 0x3001, 0x3001, 0x2401, 0x2401,
    0x1c01, 0x1c01, 0x1601, 0x1601, 0x5601, 0x5601, 0x5401, 0x5401,
    0x5101, 0x5101, 0x4801, 0x4801, 0x3801, 0x3801, 0x3401, 0x3401,
    0x3001, 0x3001, 0x2801, 0x2801, 0x2401, 0x2401, 0x2201, 0x2201,
    0x1c01, 0x1c01, 0x1801, 0x1801, 0x1601, 0x1601, 0x1401, 0x1401,
    0x1201, 0x1201, 0x1101, 0x1101, 0x0ac1, 0x0ac1, 0x09c1, 0x09c1,
    0x08a1, 0x08a1, 0x0521, 0x0521, 0x0441, 0x0441, 0x02a1, 0x02a1,
    0x0221, 0x0221, 0x0141, 0x0141, 0x0111, 0x0111, 0x0085, 0x0085,
    0x0049, 0x0049, 0x0025, 0x0025, 0x0015, 0x0015, 0x0009, 0x0009,
    0x0005, 0x0005, 0x0001, 0x0001, 0x5601, 0x5601
};
const uint8_t ff_mqc_nlps[2 * 47] = {
     3,  2, 12, 13, 18, 19, 24, 25, 58, 59, 66, 67, 13, 12, 28, 29,
    28, 29, 28, 29, 34, 35, 36, 37, 40, 41, 42, 43, 29, 28, 28, 29,
    30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 38, 39, 40, 41, 42, 43,
    44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
    60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,
    76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 92, 93
};
const uint8_t ff_mqc_nmps[2 * 47] = {
     2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 76, 77, 14, 15, 16, 17,
    18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 58, 59, 30, 31, 32, 33,
    34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65,
    66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,
    82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 90, 91, 92, 93
};

/*
 * Qe, MPS, next-MPS and next-LPS state of each context state, packed so the
 * decoder reads all four with one 32-bit load.
 */
const uint32_t ff_mqc_packed[2 * 47] = {
    0x03045601, 0x02075601, 0x0c083401, 0x0d0b3401,
    0x120c1801, 0x130f1801, 0x18100ac1, 0x19130ac1,
    0x3a140521, 0x3b170521, 0x42980221, 0x439b0221,
    0x0d1c5601, 0x0c1f5601, 0x1c205401, 0x1d235401,
    0x1c244801, 0x1d274801, 0x1c283801, 0x1d2b3801,
    0x222c3001, 0x232f3001, 0x24302401, 0x25332401,
    0x28341c01, 0x29371c01, 0x2a741601, 0x2b771601,
    0x1d3c5601, 0x1c3f5601, 0x1c405401, 0x1d435401,
    0x1e445101, 0x1f475101, 0x20484801, 0x214b4801,
    0x224c3801, 0x234f3801, 0x24503401, 0x25533401,
    0x26543001, 0x27573001, 0x26582801, 0x275b2801,
    0x285c2401, 0x295f2401, 0x2a602201, 0x2b632201,
    0x2c641c01, 0x2d671c01, 0x2e681801, 0x2f6b1801,
    0x306c1601, 0x316f1601, 0x32701401, 0x33731401,
    0x34741201, 0x35771201, 0x36781101, 0x377b1101,
    0x387c0ac1, 0x397f0ac1, 0x3a8009c1, 0x3b8309c1,
    0x3c8408a1, 0x3d8708a1, 0x3e880521, 0x3f8b0521,
    0x408c0441, 0x418f0441, 0x429002a1, 0x439302a1,
    0x44940221, 0x45970221, 0x46980141, 0x479b0141,
    0x489c0111, 0x499f0111, 0x4aa00085, 0x4ba30085,
    0x4ca40049, 0x4da70049, 0x4ea80025, 0x4fab0025,
    0x50ac0015, 0x51af0015, 0x52b00009, 0x53b30009,
    0x54b40005, 0x55b70005, 0x56b40001, 0x57b70001,
    0x5cb85601, 0x5dbb5601,
};

void ff_mqc_init_contexts(MqcState *mqc)
{
    memset(mqc->cx_states, 0, sizeof(mqc->cx_states));
    mqc->cx_states[MQC_CX_UNI] = 2 * 46;
    mqc->cx_states[MQC_CX_RL]  = 2 * 3;
    mqc->cx_states[0]          = 2 * 4;
}

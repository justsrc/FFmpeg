/*
 * MQ-coder decoder
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
 * MQ-coder decoder
 * @file
 * @author Kamil Nowosad
 */

#include "mqc.h"

void ff_mqc_initdec(MqcState *mqc, uint8_t *bp, int raw, int reset)
{
    mqc->raw = raw;
    if (reset)
        ff_mqc_init_contexts(mqc);
    mqc->bp = bp;
    if (raw) {
        mqc->c  = (*mqc->bp ^ 0xff) << 16;
        ff_mqc_bytein32_inline(mqc);
        mqc->c = mqc->c << 7;
    } else {
        mqc->c64 = (uint64_t)(*mqc->bp ^ 0xff) << 48;
        mqc->shift_count = 0;
        ff_mqc_bytein_inline(mqc);
        mqc->c64 <<= 7;
        mqc->shift_count -= 7;
    }
    mqc->a = 0x8000;
}

int ff_mqc_decode(MqcState *mqc, uint8_t *cxstate)
{
    return ff_mqc_decode_inline(mqc, cxstate);
}

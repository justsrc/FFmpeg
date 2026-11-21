/*
 * Copyright (c) 2026 Justin Jacob
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

#include "libavutil/aarch64/asm.S"

/*
 * Inverse irreversible color transform (float), FF_DWT97.
 * x0 = src0, x1 = src1, x2 = src2, w3 = csize (float samples per plane)
 *
 * Bit-exact to the C implementation as compiled on arm64 with the default
 * FP contraction: every multiply-accumulate is fused (fmla), and the
 * negations of the two subtraction terms are folded into the constants,
 * replicating the per-sample operation order of ict_float() in
 * libavcodec/jpeg2000dsp.c.
 */

function ff_ict_float_neon, export=1
        cbz             w3,  8f
        mov             w8,  #0x74bc
        movk            w8,  #0x3fb3, lsl #16
        dup             v16.4s, w8           //  1.402
        mov             w8,  #0x31cf
        movk            w8,  #0xbeb0, lsl #16
        dup             v17.4s, w8           // -0.34413
        mov             w8,  #0xd1e1
        movk            w8,  #0xbf36, lsl #16
        dup             v18.4s, w8           // -0.71414
        mov             w8,  #0xd0e5
        movk            w8,  #0x3fe2, lsl #16
        dup             v19.4s, w8           //  1.772
        bic             w4,  w3, #3
        cbz             w4,  2f
1:      ld1             {v0.4s}, [x0]
        ld1             {v1.4s}, [x1]
        ld1             {v2.4s}, [x2]
        subs            w4,  w4, #4
        mov             v3.16b, v0.16b
        fmla            v3.4s, v1.4s, v17.4s // i1 = src0 - 0.34413*src1
        fmla            v3.4s, v2.4s, v18.4s // i1 = i1     - 0.71414*src2
        mov             v4.16b, v0.16b
        fmla            v4.4s, v1.4s, v19.4s // i2 = src0 + 1.772*src1
        fmla            v0.4s, v2.4s, v16.4s // i0 = src0 + 1.402*src2
        st1             {v0.4s}, [x0], #16
        st1             {v3.4s}, [x1], #16
        st1             {v4.4s}, [x2], #16
        b.ne            1b
2:      and             w3,  w3, #3
        cbz             w3,  8f
3:      ldr             s20, [x0]
        ldr             s21, [x1]
        ldr             s22, [x2]
        fmadd           s23, s21, s17, s20   // t  = src0 - 0.34413*src1
        fmadd           s23, s22, s18, s23   // i1 = t    - 0.71414*src2
        fmadd           s24, s21, s19, s20   // i2 = src0 + 1.772*src1
        fmadd           s20, s22, s16, s20   // i0 = src0 + 1.402*src2
        str             s20, [x0], #4
        str             s23, [x1], #4
        str             s24, [x2], #4
        subs            w3,  w3, #1
        b.ne            3b
8:      ret
endfunc

/*
 * Reversible color transform (integer), FF_DWT53.
 * x0 = src0, x1 = src1, x2 = src2, w3 = csize (int32 samples per plane)
 */

function ff_rct_int_neon, export=1
        cbz             w3,  8f
        bic             w4,  w3, #3
        cbz             w4,  2f
1:      ld1             {v0.4s}, [x0]
        ld1             {v1.4s}, [x1]
        ld1             {v2.4s}, [x2]
        subs            w4,  w4, #4
        add             v3.4s, v1.4s, v2.4s  // src1 + src2
        sshr            v3.4s, v3.4s, #2
        sub             v3.4s, v0.4s, v3.4s  // i1
        add             v0.4s, v3.4s, v2.4s  // i0
        add             v2.4s, v3.4s, v1.4s  // i2
        st1             {v0.4s}, [x0], #16
        st1             {v3.4s}, [x1], #16
        st1             {v2.4s}, [x2], #16
        b.ne            1b
2:      and             w3,  w3, #3
        cbz             w3,  8f
3:      ldr             w5,  [x0]
        ldr             w6,  [x1]
        ldr             w7,  [x2]
        add             w9,  w6, w7
        asr             w9,  w9, #2
        sub             w5,  w5, w9          // i1
        add             w9,  w5, w7          // i0
        add             w7,  w5, w6          // i2
        str             w9,  [x0], #4
        str             w5,  [x1], #4
        str             w7,  [x2], #4
        subs            w3,  w3, #1
        b.ne            3b
8:      ret
endfunc

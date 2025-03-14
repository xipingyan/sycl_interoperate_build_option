
#ifndef GPU_INTEL_OCL_GENERIC_VECTOR_OPS_H
#define GPU_INTEL_OCL_GENERIC_VECTOR_OPS_H
typedef half __attribute__((ext_vector_type(1))) half1;
typedef uint __attribute__((ext_vector_type(1))) uint1;
typedef float __attribute__((ext_vector_type(1))) float1;
float1 __attribute__((overloadable)) vmad(float1 a, float1 b, float1 c) {
 c[0] = mad(a[0], b[0], c[0]);
 return c;
}
float2 __attribute__((overloadable)) vmad(float2 a, float2 b, float2 c) {
 return mad(a, b, c);
}
float4 __attribute__((overloadable)) vmad(float4 a, float4 b, float4 c) {
 return mad(a, b, c);
}
float8 __attribute__((overloadable)) vmad(float8 a, float8 b, float8 c) {
 return mad(a, b, c);
}
float16 __attribute__((overloadable)) vmad(float16 a, float16 b, float16 c) {
 return mad(a, b, c);
}
float1 __attribute__((overloadable)) native_vrecip(float1 x) {
 x[0] = native_recip(x[0]);
 return x;
}
float2 __attribute__((overloadable)) native_vrecip(float2 x) {
 return native_recip(x);
}
float4 __attribute__((overloadable)) native_vrecip(float4 x) {
 return native_recip(x);
}
float8 __attribute__((overloadable)) native_vrecip(float8 x) {
 return native_recip(x);
}
float16 __attribute__((overloadable)) native_vrecip(float16 x) {
 return native_recip(x);
}
float1 __attribute__((overloadable)) native_vexp2(float1 x) {
 x[0] = native_exp2(x[0]);
 return x;
}
float2 __attribute__((overloadable)) native_vexp2(float2 x) {
 return native_exp2(x);
}
float4 __attribute__((overloadable)) native_vexp2(float4 x) {
 return native_exp2(x);
}
float8 __attribute__((overloadable)) native_vexp2(float8 x) {
 return native_exp2(x);
}
float16 __attribute__((overloadable)) native_vexp2(float16 x) {
 return native_exp2(x);
}
#endif

#ifndef GPU_OCL_SDPA_UTILS_H
#define GPU_OCL_SDPA_UTILS_H
#define _4D_OFF(tag, x0, x1, x2, x3) (((x0) % tag##_B0) * tag##_SB0 + ((x0) / tag##_B0) * tag##_S0 + ((x1) % tag##_B1) * tag##_SB1 + ((x1) / tag##_B1) * tag##_S1 + ((x2) % tag##_B2) * tag##_SB2 + ((x2) / tag##_B2) * tag##_S2 + ((x3) % tag##_B3) * tag##_SB3 + ((x3) / tag##_B3) * tag##_S3)
#define QRY_OFF(x0, x1, x2, x3) _4D_OFF(QRY, x0, x1, x2, x3)
#define KEY_OFF(x0, x1, x2, x3) _4D_OFF(KEY, x0, x1, x2, x3)
#define VAL_OFF(x0, x1, x2, x3) _4D_OFF(VAL, x0, x1, x2, x3)
#define MSK_OFF(x0, x1, x2, x3) _4D_OFF(MSK, x0, x1, x2, x3)
#define DST_OFF(x0, x1, d, h, w) (((x0) % DST_B0) * DST_SB0 + ((x0) / DST_B0) * DST_S0 + ((x1) % DST_B1) * DST_SB1 + ((x1) / DST_B1) * DST_S1)
#endif

#ifndef GPU_OCL_TILE_OPS_H
#define GPU_OCL_TILE_OPS_H
float __builtin_IB_atomic_max_local_f32(__local float *, float);
__attribute__((overloadable)) float local_atomic_max(local float *p, float v) {
 return __builtin_IB_atomic_max_local_f32(p, v);
}
__attribute__((overloadable)) half local_atomic_max(
 local half *p, half v) {
 return v;
}
__attribute__((overloadable)) uint local_atomic_max(local uint *p, uint v) {
 return atomic_max(p, v);
}
__attribute__((overloadable)) int local_atomic_max(local int *p, int v) {
 return atomic_max(p, v);
}
#define DEF_BLOCK_LOAD_STORE(type, itype, suffix, n) __attribute__((overloadable)) type##n block_load( const global type *p, int vlen) __attribute__((enable_if(vlen == n, "wrong vector length"))) { return as_##type##n( intel_sub_group_block_read##suffix##n((global void *)p)); } __attribute__((overloadable)) void block_store( global type *p, type##n v) { intel_sub_group_block_write##suffix##n( (global itype *)p, as_##itype##n(v)); }
#define DEF_BLOCK_LOAD_STORE1(type, itype, suffix) __attribute__((overloadable)) type##1 block_load(const global type *p, int vlen) __attribute__( (enable_if(vlen == 1, "wrong vector length"))) { type##1 x; x[0] = as_##type( intel_sub_group_block_read##suffix((global void *)p)); return x; } __attribute__((overloadable)) void block_store( global type *p, type##1 v) { intel_sub_group_block_write##suffix( (global itype *)p, as_##itype(v[0])); }
#define DEF_BLOCK_LOAD_STORE16(type, itype, suffix) __attribute__((overloadable)) type##16 block_load(const global type *p, int vlen) __attribute__( (enable_if(vlen == 16, "wrong vector length"))) { type##16 x; x.s01234567 = as_##type##8( intel_sub_group_block_read8##suffix((global void *)p)); x.s89abcdef = as_##type##8( intel_sub_group_block_read8##suffix((global void *)(p + 8 * get_sub_group_size()))); return x; } __attribute__((overloadable)) void block_store( global type *p, type##16 v) { intel_sub_group_block_write8##suffix( (global itype *)p, as_##itype##8(v.s01234567)); intel_sub_group_block_write8##suffix( (global itype *)(p + 8 * get_sub_group_size()), as_##itype##8(v.s89abcdef)); }
DEF_BLOCK_LOAD_STORE1(half, ushort, _us)
DEF_BLOCK_LOAD_STORE(half, ushort, _us, 2)
DEF_BLOCK_LOAD_STORE(half, ushort, _us, 4)
DEF_BLOCK_LOAD_STORE(half, ushort, _us, 8)
DEF_BLOCK_LOAD_STORE(half, ushort, _us, 16)
DEF_BLOCK_LOAD_STORE1(uint, uint, )
DEF_BLOCK_LOAD_STORE(uint, uint, , 2)
DEF_BLOCK_LOAD_STORE(uint, uint, , 4)
DEF_BLOCK_LOAD_STORE(uint, uint, , 8)
DEF_BLOCK_LOAD_STORE16(uint, uint, )
#define DEF_BLOCK2D_LOAD_STORE(type, itype, vl, SG, suffix, BR, BC) itype##vl __builtin_IB_subgroup_block_read_flat_##suffix( long, int, int, int, int2); void __builtin_IB_subgroup_block_write_flat_##suffix( long, int, int, int, int2, itype##vl); __attribute__((overloadable)) type##vl block2d_load(const global type *p, int w, int h, int ld, int x, int y, int br, int bc, int sg) __attribute__((enable_if(br == BR, "wrong #rows"))) __attribute__((enable_if(bc == BC, "wrong #columns"))) __attribute__( (enable_if(sg == SG, "wrong subgroup size"))) { ulong pp = as_long(p); ulong prem = pp & 0x3F; pp &= ~0x3F; x += (prem / sizeof(type)); w += prem; int2 coord = {x, y}; return as_##type##vl(__builtin_IB_subgroup_block_read_flat_##suffix( pp, w - 1, h - 1, ld - 1, coord)); } __attribute__((overloadable)) void block2d_store(type##vl v, global type *p, int w, int h, int ld, int x, int y, int br, int bc, int sg) __attribute__((enable_if(br == BR, "wrong #rows"))) __attribute__((enable_if(bc == BC, "wrong #columns"))) __attribute__( (enable_if(sg == SG, "wrong subgroup size"))) { ulong pp = as_long(p); ulong prem = pp & 0x3F; pp &= ~0x3F; x += (prem / sizeof(type)); w += prem; int2 coord = {x, y}; __builtin_IB_subgroup_block_write_flat_##suffix( pp, w - 1, h - 1, ld - 1, coord, as_##itype##vl(v)); }
DEF_BLOCK2D_LOAD_STORE(half, ushort, 8, 16, u16_m8k16v1, 16, 8)
DEF_BLOCK2D_LOAD_STORE(half, ushort, 8, 16, u16_m4k32v1, 32, 4)
DEF_BLOCK2D_LOAD_STORE(half, ushort, 16, 16, u16_m8k32v1, 32, 8)
#define tile_fill(t, v) do { _Pragma("unroll") for (int i = 0; i < sizeof(t.x) / sizeof(t.x[0]); i++) t.x[i] = v; } while (0)
#define tile_elementwise(t, f) do { _Pragma("unroll") for (int i = 0; i < sizeof(t.x) / sizeof(t.x[0]); i++) t.x[i] = f(t.x[i]); } while (0)
#define tile_elementwise_s(t, f) do { _Pragma("unroll") for (int i = 0; i < sizeof(t.x) / sizeof(t.x[0]); i++) { _Pragma("unroll") for (int s = 0; s < sizeof(t.x[0]) / sizeof(t.x[0][0]); s++) t.x[i][s] = f(t.x[i][s]); } } while (0)
#define tile_binary(t, t2, f) do { _Pragma("unroll") for (int i = 0; i < sizeof(t.x) / sizeof(t.x[0]); i++) t.x[i] = f(t.x[i], t2.x[i]); } while (0)
#define tile_copy(t, t_new) do { _Pragma("unroll") for (int i = 0; i < sizeof(t.x) / sizeof(t.x[0]); i++) t_new.x[i] = __builtin_convertvector(t.x[i], __typeof__(t_new.x[i])); } while (0)
#define tile_copy_to_half2(t, t_new) do { _Pragma("unroll") for (int i = 0; i < sizeof(t.x) / sizeof(t.x[0]); i++) { _Pragma("unroll") for (int s = 0; s < sizeof(t.x[0]) / sizeof(t.x[0][0]) / 2; s++) { half2 v = {t.x[i][2 * s], t.x[i][2 * s + 1]}; t_new.x[i][s] = as_uint(v); } } } while (0)
#define tile_access(t, i0, j, sg, br, bc, nbr) (t).x[(i0) / (br) + (nbr) * ((j) / (bc))] [((i0) % (br)) / (sg) + ((j) % (bc)) * ((br) / (sg))]
#define xlane_tile_access(t, i, j, sg, br, bc, nbr) sub_group_broadcast(tile_access(t, i, j, sg, br, bc, nbr), i % sg)
#define DECLARE_2D_TILE_OPS(tile_type, element_type, sg, br, bc, nbr, nbc) __attribute__((overloadable)) void tile_load_full(tile_type *t, const global element_type *ptr, int ld, int offset_r,  int offset_c) { ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int j = 0; j < bc * nbc; j++, ptr += ld) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { int i = i0 + get_sub_group_local_id(); tile_access(*t, i0, j, sg, br, bc, nbr) = ptr[i]; } } } __attribute__((overloadable)) void tile_load_full(tile_type *t, const local element_type *ptr, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int j = 0; j < bc * nbc; j++, ptr += ld) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { int i = i0 + get_sub_group_local_id(); tile_access(*t, i0, j, sg, br, bc, nbr) = ptr[i]; } } } __attribute__((overloadable)) void tile_load(tile_type *t, const global element_type *ptr, int m, int n, int ld, int offset_r, int offset_c) { if (m >= offset_r + br * nbr && n >= offset_c + bc * nbc) { tile_load_full(t, ptr, ld, offset_r, offset_c); return; } ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int j = 0; j < bc * nbc; j++, ptr += ld) { if (offset_c + j < n) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { int i = i0 + get_sub_group_local_id(); if (offset_r + i < m) tile_access(*t, i0, j, sg, br, bc, nbr) = ptr[i]; } } } } __attribute__((overloadable)) void tile_load(tile_type *t, const global element_type *ptr, int m, int n, int offset_r, int offset_c) { tile_load(t, ptr, m, n, m, offset_r, offset_c); } __attribute__((overloadable)) void tile_load_t_full(tile_type *t, const global element_type *ptr, int ld, int offset_r, int offset_c) { ptr += ld * offset_r + offset_c; _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg, ptr += ld*sg) { _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { tile_access(*t, i0, j, sg, br, bc, nbr) = ptr[get_sub_group_local_id() * ld + j]; } } } __attribute__((overloadable)) void tile_load_t(tile_type *t, const global element_type *ptr, int m, int n, int ld, int offset_r, int offset_c) { if (m >= offset_r + br * nbr && n >= offset_c + bc * nbc) { tile_load_t_full(t, ptr, ld, offset_r, offset_c); return; } ptr += ld * offset_r + offset_c; _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg, ptr += ld*sg) { int i = i0 + get_sub_group_local_id(); if (offset_r + i < m) _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { if (offset_c + j < n) { tile_access(*t, i0, j, sg, br, bc, nbr) = ptr[get_sub_group_local_id() * ld + j]; } } } } __attribute__((overloadable)) void tile_load_t(tile_type *t, const global element_type *ptr, int m, int n, int offset_r, int offset_c) { tile_load(t, ptr, m, n, n, offset_r, offset_c); } __attribute__((overloadable)) void tile_store_full(tile_type t, local element_type *ptr, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int j = 0; j < bc * nbc; j++, ptr += ld) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { int i = i0 + get_sub_group_local_id(); ptr[i] = tile_access(t, i0, j, sg, br, bc, nbr); } } } __attribute__((overloadable)) void tile_store_full(tile_type t, global element_type *ptr, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int j = 0; j < bc * nbc; j++, ptr += ld) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { int i = i0 + get_sub_group_local_id(); ptr[i] = tile_access(t, i0, j, sg, br, bc, nbr); } } } __attribute__((overloadable)) void tile_store(tile_type t, global element_type *ptr, int m, int n, int ld, int offset_r, int offset_c) { if (m >= offset_r + br * nbr && n >= offset_c + bc * nbc) { tile_store_full(t, ptr, ld, offset_r, offset_c); return; } ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int j = 0; j < bc * nbc; j++, ptr += ld) { if (offset_c + j < n) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { int i = i0 + get_sub_group_local_id(); if (offset_r + i < m) ptr[i] = tile_access(t, i0, j, sg, br, bc, nbr); } } } } __attribute__((overloadable)) void tile_store(tile_type t, global element_type *ptr, int m, int n, int offset_r, int offset_c) { tile_store(t, ptr, m, n, m, offset_r, offset_c); } __attribute__((overloadable)) void tile_store_t_sys_src1(tile_type t, local element_type *ptr, int ld, int offset_r, int offset_c) { offset_c += get_sub_group_local_id(); int offset_r0 = offset_r & (sg - 1); int offset_r1 = offset_r & ~(sg - 1); ptr += offset_r0 + sg * offset_c + ld * offset_r1; _Pragma("unroll") for (int j0 = 0; j0 < br * nbr; j0 += sg, ptr += sg * sg) { _Pragma("unroll") for (int i = 0; i < bc * nbc; i++) ptr[i] = tile_access(t, j0, i, sg, br, bc, nbr); } } __attribute__((overloadable)) void tile_store_t_sys_src2(tile_type t, local element_type *ptr, int tile_n, int ld, int offset_r, int offset_c) { const int cp = 32 / sizeof(element_type); offset_c += get_sub_group_local_id(); int offset_r0 = offset_r & (cp - 1); int offset_r1 = offset_r & ~(cp - 1); ptr += offset_r0 + tile_n * offset_r1; _Pragma("unroll") for (int j0 = 0; j0 < br * nbr; j0 += sg, offset_c += sg) { int offset_c0 = offset_c & (tile_n - 1); int offset_c1 = offset_c & ~(tile_n - 1); local element_type *ptr_j = ptr + cp * offset_c0 + ld * offset_c1; _Pragma("unroll") for (int i = 0; i < bc * nbc; i++) { *ptr_j = tile_access(t, j0, i, sg, br, bc, nbr); ptr_j++; if ((~i & (cp - 1)) == 0) ptr_j += cp * (tile_n - 1); } } } __attribute__((overloadable)) void tile_atomic_max_full(tile_type t, local element_type *ptr, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int j = 0; j < bc * nbc; j++, ptr += ld) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { int i = i0 + get_sub_group_local_id(); (void)local_atomic_max( ptr + i, tile_access(t, i0, j, sg, br, bc, nbr)); } } }
#define DECLARE_2D_TILE_VREDUCE(tile_type, sg, br, bc, nbr, nbc, rtile_type, rsg, rbr, rbc, rnbr, rnbc) __attribute__((overloadable)) void tile_vreduce_add( tile_type t, rtile_type *tr) { _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { tile_access(*tr, i0, 0, rsg, rbr, rbc, rnbr) += tile_access(t, i0, j, sg, br, bc, nbr); } } } __attribute__((overloadable)) void tile_vreduce_max( tile_type t, rtile_type *tr) { _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { tile_access(*tr, i0, 0, rsg, rbr, rbc, rnbr) = max(tile_access(t, i0, j, sg, br, bc, nbr), tile_access(*tr, i0, 0, rsg, rbr, rbc, rnbr)); } } } __attribute__((overloadable)) void tile_vbroadcast_sub( tile_type *t, rtile_type tr) { _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { tile_access(*t, i0, j, sg, br, bc, nbr) -= tile_access(tr, i0, 0, rsg, rbr, rbc, rnbr); } } } __attribute__((overloadable)) void tile_vbroadcast_mul( tile_type *t, rtile_type tr) { _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { tile_access(*t, i0, j, sg, br, bc, nbr) *= tile_access(tr, i0, 0, rsg, rbr, rbc, rnbr); }  } }
#define DECLARE_2D_TILE_HREDUCE(tile_type, sg, br, bc, nbr, nbc, rtile_type, rsg, rbr, rbc, rnbr, rnbc) __attribute__((overloadable)) void tile_hbroadcast_add( tile_type *t, rtile_type tr) { _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { tile_access(*t, i0, j, sg, br, bc, nbr) += xlane_tile_access(tr, j, 0, rsg, rbr, rbc, rnbr); } } } __attribute__((overloadable)) void tile_hbroadcast_mul( tile_type *t, rtile_type tr) { _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { tile_access(*t, i0, j, sg, br, bc, nbr) *= xlane_tile_access(tr, j, 0, rsg, rbr, rbc, rnbr); } } } __attribute__((overloadable)) void tile_hbroadcast_min( tile_type *t, rtile_type tr) { _Pragma("unroll") for (int j = 0; j < bc * nbc; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { tile_access(*t, i0, j, sg, br, bc, nbr) = min( tile_access(*t, i0, j, sg, br, bc, nbr), xlane_tile_access(tr, j, 0, rsg, rbr, rbc, rnbr)); } } }
#define DECLARE_2D_TILE_RSELECT(tile_type0, sg0, br0, bc0, nbr0, nbc0, tile_type1, sg1, br1, bc1, nbr1, nbc1) __attribute__((overloadable)) void tile_rselect( tile_type0 *t0, tile_type1 t1, int idx) { _Pragma("unroll") for (int j = 0; j < bc0 * nbc0; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br0 * nbr0; i0 += sg0) { tile_access(*t0, i0, j, sg0, br0, bc0, nbr0) = tile_access(t1, i0, j, sg1, br1, bc1, nbr1); _Pragma("unroll") for (int z = 1; z < (br1 * nbr1 / br0 * nbr0); z++) if (z == idx) { tile_access(*t0, i0, j, sg0, br0, bc0, nbr0) = tile_access(t1, i0 + z * br0 * nbr0, j, sg1, br1, bc1, nbr1); } } } }
#define DECLARE_2D_TILE_COPY_REBLOCK(tile_type0, sg0, br0, bc0, nbr0, nbc0, tile_type1, sg1, br1, bc1, nbr1, nbc1) __attribute__((overloadable)) void tile_copy_reblock( tile_type0 t0, tile_type1 *t1) { _Pragma("unroll") for (int j = 0; j < bc0 * nbc0; j++) { _Pragma("unroll") for (int i0 = 0; i0 < br0 * nbr0; i0 += sg0) { tile_access(*t1, i0, j, sg1, br1, bc1, nbr1) = tile_access(t0, i0, j, sg0, br0, bc0, nbr0); } } }
#define DECLARE_2D_TILE(tile_type, element_type, sg, br, bc, nbr, nbc) typedef element_type __attribute__((ext_vector_type(br * bc / sg))) _e_##tile_type; typedef struct { _e_##tile_type x[nbr * nbc]; } tile_type; DECLARE_2D_TILE_OPS(tile_type, element_type, sg, br, bc, nbr, nbc)
#define DECLARE_2D_TILE_BLOCK_OPS( tile_type, element_type, sg, br, bc, nbr, nbc) __attribute__((overloadable)) void tile_load_block(tile_type *t, const global element_type *ptr, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int jj = 0; jj < nbc; jj++, ptr += ld * bc) { _Pragma("unroll") for (int ii = 0; ii < nbr; ii++)(t) ->x[ii + nbr * jj] = block_load(ptr + ii * br, br / SUBGROUP_SIZE); } } __attribute__((overloadable)) void tile_store_block(tile_type t, global element_type *ptr, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int jj = 0; jj < nbc; jj++, ptr += ld * bc) { _Pragma("unroll") for (int ii = 0; ii < nbr; ii++) block_store(ptr + ii * br, (t).x[ii + nbr * jj]); } } __attribute__((overloadable)) void tile_load_block(tile_type *t, const global element_type *ptr, int n, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; n -= offset_c; _Pragma("unroll") for (int jj = 0; jj < nbc; jj++, ptr += ld * bc) { if (jj < n) { _Pragma("unroll") for (int ii = 0; ii < nbr; ii++)(t) ->x[ii + nbr * jj] = block_load(ptr + ii * br, br / SUBGROUP_SIZE); } } } __attribute__((overloadable)) void tile_store_block(tile_type t, global element_type *ptr, int n, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; n -= offset_c; _Pragma("unroll") for (int jj = 0; jj < nbc; jj++, ptr += ld * bc) { if (jj < n) { _Pragma("unroll") for (int ii = 0; ii < nbr; ii++) block_store(ptr + ii * br, (t).x[ii + nbr * jj]); } } }
#define DECLARE_2D_TILE_BLOCK2D_OPS( tile_type, element_type, sg, br, bc, nbr, nbc) __attribute__((overloadable)) void tile_load_block2d(tile_type *t, const global element_type *ptr, int m, int n, int ld, int offset_r, int offset_c) { const int e = sizeof(element_type); _Pragma("unroll") for (int jj = 0; jj < nbc; jj++) { _Pragma("unroll") for (int ii = 0; ii < nbr; ii++)(t) ->x[ii + nbr * jj] = block2d_load(ptr, m * e, n, ld * e, offset_r + ii * br, offset_c + jj * bc, br, bc, sg); } } __attribute__((overloadable)) void tile_load_block2d(tile_type *t, const global element_type *ptr, int m, int n, int offset_r, int offset_c) { tile_load_block2d(t, ptr, m, n, m, offset_r, offset_c); } __attribute__((overloadable)) void tile_store_block2d(tile_type t, global element_type *ptr, int m, int n, int ld, int offset_r, int offset_c) { const int e = sizeof(element_type); _Pragma("unroll") for (int jj = 0; jj < nbc; jj++) { _Pragma("unroll") for (int ii = 0; ii < nbr; ii++) block2d_store( (t).x[ii + nbr * jj], ptr, m *e, n, ld *e, offset_r + ii * br, offset_c + jj * bc, br, bc, sg); } } __attribute__((overloadable)) void tile_store_block2d(tile_type t, const global element_type *ptr, int m, int n, int offset_r, int offset_c) { tile_store_block2d(t, ptr, m, n, m, offset_r, offset_c); }
#define DECLARE_2D_TILE_LOAD_PACKED_HALF(tile_type, sg, br, bc, nbr, nbc) __attribute__((overloadable)) void tile_load_packed_half(tile_type *t, const global half *ptr, int m, int n, int ld, int offset_r, int offset_c) { ptr += ld * offset_c + offset_r; _Pragma("unroll") for (int j = 0; j < bc * nbc; j++, ptr += ld) { if (offset_c + j < n) { _Pragma("unroll") for (int i0 = 0; i0 < br * nbr; i0 += sg) { int i = 2 * (i0 + get_sub_group_local_id()); half2 loaded = 0; if (offset_r + i < m) loaded.s0 = ptr[i]; if (offset_r + i + 1 < m) loaded.s1 = ptr[i + 1]; tile_access(*t, i0, j, sg, br, bc, nbr) = as_uint(loaded); } } } } __attribute__((overloadable)) void tile_load_packed_half(tile_type *t, const global half *ptr, int m, int n, int offset_r, int offset_c) { tile_load_packed_half(t, ptr, m, n, m, offset_r, offset_c); }
#define cooperative_prefetch_2d(ptr, r, c, ld, sg_id, n_sg, sg_size, caching) cooperative_prefetch_2d_internal((const global char *)ptr, (r) * sizeof(*(ptr)), c, (ld) * sizeof(*(ptr)), sg_id, n_sg, sg_size, caching)
#define cooperative_prefetch_2d_rem( ptr, r, c, rmax, cmax, ld, sg_id, n_sg, sg_size, caching) cooperative_prefetch_2d_internal((const global char *)ptr, (r) * sizeof(*(ptr)), c, (rmax) * sizeof(*(ptr)), cmax, (ld) * sizeof(*(ptr)), sg_id, n_sg, sg_size, caching)
enum LSC_LDCC {
 LSC_LDCC_DEFAULT = 0,
 LSC_LDCC_L1UC_L3UC = 1,
 LSC_LDCC_L1UC_L3C = 2,
 LSC_LDCC_L1C_L3UC = 3,
 LSC_LDCC_L1C_L3C = 4,
 LSC_LDCC_L1S_L3UC = 5,
 LSC_LDCC_L1S_L3C = 6,
 LSC_LDCC_L1IAR_L3C = 7,
};
extern void __builtin_IB_lsc_prefetch_global_uchar(
 const __global uchar *base, int immElemOff, enum LSC_LDCC cacheOpt);
extern void __builtin_IB_lsc_prefetch_global_uint(
 const __global uint *base, int immElemOff, enum LSC_LDCC cacheOpt);
__attribute__((overloadable)) void cooperative_prefetch_2d_internal(
 const global char *ptr, uint rbytes, uint c, uint ld_bytes, uint sg_id,
 uint n_sg, uint sg_size, enum LSC_LDCC caching) {
 const uint cl_per_col = (rbytes + 63) >> 6;
 const uint cl = cl_per_col * c;
 const uint cl_per_sg = (cl + n_sg - 1) / n_sg;
 const uint cl_iters = (cl_per_sg + sg_size - 1) / sg_size;
#pragma unroll
 for (uint ii_cl = 0; ii_cl < cl_iters; ii_cl++) {
 uint i_cl = ii_cl + (sg_id * cl_per_sg) + get_sub_group_local_id();
 uint r_cl = i_cl % cl_per_col;
 uint c_cl = i_cl / cl_per_col;
 if (i_cl < cl) {
 __builtin_IB_lsc_prefetch_global_uint(
 (const global uint *)(ptr + r_cl * 64 + c_cl * ld_bytes), 0,
 caching);
 }
 }
}
__attribute__((overloadable)) void cooperative_prefetch_2d_internal(
 const global char *ptr, uint rbytes, uint c, uint rbytes_max,
 uint c_max, uint ld_bytes, uint sg_id, uint n_sg, uint sg_size,
 enum LSC_LDCC caching) {
 const uint cl_per_col = (rbytes_max + 63) >> 6;
 const uint cl = cl_per_col * c_max;
 const uint cl_per_sg = (cl + n_sg - 1) / n_sg;
 const uint cl_iters = (cl_per_sg + sg_size - 1) / sg_size;
 const uint max_off = rbytes - 1 + (c - 1) * ld_bytes;
#pragma unroll
 for (uint ii_cl = 0; ii_cl < cl_iters; ii_cl++) {
 uint i_cl = ii_cl + (sg_id * cl_per_sg) + get_sub_group_local_id();
 uint r_cl = i_cl % cl_per_col;
 uint c_cl = i_cl / cl_per_col;
 uint pf_off = min(r_cl * 64 + c_cl * ld_bytes, max_off);
 if (i_cl < cl) {
 __builtin_IB_lsc_prefetch_global_uchar(
 (const global uchar *)(ptr + pf_off), 0, caching);
 }
 }
}
#endif

#ifdef intel_convert_as_bfloat16_float
#define _convert_as_bfloat16_float(val) intel_convert_as_bfloat16_float(val)
#else
inline float _convert_as_bfloat16_float(ushort source) {
 uint u = 0;
 if ( (source>>15) ) {
 u = 1 << 31;
 }
 u += ( ( (source >> 7) & 0b11111111)) << 23;
 u += (source & 0b1111111) << 16;
 float* f = &u;
 return *f;
}
#endif
#ifdef intel_convert_bfloat16_as_ushort
#define _convert_bfloat16_as_ushort(val) intel_convert_bfloat16_as_ushort(val)
#else
inline ushort _convert_bfloat16_as_ushort(float source) {
 uint* in = &source;
 ushort u = 0;
 if ( (*in>>31) ) {
 u = 1 << 15;
 }
 u += ( ( (*in >> 23) & 0b11111111)) << 7;
 u += (*in >> 16) & 0b1111111;
 return u;
}
#endif

#if defined(cl_khr_fp16)
#pragma OPENCL EXTENSION cl_khr_fp16 : enable
#endif
#if !defined(cl_intel_subgroups) && defined(cl_khr_subgroups)
#pragma OPENCL EXTENSION cl_khr_subgroups : enable
#endif
#define __CAT(x, y) x##y
#define CAT(x, y) __CAT(x, y)
#define OFFSET_GLOBAL_PTR(elem_type, ptr, byte_offset) ((__global elem_type*)((__global char*)(ptr) + (byte_offset)))
#define MULTIPLY_OFFSET(elem_type, byte_offset) ((byte_offset) * sizeof(elem_type))
#if OPT_HINTS_SUPPORTED
# define ASSUME_HINT(x) __builtin_assume(x)
#else
# define ASSUME_HINT(x) do { } while (0)
#endif
#define unroll_for __attribute__((opencl_unroll_hint)) for
#define CEIL_DIV(a, b) (((a) + (b) - 1)/(b))
#define ALIGN(a, b) (CEIL_DIV(a, b) * (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define CLAMP(v,l,u) MAX((l),MIN((v),(u)))
#define MAKE_VECTOR_TYPE_IMPL_1(elem_type) elem_type
#define MAKE_VECTOR_TYPE_IMPL_2(elem_type) CAT(elem_type, 2)
#define MAKE_VECTOR_TYPE_IMPL_3(elem_type) CAT(elem_type, 3)
#define MAKE_VECTOR_TYPE_IMPL_4(elem_type) CAT(elem_type, 4)
#define MAKE_VECTOR_TYPE_IMPL_8(elem_type) CAT(elem_type, 8)
#define MAKE_VECTOR_TYPE_IMPL_16(elem_type) CAT(elem_type, 16)
#define MAKE_VECTOR_TYPE(elem_type, size) CAT(MAKE_VECTOR_TYPE_IMPL_, size)(elem_type)
#define AS_TYPE(type, val) CAT(as_, type)(val)
#define TYPE_SIZE_uchar 1
#define TYPE_SIZE_char 1
#define TYPE_SIZE_ushort 2
#define TYPE_SIZE_short 2
#define TYPE_SIZE_half 2
#define TYPE_SIZE_int 4
#define TYPE_SIZE_uint 4
#define TYPE_SIZE_float 4
#define TYPE_SIZE_ulong 8
#define TYPE_SIZE_long 8
#define TYPE_SIZE(type) CAT(TYPE_SIZE_, type)
#ifdef cl_intel_required_subgroup_size
#define REQD_SUB_GROUP_SIZE(sg_size) __attribute__((intel_reqd_sub_group_size(sg_size)))
#else
#define REQD_SUB_GROUP_SIZE(sg_size)
#endif

#define GET_DATA_INDEX(prefix, b, f, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (f)*CAT(prefix, _FEATURE_PITCH) + (b)*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_RAW(prefix, i0, i1, i2, i3) CAT(prefix, _OFFSET) + (i0)*CAT(prefix, _PITCHES)[0] + (i1)*CAT(prefix, _PITCHES)[1] + (i2)*CAT(prefix, _PITCHES)[2] + (i3)*CAT(prefix, _PITCHES)[3]
#define GET_DATA_INDEX_SAFE(prefix, b, f, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (f % CAT(prefix, _FEATURE_NUM))*CAT(prefix, _FEATURE_PITCH) + (b % CAT(prefix, _BATCH_NUM ))*CAT(prefix, _BATCH_PITCH)
 #define GET_DATA_INDEX_5D(prefix, b, f, z, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (z)*CAT(prefix, _Z_PITCH) + (f)*CAT(prefix, _FEATURE_PITCH) + (b)*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_5D_RAW(prefix, i0, i1, i2, i3, i4) CAT(prefix, _OFFSET) + (i0)*CAT(prefix, _PITCHES)[0] + (i1)*CAT(prefix, _PITCHES)[1] + (i2)*CAT(prefix, _PITCHES)[2] + (i3)*CAT(prefix, _PITCHES)[3] + (i4)*CAT(prefix, _PITCHES)[4]
#define GET_DATA_INDEX_5D_SAFE(prefix, b, f, z, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (z % CAT(prefix, _SIZE_Z ))*CAT(prefix, _Z_PITCH) + (f % CAT(prefix, _FEATURE_NUM))*CAT(prefix, _FEATURE_PITCH) + (b % CAT(prefix, _BATCH_NUM ))*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_6D(prefix, b, f, w, z, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (z)*CAT(prefix, _Z_PITCH) + (w)*CAT(prefix, _W_PITCH) + (f)*CAT(prefix, _FEATURE_PITCH) + (b)*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_6D_SAFE(prefix, b, f, w, z, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (z % CAT(prefix, _SIZE_Z ))*CAT(prefix, _Z_PITCH) + (w % CAT(prefix, _SIZE_W ))*CAT(prefix, _W_PITCH) + (f % CAT(prefix, _FEATURE_NUM))*CAT(prefix, _FEATURE_PITCH) + (b % CAT(prefix, _BATCH_NUM ))*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_6D_RAW(prefix, i0, i1, i2, i3, i4, i5) CAT(prefix, _OFFSET) + (i0)*CAT(prefix, _PITCHES)[0] + (i1)*CAT(prefix, _PITCHES)[1] + (i2)*CAT(prefix, _PITCHES)[2] + (i3)*CAT(prefix, _PITCHES)[3] + (i4)*CAT(prefix, _PITCHES)[4] + (i5)*CAT(prefix, _PITCHES)[5]
#define GET_DATA_INDEX_7D(prefix, b, f, u, w, z, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (z)*CAT(prefix, _Z_PITCH) + (w)*CAT(prefix, _W_PITCH) + (u)*CAT(prefix, _U_PITCH) + (f)*CAT(prefix, _FEATURE_PITCH) + (b)*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_7D_SAFE(prefix, b, f, u, w, z, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (z % CAT(prefix, _SIZE_Z ))*CAT(prefix, _Z_PITCH) + (w % CAT(prefix, _SIZE_W ))*CAT(prefix, _W_PITCH) + (u % CAT(prefix, _SIZE_U ))*CAT(prefix, _U_PITCH) + (f % CAT(prefix, _FEATURE_NUM))*CAT(prefix, _FEATURE_PITCH) + (b % CAT(prefix, _BATCH_NUM ))*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_7D_RAW(prefix, i0, i1, i2, i3, i4, i5, i6) CAT(prefix, _OFFSET) + (i0)*CAT(prefix, _PITCHES)[0] + (i1)*CAT(prefix, _PITCHES)[1] + (i2)*CAT(prefix, _PITCHES)[2] + (i3)*CAT(prefix, _PITCHES)[3] + (i4)*CAT(prefix, _PITCHES)[4] + (i5)*CAT(prefix, _PITCHES)[5] + (i6)*CAT(prefix, _PITCHES)[6]
#define GET_DATA_INDEX_8D(prefix, b, f, v, u, w, z, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (z)*CAT(prefix, _Z_PITCH) + (w)*CAT(prefix, _W_PITCH) + (u)*CAT(prefix, _U_PITCH) + (v)*CAT(prefix, _V_PITCH) + (f)*CAT(prefix, _FEATURE_PITCH) + (b)*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_8D_SAFE(prefix, b, f, v, u, w, z, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (z % CAT(prefix, _SIZE_Z ))*CAT(prefix, _Z_PITCH) + (w % CAT(prefix, _SIZE_W ))*CAT(prefix, _W_PITCH) + (u % CAT(prefix, _SIZE_U ))*CAT(prefix, _U_PITCH) + (v % CAT(prefix, _SIZE_V ))*CAT(prefix, _V_PITCH) + (f % CAT(prefix, _FEATURE_NUM))*CAT(prefix, _FEATURE_PITCH) + (b % CAT(prefix, _BATCH_NUM ))*CAT(prefix, _BATCH_PITCH)
#define GET_DATA_INDEX_8D_RAW(prefix, i0, i1, i2, i3, i4, i5, i6, i7) CAT(prefix, _OFFSET) + (i0)*CAT(prefix, _PITCHES)[0] + (i1)*CAT(prefix, _PITCHES)[1] + (i2)*CAT(prefix, _PITCHES)[2] + (i3)*CAT(prefix, _PITCHES)[3] + (i4)*CAT(prefix, _PITCHES)[4] + (i5)*CAT(prefix, _PITCHES)[5] + (i6)*CAT(prefix, _PITCHES)[6] + (i7)*CAT(prefix, _PITCHES)[7]
#define GET_DATA_BS_FYX_BSV8_INDEX(prefix, b, f, y, x, sub_group_size) CAT(prefix, _OFFSET) + ((b) % (sub_group_size)) + (sub_group_size)*( (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (f)*CAT(prefix, _FEATURE_PITCH) + ((b) / (sub_group_size))*CAT(prefix, _BATCH_PITCH) )
inline uint get_b_fs_yx_fsv_index(uint b, uint f, uint y, uint x,
 uint x_size, uint y_size, uint f_size, uint b_size,
 uint b_pad_before, uint b_pad_after,
 uint f_pad_before, uint f_pad_after,
 uint y_pad_before, uint y_pad_after,
 uint x_pad_before, uint x_pad_after, uint alignment) {
 const uint feature = f + f_pad_before;
 const uint fs = feature / alignment;
 const uint fsv = feature % alignment;
 const uint x_pitch = alignment;
 const uint y_pitch = x_pitch * (x_pad_before + x_size + x_pad_after);
 const uint total_f_size = f_pad_before + f_size + f_pad_after;
 const uint fs_pitch = y_pitch * (y_pad_before + y_size + y_pad_after);
 const uint b_pitch = fs_pitch * ((total_f_size + alignment - 1) / alignment);
 const uint output_offset = (b_pad_before + b) * b_pitch +
 fs * fs_pitch +
 (y_pad_before + y) * y_pitch +
 (x_pad_before + x) * x_pitch
 + fsv;
 return output_offset;
}
inline uint get_b_fs_yx_fsv_index_safe(uint b, uint f, uint y, uint x,
 uint x_size, uint y_size, uint f_size, uint b_size,
 uint b_pad_before, uint b_pad_after,
 uint f_pad_before, uint f_pad_after,
 uint y_pad_before, uint y_pad_after,
 uint x_pad_before, uint x_pad_after, uint alignment) {
 const uint f_mod = f_pad_before + (f % f_size);
 const uint fs = f_mod / alignment;
 const uint fsv = f_mod % alignment;
 const uint x_pitch = alignment;
 const uint y_pitch = x_pitch * (x_pad_before + x_size + x_pad_after);
 const uint total_f_size = f_pad_before + f_size + f_pad_after;
 const uint fs_pitch = y_pitch * (y_pad_before + y_size + y_pad_after);
 const uint b_pitch = fs_pitch * ((total_f_size + alignment - 1) / alignment);
 const uint output_offset = (b_pad_before + (b % b_size)) * b_pitch +
 fs * fs_pitch +
 (y_pad_before + (y % y_size)) * y_pitch +
 (x_pad_before + (x % x_size)) * x_pitch
 + fsv;
 return output_offset;
}
#define GET_DATA_B_FS_YX_FSV16_INDEX(prefix, b, f, y, x) get_b_fs_yx_fsv_index( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16)
#define GET_DATA_B_FS_YX_FSV16_INDEX_SAFE(prefix, b, f, y, x) get_b_fs_yx_fsv_index_safe( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16)
#define GET_DATA_B_FS_YX_FSV2_INDEX(prefix, b, f, y, x) get_b_fs_yx_fsv_index( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 2)
#define GET_DATA_B_FS_YX_FSV2_INDEX_SAFE(prefix, b, f, y, x) get_b_fs_yx_fsv_index_safe( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 2)
#define GET_DATA_B_FS_YX_FSV4_INDEX(prefix, b, f, y, x) get_b_fs_yx_fsv_index( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 4)
#define GET_DATA_B_FS_YX_FSV4_INDEX_SAFE(prefix, b, f, y, x) get_b_fs_yx_fsv_index_safe( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 4)
#define GET_DATA_B_FS_YX_FSV8_INDEX(prefix, b, f, y, x) get_b_fs_yx_fsv_index( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8)
#define GET_DATA_B_FS_YX_FSV8_INDEX_SAFE(prefix, b, f, y, x) get_b_fs_yx_fsv_index_safe( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8)
#define GET_DATA_B_FS_YX_FSV32_INDEX(prefix, b, f, y, x) get_b_fs_yx_fsv_index( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32)
#define GET_DATA_B_FS_YX_FSV32_INDEX_SAFE(prefix, b, f, y, x) get_b_fs_yx_fsv_index_safe( b, f, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32)
inline uint get_fs_b_yx_fsv32_index(uint b, uint f, uint y, uint x,
 uint x_pad_before, uint x_size, uint x_pad_after,
 uint y_pad_before, uint y_size, uint y_pad_after,
 uint f_pad_before,
 uint size_b)
{
 const uint feature_tile_size = 32;
 const uint x_total_size = x_pad_before + x_size + x_pad_after;
 const uint y_total_size = y_pad_before + y_size + y_pad_after;
 const uint real_x = x + x_pad_before;
 const uint real_y = y + y_pad_before;
 const uint real_f = f + f_pad_before;
 const uint x_pitch = feature_tile_size;
 const uint y_pitch = x_pitch * x_total_size;
 const uint b_pitch = y_pitch * y_total_size;
 const uint f_tile_pitch = b_pitch * size_b;
 const uint feature_tile_number = real_f / feature_tile_size;
 const uint feature_local_number = real_f % feature_tile_size;
 size_t index = 0;
 index += feature_tile_number * f_tile_pitch;
 index += b * b_pitch;
 index += real_y * y_pitch;
 index += real_x * x_pitch;
 index += feature_local_number;
 return index;
}
inline uint get_fs_b_yx_fsv32_index_safe(uint b, uint f, uint y, uint x,
 uint x_pad_before, uint x_size, uint x_pad_after,
 uint y_pad_before, uint y_size, uint y_pad_after,
 uint f_pad_before, uint f_size,
 uint size_b)
{
 const uint feature_tile_size = 32;
 const uint x_total_size = x_pad_before + x_size + x_pad_after;
 const uint y_total_size = y_pad_before + y_size + y_pad_after;
 const uint real_x = (x % x_size) + x_pad_before;
 const uint real_y = (y % y_size) + y_pad_before;
 const uint real_f = (f % f_size) + f_pad_before;
 const uint x_pitch = feature_tile_size;
 const uint y_pitch = x_pitch * x_total_size;
 const uint b_pitch = y_pitch * y_total_size;
 const uint f_tile_pitch = b_pitch * size_b;
 const uint feature_tile_number = real_f / feature_tile_size;
 const uint feature_local_number = real_f % feature_tile_size;
 size_t index = 0;
 index += feature_tile_number * f_tile_pitch;
 index += b * b_pitch;
 index += real_y * y_pitch;
 index += real_x * x_pitch;
 index += feature_local_number;
 return index;
}
#define GET_DATA_FS_B_YX_FSV32_INDEX(prefix, b, f, y, x) get_fs_b_yx_fsv32_index( b, f, y, x, CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _BATCH_NUM))
#define GET_DATA_FS_B_YX_FSV32_INDEX_SAFE(prefix, b, f, y, x) get_fs_b_yx_fsv32_index_safe( b, f, y, x, CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM))
#define GET_DATA_B_FS_ZYX_FSV2_INDEX(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 2)
#define GET_DATA_B_FS_ZYX_FSV2_INDEX_SAFE(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 2)
#define GET_DATA_B_FS_ZYX_FSV4_INDEX(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 4)
#define GET_DATA_B_FS_ZYX_FSV4_INDEX_SAFE(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 4)
#define GET_DATA_B_FS_ZYX_FSV8_INDEX(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8)
#define GET_DATA_B_FS_ZYX_FSV8_INDEX_SAFE(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8)
#define GET_DATA_B_FS_ZYX_FSV16_INDEX(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16)
#define GET_DATA_B_FS_ZYX_FSV16_INDEX_SAFE(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y),  CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16)
#define GET_DATA_B_FS_ZYX_FSV32_INDEX(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32)
#define GET_DATA_B_FS_ZYX_FSV32_INDEX_SAFE(prefix, b, f, z, y, x) get_b_fs_zyx_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32)
inline uint get_b_fs_zyx_fsv_index(uint b, uint f, uint z, uint y, uint x,
 uint x_size, uint y_size, uint z_size, uint f_size,
 uint b_pad_before, uint b_pad_after,
 uint f_pad_before, uint f_pad_after,
 uint z_pad_before, uint z_pad_after,
 uint y_pad_before, uint y_pad_after,
 uint x_pad_before, uint x_pad_after,
 uint alignment)
{
 const uint feature = f + f_pad_before;
 const uint fs = feature / alignment;
 const uint fsv = feature % alignment;
 const uint x_pitch = alignment;
 const uint y_pitch = x_pitch * (x_pad_before + x_size + x_pad_after);
 const uint z_pitch = y_pitch * (y_pad_before + y_size + y_pad_after);
 const uint fs_pitch = z_pitch * (z_pad_before + z_size + z_pad_after);
 const uint total_f_size = f_pad_before + f_size + f_pad_after;
 const uint b_pitch = fs_pitch * ((total_f_size + alignment - 1) / alignment);
 const uint output_offset = (b_pad_before + b) * b_pitch +
 fs * fs_pitch +
 (z_pad_before + z) * z_pitch +
 (y_pad_before + y) * y_pitch +
 (x_pad_before + x) * x_pitch
 + fsv;
 return output_offset;
}
inline uint get_b_fs_zyx_fsv_index_safe(uint b, uint f, uint z, uint y, uint x,
 uint x_size, uint y_size, uint z_size, uint f_size,
 uint b_pad_before, uint b_pad_after,
 uint f_pad_before, uint f_pad_after,
 uint z_pad_before, uint z_pad_after,
 uint y_pad_before, uint y_pad_after,
 uint x_pad_before, uint x_pad_after,
 uint alignment) {
 const uint f_mod = f_pad_before + (f % f_size);
 const uint fs = f_mod / alignment;
 const uint fsv = f_mod % alignment;
 const uint x_pitch = alignment;
 const uint y_pitch = x_pitch * (x_pad_before + x_size + x_pad_after);
 const uint z_pitch = y_pitch * (y_pad_before + y_size + y_pad_after);
 const uint fs_pitch = z_pitch * (z_pad_before + z_size + z_pad_after);
 const uint total_f_size = f_pad_before + f_size + f_pad_after;
 const uint b_pitch = fs_pitch * ((total_f_size + alignment - 1) / alignment);
 const uint output_offset = (b_pad_before + b) * b_pitch +
 fs * fs_pitch +
 (z_pad_before + (z % z_size)) * z_pitch +
 (y_pad_before + (y % y_size)) * y_pitch +
 (x_pad_before + (x % x_size)) * x_pitch
 + fsv;
 return output_offset;
}
inline uint get_bs_fs_zyx_bsv_fsv_index_safe(uint b, uint f, uint z, uint y, uint x,
 uint x_size, uint y_size, uint z_size, uint f_size, uint b_size,
 uint b_pad_before, uint b_pad_after,
 uint f_pad_before, uint f_pad_after,
 uint z_pad_before, uint z_pad_after,
 uint y_pad_before, uint y_pad_after,
 uint x_pad_before, uint x_pad_after, uint alignmentB, uint alignmentF) {
 const uint b_mod = b_pad_before + (b % b_size);
 const uint f_mod = f_pad_before + (f % f_size);
 const uint bs = b_mod / alignmentB;
 const uint bsv = b_mod % alignmentB;
 const uint fs = f_mod / alignmentF;
 const uint fsv = f_mod % alignmentF;
 const uint x_pitch = alignmentF * alignmentB;
 const uint y_pitch = x_pitch * (x_pad_before + x_size + x_pad_after);
 const uint z_pitch = y_pitch * (y_pad_before + y_size + y_pad_after);
 const uint total_f_size = f_pad_before + f_size + f_pad_after;
 const uint fs_pitch = z_pitch * (z_pad_before + z_size + z_pad_after);
 const uint bs_pitch = fs_pitch * ((total_f_size + alignmentF - 1) / alignmentF);
 const uint output_offset = bs * bs_pitch +
 fs * fs_pitch +
 (z_pad_before + (z % z_size)) * z_pitch +
 (y_pad_before + (y % y_size)) * y_pitch +
 (x_pad_before + (x % x_size)) * x_pitch +
 (bsv * alignmentF)
 + fsv;
 return output_offset;
}
inline uint get_bs_fs_zyx_bsv_fsv_index(uint b, uint f, uint z, uint y, uint x,
 uint x_size, uint y_size, uint z_size, uint f_size,
 uint b_pad_before, uint b_pad_after,
 uint f_pad_before, uint f_pad_after,
 uint z_pad_before, uint z_pad_after,
 uint y_pad_before, uint y_pad_after,
 uint x_pad_before, uint x_pad_after,
 uint b_alignment, uint f_alignment) {
 const uint feature = f + f_pad_before;
 const uint fs = feature / f_alignment;
 const uint fsv = feature % f_alignment;
 const uint bs = (b + b_pad_before) / b_alignment;
 const uint bsv = (b + b_pad_before) % b_alignment;
 const uint bsv_pitch = f_alignment;
 const uint x_pitch = bsv_pitch * b_alignment;
 const uint y_pitch = x_pitch * (x_pad_before + x_size + x_pad_after);
 const uint z_pitch = y_pitch * (y_pad_before + y_size + y_pad_after);
 const uint fs_pitch = z_pitch * (z_pad_before + z_size + z_pad_after);
 const uint total_f_size = f_pad_before + f_size + f_pad_after;
 const uint bs_pitch = fs_pitch * ((total_f_size + f_alignment - 1) / f_alignment);
 const uint output_offset = bs * bs_pitch +
 fs * fs_pitch +
 (z_pad_before + z) * z_pitch +
 (y_pad_before + y) * y_pitch +
 (x_pad_before + x) * x_pitch +
 bsv * bsv_pitch
 + fsv;
 return output_offset;
}
#define GET_DATA_BS_FS_YX_BSV16_FSV16_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 16)
#define GET_DATA_BS_FS_YX_BSV16_FSV32_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 32)
#define GET_DATA_BS_FS_ZYX_BSV32_FSV32_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x,  CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32, 32)
#define GET_DATA_BS_FS_YX_BSV32_FSV32_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32, 32)
#define GET_DATA_BS_FS_YX_BSV4_FSV4_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 4, 4)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV4_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 4)
#define GET_DATA_BS_FS_YX_BSV16_FSV4_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 4)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV8_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 8)
#define GET_DATA_BS_FS_YX_BSV16_FSV8_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 8)
#define GET_DATA_BS_FS_ZYX_BSV8_FSV4_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8, 4)
#define GET_DATA_BS_FS_YX_BSV8_FSV4_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8, 4)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV2_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 2)
#define GET_DATA_BS_FS_YX_BSV16_FSV2_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 2)
#define GET_DATA_BS_FS_ZYX_BSV8_FSV2_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y),  CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8, 2)
#define GET_DATA_BS_FS_YX_BSV8_FSV2_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8, 2)
#define GET_DATA_BS_FS_YX_BSV4_FSV2_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 4, 2)
#define GET_DATA_BS_FS_ZYX_BSV32_FSV16_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32, 16)
#define GET_DATA_BS_FS_YX_BSV32_FSV16_INDEX(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32, 16)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV32_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 32)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV16_INDEX(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 16)
#define GET_DATA_BS_FS_YX_BSV16_FSV16_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 16)
#define GET_DATA_BS_FS_ZYX_BSV32_FSV32_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32, 32)
#define GET_DATA_BS_FS_YX_BSV32_FSV32_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32, 32)
#define GET_DATA_BS_FS_YX_BSV4_FSV4_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 4, 4)
#define GET_DATA_BS_FS_YX_BSV16_FSV4_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x,  CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 4)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV4_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 4)
#define GET_DATA_BS_FS_YX_BSV16_FSV8_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 8)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV8_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 8)
#define GET_DATA_BS_FS_YX_BSV8_FSV4_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8, 4)
#define GET_DATA_BS_FS_ZYX_BSV8_FSV4_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8, 4)
#define GET_DATA_BS_FS_YX_BSV16_FSV2_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 2)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV2_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 2)
#define GET_DATA_BS_FS_YX_BSV8_FSV2_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8, 2)
#define GET_DATA_BS_FS_ZYX_BSV8_FSV2_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 8, 2)
#define GET_DATA_BS_FS_YX_BSV4_FSV2_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z),  CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 4, 2)
#define GET_DATA_BS_FS_ZYX_BSV32_FSV16_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32, 16)
#define GET_DATA_BS_FS_YX_BSV32_FSV16_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 32, 16)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV32_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 32)
#define GET_DATA_BS_FS_YX_BSV16_FSV32_INDEX_SAFE(prefix, b, f, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 32)
#define GET_DATA_BS_FS_ZYX_BSV16_FSV16_INDEX_SAFE(prefix, b, f, z, y, x) get_bs_fs_zyx_bsv_fsv_index_safe( b, f, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _FEATURE_NUM), CAT(prefix, _BATCH_NUM), CAT(prefix, _PAD_BEFORE_BATCH_NUM), CAT(prefix, _PAD_AFTER_BATCH_NUM), CAT(prefix, _PAD_BEFORE_FEATURE_NUM), CAT(prefix, _PAD_AFTER_FEATURE_NUM), CAT(prefix, _PAD_BEFORE_SIZE_Z), CAT(prefix, _PAD_AFTER_SIZE_Z), CAT(prefix, _PAD_BEFORE_SIZE_Y), CAT(prefix, _PAD_AFTER_SIZE_Y), CAT(prefix, _PAD_BEFORE_SIZE_X), CAT(prefix, _PAD_AFTER_SIZE_X), 16, 16)

#define GET_FILTER_OS_IS_YX_ISV_OSV_INDEX(prefix, o, i, y, x, osv, isv) get_os_is_zyx_isv_osv_index( o, i, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 1, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), osv, isv )
#define GET_FILTER_IS_OS_YX_OSV_ISV_INDEX(prefix, o, i, y, x, osv, isv) get_os_is_zyx_isv_osv_index( i, o, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 1, CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), isv, osv )
#define GET_FILTER_IS_OS_YX_ISV_OSV_INDEX(prefix, o, i, y, x, osv, isv) get_is_os_zyx_isv_osv_index( o, i, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 1, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), osv, isv )
#define GET_FILTER_OS_IS_ZYX_ISV_OSV_INDEX(prefix, o, i, z, y, x, osv, isv) get_os_is_zyx_isv_osv_index( o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), osv, isv )
#define GET_FILTER_IS_OS_ZYX_ISV16_OSV16_INDEX(prefix, o, i, z, y, x, sub_group_size) CAT(prefix, _OFFSET) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*(sub_group_size)*CAT(prefix, _X_PITCH) + (y)*(sub_group_size)*CAT(prefix, _Y_PITCH) + (z)*(sub_group_size)*CAT(prefix, _Z_PITCH) + ((i) % (sub_group_size)) + ((o) / (sub_group_size))*(sub_group_size)*CAT(prefix, _OFM_PITCH) + ((i) / (sub_group_size))*CAT(prefix, _IFM_PITCH) )
#define GET_FILTER_IS_OS_YX_ISV16_OSV16_INDEX(prefix, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*(sub_group_size)*CAT(prefix, _X_PITCH) + (y)*(sub_group_size)*CAT(prefix, _Y_PITCH) + ((i) % (sub_group_size)) + ((o) / (sub_group_size))*(sub_group_size)*CAT(prefix, _OFM_PITCH) + ((i) / (sub_group_size))*CAT(prefix, _IFM_PITCH) )
#define GET_FILTER_OS_IS_YX_ISV8_OSV16_ISV2_INDEX(prefix, o, i, y, x, sub_group_size) get_os_is_zyx_isv8_osv16_isv2_index( 0, o, i, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _GROUPS_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), CAT(prefix, _OFFSET) )
#define GET_FILTER_OS_IS_ZYX_ISV8_OSV16_ISV2_INDEX(prefix, o, i, z, y, x, sub_group_size) get_os_is_zyx_isv8_osv16_isv2_index( 0, o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _GROUPS_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), CAT(prefix, _OFFSET) )
inline uint get_os_is_zyx_isv_osv_index(uint o, uint i, uint z, uint y, uint x,
 uint x_size, uint y_size, uint z_size, uint i_size, uint o_size, uint osv_size, uint isv_size)
{
 const uint isv = i % isv_size;
 const uint osv = o % osv_size;
 const uint is = i / isv_size;
 const uint os = o / osv_size;
 const uint x_pitch = osv_size * isv_size;
 const uint y_pitch = x_pitch * x_size;
 const uint z_pitch = y_pitch * y_size;
 const uint is_pitch = z_pitch * z_size;
 const uint os_pitch = is_pitch * ((i_size + isv_size - 1) / isv_size);
 const uint output_offset =
 osv +
 isv * osv_size +
 x * x_pitch +
 y * y_pitch +
 z * z_pitch +
 is * is_pitch +
 os * os_pitch;
 return output_offset;
}
inline uint get_is_os_zyx_isv_osv_index(uint o, uint i, uint z, uint y, uint x,
 uint x_size, uint y_size, uint z_size, uint i_size, uint o_size, uint osv_size, uint isv_size)
{
 const uint isv = i % isv_size;
 const uint osv = o % osv_size;
 const uint is = i / isv_size;
 const uint os = o / osv_size;
 const uint x_pitch = osv_size * isv_size;
 const uint y_pitch = x_pitch * x_size;
 const uint z_pitch = y_pitch * y_size;
 const uint os_pitch = z_pitch * z_size;
 const uint is_pitch = os_pitch * ((o_size + osv_size - 1) / osv_size);
 const uint output_offset =
 osv +
 isv * osv_size +
 x * x_pitch +
 y * y_pitch +
 z * z_pitch +
 os * os_pitch +
 is * is_pitch;
 return output_offset;
}
inline uint get_os_is_zyx_osv_isv_index(uint o, uint i, uint z, uint y, uint x,
 uint x_size, uint y_size, uint z_size, uint i_size, uint o_size, uint osv_size, uint isv_size)
{
 const uint isv = i % isv_size;
 const uint osv = o % osv_size;
 const uint is = i / isv_size;
 const uint os = o / osv_size;
 const uint x_pitch = osv_size * isv_size;
 const uint y_pitch = x_pitch * x_size;
 const uint z_pitch = y_pitch * y_size;
 const uint is_pitch = z_pitch * z_size;
 const uint os_pitch = is_pitch * ((i_size + isv_size - 1) / isv_size);
 const uint output_offset =
 isv +
 osv * isv_size +
 x * x_pitch +
 y * y_pitch +
 z * z_pitch +
 is * is_pitch +
 os * os_pitch;
 return output_offset;
}
inline uint get_g_os_is_zyx_osv_isv_index(uint g, uint o, uint i, uint z, uint y, uint x,
 uint x_size, uint y_size, uint z_size, uint i_size, uint o_size, uint osv_size, uint isv_size)
{
 const uint isv = i % isv_size;
 const uint osv = o % osv_size;
 const uint is = i / isv_size;
 const uint os = o / osv_size;
 const uint x_pitch = osv_size * isv_size;
 const uint y_pitch = x_pitch * x_size;
 const uint z_pitch = y_pitch * y_size;
 const uint is_pitch = z_pitch * z_size;
 const uint os_pitch = is_pitch * ((i_size + isv_size - 1) / isv_size);
 const uint g_pitch = os_pitch * ((o_size + osv_size - 1) / osv_size);
 const uint output_offset =
 isv +
 osv * isv_size +
 x * x_pitch +
 y * y_pitch +
 z * z_pitch +
 is * is_pitch +
 os * os_pitch +
 g * g_pitch;
 return output_offset;
}
#define GET_FILTER_G_OS_IS_ZYX_OSV16_ISV16_INDEX(prefix, g, o, i, z, y, x) get_g_os_is_zyx_osv_isv_index( g, o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), 16, 16)
#define GET_FILTER_OS_IS_YX_OSV16_ISV16_INDEX(prefix, o, i, y, x) get_os_is_zyx_osv_isv_index( o, i, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 1, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), 16, 16)
#define GET_FILTER_OS_IS_ZYX_OSV16_ISV16_INDEX(prefix, o, i, z, y, x) get_os_is_zyx_osv_isv_index( o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), 16, 16)
#define GET_FILTER_OS_IS_ZYX_OSV32_ISV16_INDEX(prefix, o, i, z, y, x) get_os_is_zyx_osv_isv_index( o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), 32, 16)
#define GET_FILTER_OS_IS_ZYX_OSV64_ISV16_INDEX(prefix, o, i, z, y, x) get_os_is_zyx_osv_isv_index( o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), 64, 16)
#define GET_FILTER_G_OS_IS_YX_ISV8_OSV16_ISV2_INDEX(prefix, g, o, i, y, x, sub_group_size) get_os_is_zyx_isv8_osv16_isv2_index( g, o, i, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _GROUPS_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), CAT(prefix, _OFFSET) )
#define GET_FILTER_G_OS_IS_ZYX_ISV8_OSV16_ISV2_INDEX(prefix, g, o, i, z, y, x, sub_group_size) get_os_is_zyx_isv8_osv16_isv2_index( g, o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _GROUPS_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), CAT(prefix, _OFFSET) )
inline uint get_os_is_zyx_isv8_osv16_isv2_index(uint g, uint o, uint i, uint z, uint y, uint x, uint x_size, uint y_size, uint z_size,
 uint g_size, uint o_size, uint i_size, uint offset)
{
 const uint group_offset = g * o_size * i_size * z_size * y_size * x_size;
 const uint xyz_offset = (x + y * x_size + z * x_size * y_size)* 8*16*2;
 const uint i2_val = i % 2;
 const uint i2_slice = i / 2;
 const uint i8_v = i2_slice % 8;
 const uint i8_s = i2_slice / 8;
 const uint i2_offset = i2_val;
 const uint o_offset = (o % 16)*2 + (o / 16) * 16 * i_size * x_size * y_size * z_size;
 const uint i8_offset = 8*16*2* x_size*y_size*z_size * i8_s + 16*2*i8_v;
 const size_t idx = offset + group_offset + xyz_offset + i2_offset + i8_offset + o_offset;
 return idx;
}
inline uint get_os_zyxi_osv16_index(uint o, uint i, uint z, uint y, uint x, uint i_size, uint o_size, uint x_size, uint y_size, uint z_size)
{
 const size_t idx = o%16 + (o / 16)*i_size*x_size*y_size*z_size*16 +
 16 *(i+ x*i_size + y*i_size*x_size + z*i_size*x_size*y_size);
 return idx;
}
#define GET_FILTER_OS_ZYXI_OSV16(prefix, o, i, z, y, x) get_os_zyxi_osv16_index( o, i, z, y, x, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z))
#define GET_FILTER_GOIYX(prefix, g, o, i, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + (o)*CAT(prefix, _OFM_PITCH) + (g)*CAT(prefix, _GROUPS_PITCH)
#define GET_FILTER_GIOYX(prefix, g, o, i, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + (o)*CAT(prefix, _OFM_PITCH) + (g)*CAT(prefix, _GROUPS_PITCH)
#define GET_FILTER_GIOYX_SAFE(prefix, g, o, i, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (i % CAT(prefix, _IFM_NUM))*CAT(prefix, _IFM_PITCH) + (o % CAT(prefix, _OFM_NUM))*CAT(prefix, _OFM_PITCH) + (g % CAT(prefix, _GROUPS_NUM))*CAT(prefix, _GROUPS_PITCH)
#define GET_FILTER_GOIYX_SAFE(prefix, g, o, i, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (i % CAT(prefix, _IFM_NUM))*CAT(prefix, _IFM_PITCH) + (o % CAT(prefix, _OFM_NUM))*CAT(prefix, _OFM_PITCH) + (g % CAT(prefix, _GROUPS_NUM))*CAT(prefix, _GROUPS_PITCH)
#define GET_FILTER_INDEX(prefix, g, o, i, y, x) GET_FILTER_GOIYX(prefix, g, o, i, y, x)
#define GET_FILTER_INDEX_SAFE(prefix, g, o, i, y, x) GET_FILTER_GOIYX_SAFE(prefix, g, o, i, y, x)
#define GET_FILTER_GOIZYX(prefix, g, o, i, z, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (z)*CAT(prefix, _Z_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + (o)*CAT(prefix, _OFM_PITCH) + (g)*CAT(prefix, _GROUPS_PITCH)
#define GET_FILTER_GOIZYX_SAFE(prefix, g, o, i, z, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (z % CAT(prefix, _SIZE_Z ))*CAT(prefix, _Z_PITCH) + (i % CAT(prefix, _IFM_NUM))*CAT(prefix, _IFM_PITCH) + (o % CAT(prefix, _OFM_NUM))*CAT(prefix, _OFM_PITCH) + (g % CAT(prefix, _GROUPS_NUM))*CAT(prefix, _GROUPS_PITCH)
#define GET_FILTER_GIOZYX(prefix, g, o, i, z, y, x) CAT(prefix, _OFFSET) + (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (z)*CAT(prefix, _Z_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + (o)*CAT(prefix, _OFM_PITCH) + (g)*CAT(prefix, _GROUPS_PITCH)
#define GET_FILTER_GIOZYX_SAFE(prefix, g, o, i, z, y, x) CAT(prefix, _OFFSET) + (x % CAT(prefix, _SIZE_X ))*CAT(prefix, _X_PITCH) + (y % CAT(prefix, _SIZE_Y ))*CAT(prefix, _Y_PITCH) + (z % CAT(prefix, _SIZE_Z ))*CAT(prefix, _Z_PITCH) + (i % CAT(prefix, _IFM_NUM))*CAT(prefix, _IFM_PITCH) + (o % CAT(prefix, _OFM_NUM))*CAT(prefix, _OFM_PITCH) + (g % CAT(prefix, _GROUPS_NUM))*CAT(prefix, _GROUPS_PITCH)
#define GET_FILTER_INDEX_5D(prefix, g, o, i, z, y, x) GET_FILTER_GOIZYX(prefix, g, o, i, z, y, x)
#define GET_FILTER_INDEX_5D_SAFE(prefix, g, o, i, z, y, x) GET_FILTER_GOIZYX_SAFE(prefix, g, o, i, z, y, x)
#define GET_FILTER_OS_IYX_OSV_INDEX(prefix, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) +  ((o) % (sub_group_size)) + (sub_group_size)*( (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + ((o) / (sub_group_size))*CAT(prefix, _OFM_PITCH) )
#define GET_FILTER_OS_IYX_OSV_INDEX_INT4_PACKED(prefix, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + ((o) / (sub_group_size))*(CAT(prefix, _OFM_PITCH)/2) )
#define GET_FILTER_OS_IS_YX_OSV_ISV_INDEX_INT4_PACKED(prefix, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + ((o) / (sub_group_size))*(CAT(prefix, _OFM_PITCH)/2) )
#define GET_FILTER_OS_IYX_OSV_ROTATE_180_INDEX(prefix, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + ((o) % (sub_group_size)) + (sub_group_size)*( (CAT(prefix, _SIZE_X ) - x - 1)*CAT(prefix, _X_PITCH) + (CAT(prefix, _SIZE_Y ) - y - 1)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + ((o) / (sub_group_size))*CAT(prefix, _OFM_PITCH) )
inline uint get_gi_yxs_os_yxsv2_osv_index(uint g, uint o, uint i, uint y, uint x, uint x_size, uint g_pitch, uint i_pitch,
 uint y_pitch, uint x_pitch, uint offset, uint sub_group_size)
{
 const uint aligned_ofm_line = x_pitch;
 const uint ifm_height_pitch = (i_pitch/aligned_ofm_line);
 const uint dst_height = i*ifm_height_pitch + y*x_size + x;
 const uint base_filter_index = y*x_size + x;
 const uint aligned_height = dst_height & 0xfffffffe;
 const uint base_filter_odd = (base_filter_index & 0x1);
 uint slice_id = o / sub_group_size;
 uint id_in_slice = o % sub_group_size;
 uint slice_pitch = 2*sub_group_size;
 uint offset_in_slice = (int)(sub_group_size*base_filter_odd);
 const uint in_line = (slice_pitch*slice_id + offset_in_slice + id_in_slice);
 size_t idx = offset + aligned_height*aligned_ofm_line + in_line;
 idx += g * g_pitch;
 return idx;
}
#define GET_FILTER_I_YXS_OS_YXSV2_OSV_INDEX(prefix, o, i, y, x, sub_group_size) get_gi_yxs_os_yxsv2_osv_index( 0, o, i, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _GROUPS_PITCH), CAT(prefix, _IFM_PITCH), CAT(prefix, _Y_PITCH), CAT(prefix, _X_PITCH), CAT(prefix, _OFFSET), sub_group_size)
inline uint get_giy_xs_os_xsv2_osv_index(uint g, uint o, uint i, uint y, uint x, uint x_size, uint g_pitch,
 uint i_pitch, uint y_pitch, uint x_pitch, uint offset, uint sub_group_size)
{
 const uint aligned_ofm_line = x_pitch;
 const uint ifm_height_pitch = (i_pitch/aligned_ofm_line);
 const uint aligned_x_line = y_pitch / x_pitch;
 const uint dst_height = i*ifm_height_pitch + y*aligned_x_line + x;
 const uint base_filter_index = x;
 const uint aligned_height = dst_height & 0xfffffffe;
 const uint base_filter_odd = (base_filter_index & 0x1);
 uint slice_id = o / sub_group_size;
 uint id_in_slice = o % sub_group_size;
 uint slice_pitch = 2*sub_group_size;
 uint offset_in_slice = (int)(sub_group_size*base_filter_odd);
 const bool last_line_in_base_filter = (x == (x_size - 1));
 if (last_line_in_base_filter && base_filter_odd == 0)
 {
 const uint element_in_slice = 32;
 slice_id = o / element_in_slice;
 id_in_slice = o % element_in_slice;
 slice_pitch = 2*element_in_slice;
 offset_in_slice = 0;
 }
 const uint in_line = (slice_pitch*slice_id + offset_in_slice + id_in_slice);
 size_t idx = offset + aligned_height*aligned_ofm_line + in_line;
 idx += g * g_pitch;
 return idx;
}
#define GET_FILTER_IY_XS_OS_XSV2_OSV_INDEX(prefix, o, i, y, x, sub_group_size) get_giy_xs_os_xsv2_osv_index( 0, o, i, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _GROUPS_PITCH), CAT(prefix, _IFM_PITCH), CAT(prefix, _Y_PITCH), CAT(prefix, _X_PITCH), CAT(prefix, _OFFSET), sub_group_size)
inline uint get_is_os_zyx_isa8_osv8_isv2_index(uint o, uint i, uint z, uint y, uint x, uint size_x,
 uint size_y, uint size_z, uint size_ifm, uint size_ofm, uint offset)
{
 const uint isv2_idx = i % 2;
 const uint osv_idx = o % 8;
 const uint isv1_idx = (i / 2) % 8;
 const uint is_idx = i / 16;
 const uint os_idx = o / 8;
 const uint of_8_aligned = ((size_ofm + 7) / 8);
 size_t idx = offset +
 isv2_idx +
 osv_idx * 2 +
 isv1_idx * 8 * 2 +
 x * 8 * 8 * 2 +
 y * size_x * 8 * 8 * 2 +
 z * size_y * size_x * 8 * 8 * 2 +
 os_idx * size_z * size_y * size_x * 8 * 8 * 2 +
 is_idx * of_8_aligned * size_z * size_y * size_x * 8 * 8 * 2;
 return idx;
}
inline uint get_g_os_is_zyx_isa_osv_isv_index(uint g, uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z, uint size_ifm, uint size_ofm, uint offset,
 uint isa, uint osv, uint isv)
{
 const uint isv2_idx = i % isv;
 const uint osv_idx = o % osv;
 const uint isv1_idx = (i / isv) % isa;
 const uint is_idx = i / (isa * isv);
 const uint os_idx = o / osv;
 const uint if_aligned = ((size_ifm + (isa * isv) - 1) / (isa * isv));
 const uint of_aligned = ((size_ofm + (osv - 1)) / osv);
 size_t idx = offset +
 isv2_idx +
 osv_idx * isv +
 isv1_idx * osv * isv +
 x * isa * osv * isv +
 y * size_x * isa * osv * isv +
 z * size_y * size_x * isa * osv * isv +
 is_idx * size_z * size_y * size_x * isa * osv * isv +
 os_idx * if_aligned * size_z * size_y * size_x * isa * osv * isv +
 g * of_aligned * if_aligned * size_z * size_y * size_x * isa * osv * isv;
 return idx;
}
#define GET_FILTER_G_OS_IS_ZYX_ISA_OSV_ISV_INDEX(prefix, g, o, i, z, y, x, isa, osv, isv) get_g_os_is_zyx_isa_osv_isv_index( g, o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET), isa, osv, isv)
inline uint get_g_os_is_yx_isa8_osv8_isv4_index(uint g, uint o, uint i, uint y, uint x, uint size_x,
 uint size_y, uint size_ifm, uint size_ofm, uint offset)
{
 const uint isv2_idx = i % 4;
 const uint osv_idx = o % 8;
 const uint isv1_idx = (i / 4) % 8;
 const uint is_idx = i / 32;
 const uint os_idx = o / 8;
 const uint if_32_aligned = ((size_ifm + 31) / 32);
 const uint of_8_aligned = ((size_ofm + 7) / 8);
 size_t idx = offset +
 isv2_idx +
 osv_idx * 4 +
 isv1_idx * 8 * 4 +
 x * 8 * 8 * 4 +
 y * size_x * 8 * 8 * 4 +
 is_idx * size_y * size_x * 4 * 8 * 8 +
 os_idx * if_32_aligned * size_y * size_x * 4 * 8 * 8 +
 g * of_8_aligned * if_32_aligned * size_y * size_x * 4 * 8 * 8;
 return idx;
}
#define GET_FILTER_OS_IS_YX_ISA8_OSV8_ISV4_INDEX(prefix, o, i, y, x) get_g_os_is_yx_isa8_osv8_isv4_index( 0, o, i, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET))
inline uint get_is_os_yx_isa8_osv8_isv2_index(uint o, uint i, uint y, uint x, uint size_x,
 uint size_y, uint size_ifm, uint size_ofm, uint offset)
{
	const uint isv2_idx = i % 2;
	const uint osv_idx = o % 8;
	const uint isv1_idx = (i / 2) % 8;
	const uint is_idx = i / 16;
	const uint os_idx = o / 8;
 const uint of_8_aligned = ((size_ofm + 7) / 8);
	size_t idx = offset +
 isv2_idx +
 osv_idx * 2 +
 isv1_idx * 8 * 2 +
 x * 8 * 8 * 2 +
 y * size_x * 8 * 8 * 2 +
 os_idx * size_y * size_x * 2 * 8 * 8 +
 is_idx * of_8_aligned * size_y * size_x * 2 * 8 * 8;
 return idx;
}
inline uint get_is_os_yx_isa8_osv8_isv4_index(uint o, uint i, uint y, uint x, uint size_x,
 uint size_y, uint size_ifm, uint size_ofm, uint offset)
{
	const uint isv2_idx = i % 4;
	const uint osv_idx = o % 8;
	const uint isv1_idx = (i / 4) % 8;
	const uint is_idx = i / 32;
	const uint os_idx = o / 8;
 const uint of_8_aligned = ((size_ofm + 7) / 8);
	size_t idx = offset +
 isv2_idx +
 osv_idx * 4 +
 isv1_idx * 8 * 4 +
 x * 8 * 8 * 4 +
 y * size_x * 8 * 8 * 4 +
 os_idx * size_y * size_x * 4 * 8 * 8 +
 is_idx * of_8_aligned * size_y * size_x * 4 * 8 * 8;
 return idx;
}
inline uint get_is_os_yx_osa8_isv16_osv4_index(uint o, uint i, uint y, uint x, uint size_x,
 uint size_y, uint size_ifm, uint size_ofm, uint offset)
{
	const uint osv2_idx = o % 4;
	const uint isv_idx = i % 16;
	const uint osv1_idx = (o / 4) % 8;
	const uint os_idx = o / 32;
	const uint is_idx = i / 16;
 const uint of_32_aligned = ((size_ofm + 31) / 32);
	size_t idx = offset +
 osv2_idx +
 isv_idx * 4 +
 osv1_idx * 16 * 4 +
 x * 8 * 16 * 4 +
 y * size_x * 8 * 16 * 4 +
 os_idx * size_y * size_x * 4 * 16 * 8 +
 is_idx * of_32_aligned * size_y * size_x * 4 * 16 * 8;
 return idx;
}
inline uint get_os_is_zyx_isa8_osv8_isv4_index(uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z,
 uint size_ifm, uint size_ofm, uint offset)
{
 const uint ifm_slices = (size_ifm + 31)/32;
 const uint isv2_idx = i % 4;
 const uint osv_idx = o % 8;
 const uint isv1_idx = (i / 4) % 8;
 const uint is_idx = i / 32;
 const uint os_idx = o / 8;
 size_t idx = offset + isv2_idx + 4 * (osv_idx + 8 * isv1_idx);
 idx += x * 4 * 8 * 8;
 idx += y * size_x * 4 * 8 * 8;
 idx += z * size_y * size_x * 4 * 8 * 8;
 idx += is_idx * size_z * size_y * size_x * 4 * 8 * 8;
 idx += os_idx * ifm_slices * size_z * size_y * size_x * 4 * 8 * 8;
 return idx;
}
#define GET_FILTER_OS_IS_ZYX_ISA8_OSV8_ISV4_INDEX(prefix, o, i, z, y, x) get_os_is_zyx_isa8_osv8_isv4_index( o, i, z, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET))
inline uint get_os_is_yx_isa8_osv16_isv4_index(uint o, uint i, uint y, uint x, uint size_x, uint size_y, uint size_ifm, uint size_ofm, uint offset)
{
 const uint f_32_aligned = ((size_ifm + 31)/32) * 32;
 const uint isv2_idx = i % 4;
 const uint osv_idx = o % 16;
 const uint isv1_idx = (i / 4) % 8;
 const uint is_idx = i / 32;
 const uint os_idx = o / 16;
 size_t idx = offset + isv2_idx + 4 * (osv_idx + 16 * isv1_idx);
 idx += x * 4 * 8 * 16;
 idx += y * size_x * 4 * 8 * 16;
 idx += is_idx * size_y * size_x * 4 * 8 * 16;
 idx += os_idx * (f_32_aligned/32) * size_y * size_x * 4 * 8 * 16;
 return idx;
}
#define GET_FILTER_OS_IS_YX_ISA8_OSV16_ISV4_INDEX(prefix, o, i, y, x) get_os_is_yx_isa8_osv16_isv4_index( o, i, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET))
inline uint get_os_is_zyx_isa8_osv16_isv4_index(uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z,
 uint size_ifm, uint size_ofm, uint offset)
{
 const uint ifm_slices = (size_ifm + 31)/32;
 const uint isv2_idx = i % 4;
 const uint osv_idx = o % 16;
 const uint isv1_idx = (i / 4) % 8;
 const uint is_idx = i / 32;
 const uint os_idx = o / 16;
 size_t idx = offset + isv2_idx + 4 * (osv_idx + 16 * isv1_idx);
 idx += x * 4 * 8 * 16;
 idx += y * size_x * 4 * 8 * 16;
 idx += z * size_y * size_x * 4 * 8 * 16;
 idx += is_idx * size_z * size_y * size_x * 4 * 8 * 16;
 idx += os_idx * ifm_slices * size_z * size_y * size_x * 4 * 8 * 16;
 return idx;
}
#define GET_FILTER_OS_IS_ZYX_ISA8_OSV16_ISV4_INDEX(prefix, o, i, z, y, x) get_os_is_zyx_isa8_osv16_isv4_index( o, i, z, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET))
inline uint get_os_is_yx_osa4_isa8_osv8_isv4_swizzled_by_4_index(uint o, uint i, uint y, uint x, uint size_x, uint size_y, uint size_ifm, uint size_ofm, uint offset)
{
 const uint o_swizzled = (o % 4) * 8 + ((o % 32) / 4) + (o / 32) * 32;
 const uint isv_idx = i % 4;
 const uint isa_idx = (i / 4) % 8;
 const uint is_idx = (i / 32);
 const uint osv_idx = o_swizzled % 8;
 const uint osa_idx = (o_swizzled / 8) % 4;
 const uint os_idx = (o / 32);
 const uint f_32_aligned = ((size_ifm + 31)/32);
 size_t idx = offset +
 isv_idx +
 osv_idx * 4 +
 isa_idx * 8 * 4 +
 osa_idx * 8 * 32 +
 x * 32 * 32 +
 y * size_x * 32 * 32 +
 is_idx * 32 * 32 * size_x * size_y +
 os_idx * 32 * 32 * f_32_aligned * size_x * size_y;
 return idx;
}
inline uint get_os_is_zyx_osa4_isa8_osv8_isv4_swizzled_by_4_index(uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z,
 uint size_ifm, uint size_ofm, uint offset)
{
 const uint o_swizzled = (o % 4) * 8 + ((o % 32) / 4) + (o / 32) * 32;
 const uint isv_idx = i % 4;
 const uint isa_idx = (i / 4) % 8;
 const uint is_idx = (i / 32);
 const uint osv_idx = o_swizzled % 8;
 const uint osa_idx = (o_swizzled / 8) % 4;
 const uint os_idx = (o / 32);
 const uint f_32_aligned = ((size_ifm + 31)/32);
 size_t idx = offset +
 isv_idx +
 osv_idx * 4 +
 isa_idx * 8 * 4 +
 osa_idx * 8 * 32 +
 x * 32 * 32 +
 y * size_x * 32 * 32 +
 z * size_x * size_y * 32 * 32 +
 is_idx * 32 * 32 * size_x * size_y * size_z +
 os_idx * 32 * 32 * f_32_aligned * size_x * size_y * size_z;
 return idx;
}
inline uint get_g_is_os_yx_osa4_isa8_osv8_isv4(uint g, uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z, uint size_ifm, uint size_ofm, uint offset)
{
 const uint isv_idx = i % 4;
 const uint isa_idx = (i / 4) % 8;
 const uint is_idx = (i / 32);
 const uint osv_idx = o % 8;
 const uint osa_idx = (o / 8) % 4;
 const uint os_idx = (o / 32);
 const uint ifm_32_aligned = ((size_ifm + 31) / 32);
 const uint ofm_32_aligned = ((size_ofm + 31) / 32);
 size_t idx = offset +
 isv_idx +
 osv_idx * 4 +
 isa_idx * 8 * 4 +
 osa_idx * 8 * 32 +
 x * 32 * 32 +
 y * size_x * 32 * 32 +
 z * size_y * size_x * 32 * 32 +
 os_idx * 32 * 32 * size_x * size_y * size_z +
 is_idx * 32 * 32 * ofm_32_aligned * size_x * size_y * size_z +
 g * 32 * 32 * ifm_32_aligned * ofm_32_aligned * size_x * size_y * size_z;
 return idx;
}
inline uint get_g_os_is_yx_osa4_isa8_osv8_isv4(uint g, uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z, uint size_ifm, uint size_ofm, uint offset)
{
 const uint isv_idx = i % 4;
 const uint isa_idx = (i / 4) % 8;
 const uint is_idx = (i / 32);
 const uint osv_idx = o % 8;
 const uint osa_idx = (o / 8) % 4;
 const uint os_idx = (o / 32);
 const uint ifm_32_aligned = ((size_ifm + 31)/32);
 const uint ofm_32_aligned = ((size_ofm + 31)/32);
 size_t idx = offset +
 isv_idx +
 osv_idx * 4 +
 isa_idx * 8 * 4 +
 osa_idx * 8 * 32 +
 x * 32 * 32 +
 y * size_x * 32 * 32 +
 z * size_y * size_x * 32 * 32 +
 is_idx * 32 * 32 * size_x * size_y * size_z +
 os_idx * 32 * 32 * ifm_32_aligned * size_x * size_y * size_z +
 g * 32 * 32 * ifm_32_aligned * ofm_32_aligned * size_x * size_y * size_z;
 return idx;
}
inline uint get_g_os_is_yx_osa4_isa8_osv8_isv2(uint g, uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z, uint size_ifm, uint size_ofm, uint offset)
{
 const uint isv_idx = i % 2;
 const uint isa_idx = (i / 2) % 8;
 const uint is_idx = (i / 16);
 const uint osv_idx = o % 8;
 const uint osa_idx = (o / 8) % 4;
 const uint os_idx = (o / 32);
 const uint ifm_16_aligned = ((size_ifm + 15)/16);
 const uint ofm_32_aligned = ((size_ofm + 31)/32);
 size_t idx = offset +
 isv_idx +
 osv_idx * 2 +
 isa_idx * 8 * 2 +
 osa_idx * 8 * 16 +
 x * 32 * 16 +
 y * size_x * 32 * 16 +
 z * size_y * size_x * 32 * 16 +
 is_idx * 32 * 16 * size_x * size_y * size_z +
 os_idx * 32 * 16 * ifm_16_aligned * size_x * size_y * size_z +
 g * 32 * 16 * ifm_16_aligned * ofm_32_aligned * size_x * size_y * size_z;
 return idx;
}
inline uint get_g_os_is_yx_osa2_isa8_osv8_isv2(uint g, uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z, uint size_ifm, uint size_ofm, uint offset)
{
 const uint isv_idx = i % 2;
 const uint isa_idx = (i / 2) % 8;
 const uint is_idx = (i / 16);
 const uint osv_idx = o % 8;
 const uint osa_idx = (o / 8) % 2;
 const uint os_idx = (o / 16);
 const uint ifm_16_aligned = ((size_ifm + 15)/16);
 const uint ofm_16_aligned = ((size_ofm + 15)/16);
 size_t idx = offset +
 isv_idx +
 osv_idx * 2 +
 isa_idx * 8 * 2 +
 osa_idx * 8 * 16 +
 x * 16 * 16 +
 y * size_x * 16 * 16 +
 z * size_y * size_x * 16 * 16 +
 is_idx * 16 * 16 * size_x * size_y * size_z +
 os_idx * 16 * 16 * ifm_16_aligned * size_x * size_y * size_z +
 g * 16 * 16 * ifm_16_aligned * ofm_16_aligned * size_x * size_y * size_z;
 return idx;
}
inline uint get_g_is_os_yx_isa2_osa8_isv8_osv2(uint g, uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z, uint size_ifm, uint size_ofm, uint offset)
{
 return get_g_os_is_yx_osa2_isa8_osv8_isv2(g, i, o, z, y, x, size_x, size_y, size_z, size_ofm, size_ifm, offset);
}
inline uint get_g_is_os_yx_isa4_osa8_isv8_osv4(uint g, uint o, uint i, uint z, uint y, uint x,
 uint size_x, uint size_y, uint size_z, uint size_ifm, uint size_ofm, uint offset)
{
 return get_g_os_is_yx_osa4_isa8_osv8_isv4(g, i, o, z, y, x, size_x, size_y, size_z, size_ofm, size_ifm, offset);
}
#define GET_FILTER_OS_IS_YX_OSA4_ISA8_OSV8_ISV4_INDEX(prefix, o, i, y, x) get_g_os_is_yx_osa4_isa8_osv8_isv4( 0, o, i, 0, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 1, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET))
#define GET_FILTER_OS_IS_ZYX_OSA4_ISA8_OSV8_ISV4_INDEX(prefix, o, i, z, y, x) get_g_os_is_yx_osa4_isa8_osv8_isv4( 0, o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET))
#define GET_FILTER_OS_IS_YX_OSA4_ISA8_OSV8_ISV4_SWIZZLED_BY_4_INDEX(prefix, o, i, y, x) get_os_is_yx_osa4_isa8_osv8_isv4_swizzled_by_4_index( o, i, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET))
#define GET_FILTER_OS_IS_ZYX_OSA4_ISA8_OSV8_ISV4_SWIZZLED_BY_4_INDEX(prefix, o, i, z, y, x) get_os_is_zyx_osa4_isa8_osv8_isv4_swizzled_by_4_index( o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _OFFSET))
inline uint get_is_o32_yx_isv32_swizzled_by_4_index(uint o, uint i, uint y, uint x, uint i_size, uint o_size, uint x_size, uint y_size)
{
 const uint o_aligned_to_32 = ((o_size + 31) / 32) * 32;
 const uint o_swizzled = (o % 4) * 8 + ((o % 32) / 4) + (o / 32) * 32;
 const uint i_aligned_to_32 = ((i_size + 31) / 32) * 32;
 const uint i_val = i % 32;
 const uint i_slice = i / 32;
 const size_t idx = i_val + 32* (x + x_size * (y + y_size * (o_swizzled + o_aligned_to_32 * i_slice) ) );
 return idx;
}
#define GET_FILTER_G_OS_IS_YX_OSV16_ISV4_INDEX(prefix, g, o, i, y, x)  get_g_os_is_yx_osv_isv( g, o, i, y, x, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 16, 4)
inline uint get_g_os_is_yx_osv_isv(uint g, uint o, uint i, uint y, uint x,
 uint i_size,
 uint o_size,
 uint x_size,
 uint y_size,
 uint osv_size,
 uint isv_size)
{
 return get_g_os_is_zyx_osv_isv_index(g, o, i, 0, y, x,
 x_size, y_size, 1, i_size, o_size, osv_size, isv_size);
}
#define GET_FILTER_OS_IS_YX_OSV8_ISV4_INDEX(prefix, o, i, y, x) get_g_os_is_yx_osv_isv( 0, o, i, y, x, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 8, 4)
#define GET_FILTER_OS_IS_YX_OSV16_ISV4_INDEX(prefix, o, i, y, x) get_g_os_is_yx_osv_isv( 0, o, i, y, x, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 16, 4)
#define GET_FILTER_OS_IS_YX_OSV32_ISV4_INDEX(prefix, o, i, y, x) get_g_os_is_yx_osv_isv( 0, o, i, y, x, CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 32, 4)
#define GET_FILTER_OS_IS_ZYX_OSV32_ISV4_INDEX(prefix, o, i, z, y, x) get_os_is_zyx_osv_isv_index( o, i, z, y, x, CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_Z), CAT(prefix, _IFM_NUM), CAT(prefix, _OFM_NUM), 32, 4)
#define GET_FILTER_OS_IS_YX_OSV32_ISV4_SWIZZLED_BY_2_INDEX(prefix, o, i, y, x) get_os_is_yx_osv32_isv4_swizzled_by_2( o, i, y, x, CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), CAT(prefix, _SIZE_Y), CAT(prefix, _SIZE_X))
inline uint get_os_is_yx_osv32_isv4_swizzled_by_2(uint o, uint i, uint y, uint x,
 uint o_size,
 uint i_size,
 uint y_size,
 uint x_size)
{
 const uint osv = 32;
 const uint os = o / osv;
 const uint ofm_block = (o % osv) % 2;
 const uint ofm_in_block = (o % osv) / 2;
 const uint tile = 4;
 const uint ifm_aligned = ((i_size + tile - 1) / tile) * tile;
 const uint ifm_tile = i / tile;
 const uint id = i - ifm_tile * tile;
 uint idx = os * ifm_aligned * y_size * x_size * osv
 + ifm_tile * y_size * x_size * osv * tile
 + y * x_size * osv * tile
 + x * osv * tile
 + ofm_block * 16 * tile
 + ofm_in_block * tile
 + id;
 return idx;
}
inline uint get_os_is_osv32_isv32_swizzled_by_4_index(uint o, uint i, uint y, uint x, uint size_x, uint size_y, uint size_ifm, uint size_ofm, uint offset)
{
 const uint size_ifm_a = ((size_ifm + 31)/32) * 32;
 const uint o_hi = o / 32;
 const uint o_lo = o % 32;
 const uint i_hi = i / 32;
 const uint i_lo = i % 32;
 const uint o_lo1 = o_lo % 4;
 const uint o_lo2 = (o_lo / 4) % 8;
 const uint i_lo1 = i_lo % 4;
 const uint i_lo2 = i_lo / 4;
 const uint idx_in_group = o_lo2 * 4 + o_lo1 * (32 * 8) + i_lo2 * 32 + i_lo1;
 const uint group_idx = o_hi * (size_ifm_a / 32) + i_hi;
 return group_idx * (32 * 32) + idx_in_group;
}
inline uint get_os_i_yxs_osv_yxsv4_index(uint o, uint i, uint y, uint x, uint i_size, uint size_x, uint size_y, uint osv) {
 const uint yxsv = 4;
 uint yx = y * size_x + x;
 uint yx_size_aligned = (size_x * size_y + yxsv - 1) / yxsv * yxsv;
 uint os_index = o / osv;
 uint yxs_index = yx / yxsv;
 uint osv_index = o % osv;
 uint yxsv_index = yx % yxsv;
 uint index = 0;
 index += yxsv_index;
 index += osv_index * yxsv;
 index += yxs_index * yxsv * osv;
 index += i * osv * yx_size_aligned;
 index += os_index * osv * yx_size_aligned * i_size;
 return index;
}
#define GET_FILTER_G_OS_IYX_OSV16(prefix, g, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + (g * CAT(prefix, _GROUPS_PITCH)) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + ((o) / (sub_group_size))*CAT(prefix, _OFM_PITCH) )
#define GET_FILTER_OS_IYX_OSV16(prefix, o, i, y, x, sub_group_size) GET_FILTER_G_OS_IYX_OSV16(prefix, 0, o, i, y, x, sub_group_size)
#define GET_FILTER_GS_OIYX_GSV16(prefix, g, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + ((g) % (sub_group_size)) + (sub_group_size)*( (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + (o)*CAT(prefix, _OFM_PITCH) + ((g) / (sub_group_size))*CAT(prefix, _GROUPS_PITCH) )
#define GET_FILTER_GS_OIZYX_GSV16(prefix, g, o, i, z, y, x, sub_group_size) CAT(prefix, _OFFSET) + ((g) % (sub_group_size)) + (sub_group_size)*( (x)*CAT(prefix, _X_PITCH) + (y)*CAT(prefix, _Y_PITCH) + (z)*CAT(prefix, _Z_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + (o)*CAT(prefix, _OFM_PITCH) + ((g) / (sub_group_size))*CAT(prefix, _GROUPS_PITCH) )
#define GET_FILTER_G_OS_IYX_OSV16_ROTATE_180(prefix, g, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + (g * CAT(prefix, _GROUPS_PITCH)) + ((o) % (sub_group_size)) + (sub_group_size)*( (CAT(prefix, _SIZE_X ) - x - 1)*CAT(prefix, _X_PITCH) + (CAT(prefix, _SIZE_Y ) - y - 1)*CAT(prefix, _Y_PITCH) + (i)*CAT(prefix, _IFM_PITCH) + ((o) / (sub_group_size))*CAT(prefix, _OFM_PITCH) )
#define GET_FILTER_G_IS_OS_ZYX_ISV16_OSV16_INDEX(prefix, g, o, i, z, y, x, sub_group_size) CAT(prefix, _OFFSET) + (g)*CAT(prefix, _GROUPS_PITCH) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*(sub_group_size)*CAT(prefix, _X_PITCH) + (y)*(sub_group_size)*CAT(prefix, _Y_PITCH) + (z)*(sub_group_size)*CAT(prefix, _Z_PITCH) + ((i) % (sub_group_size)) + ((o) / (sub_group_size))*(sub_group_size)*CAT(prefix, _OFM_PITCH) + ((i) / (sub_group_size))*CAT(prefix, _IFM_PITCH) )
#define GET_FILTER_G_IS_OS_YX_ISV16_OSV16_INDEX(prefix, g, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + (g)*CAT(prefix, _GROUPS_PITCH) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*(sub_group_size)*CAT(prefix, _X_PITCH) + (y)*(sub_group_size)*CAT(prefix, _Y_PITCH) + ((i) % (sub_group_size)) + ((o) / (sub_group_size))*(sub_group_size)*CAT(prefix, _OFM_PITCH) + ((i) / (sub_group_size))*CAT(prefix, _IFM_PITCH) )
#define GET_FILTER_G_OS_IS_ZYX_ISV16_OSV16_INDEX(prefix, g, o, i, z, y, x, sub_group_size) CAT(prefix, _OFFSET) + (g)*CAT(prefix, _GROUPS_PITCH) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*(sub_group_size)*CAT(prefix, _X_PITCH) + (y)*(sub_group_size)*CAT(prefix, _Y_PITCH) + (z)*(sub_group_size)*CAT(prefix, _Z_PITCH) + ((i) % (sub_group_size)) + ((i) / (sub_group_size))*(sub_group_size)*CAT(prefix, _IFM_PITCH) + ((o) / (sub_group_size))*CAT(prefix, _OFM_PITCH) )
#define GET_FILTER_GI_YXS_OS_YXSV2_OSV_INDEX(prefix, g, o, i, y, x, sub_group_size) get_gi_yxs_os_yxsv2_osv_index( g, o, i, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _GROUPS_PITCH), CAT(prefix, _IFM_PITCH), CAT(prefix, _Y_PITCH), CAT(prefix, _X_PITCH), CAT(prefix, _OFFSET), sub_group_size)
#define GET_FILTER_GIY_XS_OS_XSV2_OSV_INDEX(prefix, g, o, i, y, x, sub_group_size) get_giy_xs_os_xsv2_osv_index( g, o, i, y, x, CAT(prefix, _SIZE_X ), CAT(prefix, _GROUPS_PITCH), CAT(prefix, _IFM_PITCH), CAT(prefix, _Y_PITCH), CAT(prefix, _X_PITCH), CAT(prefix, _OFFSET), sub_group_size)
inline uint get_gs_oi_yxs_gsv_yxsv4_index(uint g, uint o, uint i, uint y, uint x, uint o_size, uint i_size, uint size_x, uint size_y, const uint gsv) {
 const uint yxsv = 4;
 uint yx = y * size_x + x;
 uint yx_size_aligned = (size_x * size_y + yxsv - 1) / yxsv * yxsv;
 uint gs_index = g / gsv;
 uint yxs_index = yx / yxsv;
 uint gsv_index = g % gsv;
 uint yxsv_index = yx % yxsv;
 uint index = 0;
 index += yxsv_index;
 index += gsv_index * yxsv;
 index += yxs_index * yxsv * gsv;
 index += o * i * gsv * yx_size_aligned;
 index += gs_index * gsv * yx_size_aligned * o_size * i_size;
 return index;
}
#define GET_FILTER_GS_OI_YXS_GSV4_YXSV4_INDEX(prefix, g, o, i, y, x) get_gs_oi_yxs_gsv_yxsv4_index( g, o, i, y, x, CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 4)
#define GET_FILTER_GS_OI_YXS_GSV16_YXSV4_INDEX(prefix, g, o, i, y, x) get_gs_oi_yxs_gsv_yxsv4_index( g, o, i, y, x, CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 16)
#define GET_FILTER_GS_OI_YXS_GSV32_YXSV4_INDEX(prefix, g, o, i, y, x) get_gs_oi_yxs_gsv_yxsv4_index( g, o, i, y, x, CAT(prefix, _OFM_NUM), CAT(prefix, _IFM_NUM), CAT(prefix, _SIZE_X), CAT(prefix, _SIZE_Y), 32)
#define GET_FILTER_G_OS_IS_YX_ISV16_OSV16_INDEX(prefix, g, o, i, y, x, sub_group_size) CAT(prefix, _OFFSET) + (g * CAT(prefix, _GROUPS_PITCH)) + ((o) % (sub_group_size)) + (sub_group_size)*( (x)*(sub_group_size)*CAT(prefix, _X_PITCH) + (y)*(sub_group_size)*CAT(prefix, _Y_PITCH) + ((i) % (sub_group_size)) + ((i) / (sub_group_size))*(sub_group_size)*CAT(prefix, _IFM_PITCH) + ((o) / (sub_group_size))*CAT(prefix, _OFM_PITCH) )
inline uint get_g_os_zyx_is_osv_isv_index(uint g, uint o, uint i, uint z, uint y, uint x,
 uint g_size, uint o_size, uint i_size, uint z_size, uint y_size, uint x_size,
 uint osv, uint isv) {
 uint is_size = (i_size + isv - 1) / isv;
 uint os_size = (o_size + osv - 1) / osv;
 uint isv_index = i % isv;
 uint osv_index = o % osv;
 uint is_index = i / isv;
 uint os_index = o / osv;
 uint isv_pitch = 1;
 uint osv_pitch = isv_pitch * isv;
 uint is_pitch = osv_pitch * osv;
 uint x_pitch = is_pitch * is_size;
 uint y_pitch = x_pitch * x_size;
 uint z_pitch = y_pitch * y_size;
 uint os_pitch = z_pitch * z_size;
 uint g_pitch = os_pitch * os_size;
 uint index = 0;
 index += isv_index * isv_pitch;
 index += osv_index * osv_pitch;
 index += is_index * is_pitch;
 index += x * x_pitch;
 index += y * y_pitch;
 index += z * z_pitch;
 index += os_index * os_pitch;
 index += g * g_pitch;
 return index;
}
#define GET_FILTER_G_OS_ZYX_IS_OSV_ISV_INDEX(tensor, g, o, i, z, y, x, osv, isv) get_g_os_zyx_is_osv_isv_index( g, o, i, z, y, x, CAT(tensor, _GROUPS_NUM), CAT(tensor, _OFM_NUM), CAT(tensor, _IFM_NUM), CAT(tensor, _SIZE_Z), CAT(tensor, _SIZE_Y), CAT(tensor, _SIZE_X), osv, isv)
#define GET_FILTER_G_OS_ZYX_IS_OSV16_ISV4_INDEX(tensor, g, o, i, z, y, x) GET_FILTER_G_OS_ZYX_IS_OSV_ISV_INDEX(tensor, g, o, i, z, y, x, 16, 4)
#define GET_FILTER_G_OS_ZYX_IS_OSV16_ISV16_INDEX(tensor, g, o, i, z, y, x) GET_FILTER_G_OS_ZYX_IS_OSV_ISV_INDEX(tensor, g, o, i, z, y, x, 16, 16)
#define GET_FILTER_G_OS_ZYX_IS_OSV16_ISV32_INDEX(tensor, g, o, i, z, y, x) GET_FILTER_G_OS_ZYX_IS_OSV_ISV_INDEX(tensor, g, o, i, z, y, x, 16, 32)
#define GET_FILTER_G_OS_ZYX_IS_OSV32_ISV4_INDEX(tensor, g, o, i, z, y, x) GET_FILTER_G_OS_ZYX_IS_OSV_ISV_INDEX(tensor, g, o, i, z, y, x, 32, 4)
#define GET_FILTER_G_OS_ZYX_IS_OSV32_ISV16_INDEX(tensor, g, o, i, z, y, x) GET_FILTER_G_OS_ZYX_IS_OSV_ISV_INDEX(tensor, g, o, i, z, y, x, 32, 16)
#define GET_FILTER_G_OS_ZYX_IS_OSV32_ISV32_INDEX(tensor, g, o, i, z, y, x) GET_FILTER_G_OS_ZYX_IS_OSV_ISV_INDEX(tensor, g, o, i, z, y, x, 32, 32)
inline uint get_g_os_y_is_x_osv_isv_index(uint g, uint o, uint i, uint y, uint x,
 uint x_size, uint y_size, uint i_size, uint o_size, uint osv_size, uint isv_size)
{
 const uint isv = i % isv_size;
 const uint osv = o % osv_size;
 const uint is = i / isv_size;
 const uint os = o / osv_size;
 const uint x_pitch = osv_size * isv_size;
 const uint is_pitch = x_pitch * x_size;
 const uint y_pitch = is_pitch * ((i_size + isv_size - 1) / isv_size);
 const uint os_pitch = y_pitch * y_size;
 const uint g_pitch = os_pitch * ((o_size + osv_size - 1) / osv_size);
 const uint output_offset =
 isv +
 osv * isv_size +
 x * x_pitch +
 is * is_pitch +
 y * y_pitch +
 os * os_pitch +
 g * g_pitch;
 return output_offset;
}
#define GET_FILTER_G_OS_Y_IS_X_OSV_ISV_INDEX(tensor, g, o, i, y, x, osv, isv) get_g_os_y_is_x_osv_isv_index( g, o, i, y, x, CAT(tensor, _SIZE_X), CAT(tensor, _SIZE_Y), CAT(tensor, _IFM_NUM), CAT(tensor, _OFM_NUM), osv, isv)
inline uint get_g_os_zy_is_x_osv_isv_index(uint g, uint o, uint i, uint z, uint y, uint x,
 uint o_size, uint i_size, uint z_size, uint y_size, uint x_size,
 uint osv, uint isv) {
 uint is_size = (i_size + isv - 1) / isv;
 uint os_size = (o_size + osv - 1) / osv;
 uint isv_index = i % isv;
 uint osv_index = o % osv;
 uint is_index = i / isv;
 uint os_index = o / osv;
 uint isv_pitch = 1;
 uint osv_pitch = isv_pitch * isv;
 uint x_pitch = osv_pitch * osv;
 uint is_pitch = x_pitch * x_size;
 uint y_pitch = is_pitch * is_size;
 uint z_pitch = y_pitch * y_size;
 uint os_pitch = z_pitch * z_size;
 uint g_pitch = os_pitch * os_size;
 uint index = 0;
 index += isv_index * isv_pitch;
 index += osv_index * osv_pitch;
 index += is_index * is_pitch;
 index += x * x_pitch;
 index += y * y_pitch;
 index += z * z_pitch;
 index += os_index * os_pitch;
 index += g * g_pitch;
 return index;
}
#define GET_FILTER_G_OS_ZY_IS_X_OSV_ISV_INDEX(tensor, g, o, i, z, y, x, osv, isv) get_g_os_zy_is_x_osv_isv_index( g, o, i, z, y, x, CAT(tensor, _OFM_NUM), CAT(tensor, _IFM_NUM), CAT(tensor, _SIZE_Z), CAT(tensor, _SIZE_Y), CAT(tensor, _SIZE_X), osv, isv)

inline int imad_SW(int acc, uchar4 input, char4 weight) __attribute__((overloadable)) {
 acc += input[0] * weight[0];
 acc += input[1] * weight[1];
 acc += input[2] * weight[2];
 acc += input[3] * weight[3];
 return acc;
}
inline int imad_SW(int acc, char4 input, char4 weight) __attribute__((overloadable)) {
 acc += input[0] * weight[0];
 acc += input[1] * weight[1];
 acc += input[2] * weight[2];
 acc += input[3] * weight[3];
 return acc;
}
inline int imad_SW(int acc, char4 input, uchar4 weight) __attribute__((overloadable)) {
 acc += input[0] * weight[0];
 acc += input[1] * weight[1];
 acc += input[2] * weight[2];
 acc += input[3] * weight[3];
 return acc;
}
inline int imad_SW(int acc, uchar4 input, uchar4 weight) __attribute__((overloadable)) {
 acc += input[0] * weight[0];
 acc += input[1] * weight[1];
 acc += input[2] * weight[2];
 acc += input[3] * weight[3];
 return acc;
}
#define IMAD(_O, _I, _W) imad_SW(_O, _I, _W)

typedef struct __attribute__ ((packed)) int4x2_t { char s0; } int4x2_t;
typedef struct __attribute__ ((packed)) int4x4_t { int4x2_t s0; int4x2_t s1; } int4x4_t;
typedef struct __attribute__ ((packed)) int4x8_t { int4x2_t s0; int4x2_t s1; int4x2_t s2; int4x2_t s3; } int4x8_t;
typedef struct __attribute__ ((packed)) int4x16_t { int4x2_t s0; int4x2_t s1; int4x2_t s2; int4x2_t s3; int4x2_t s4; int4x2_t s5; int4x2_t s6; int4x2_t s7; } int4x16_t;
typedef struct __attribute__ ((packed)) uint4x2_t { uchar s0; } uint4x2_t;
typedef struct __attribute__ ((packed)) uint4x4_t { uint4x2_t s0; uint4x2_t s1; } uint4x4_t;
typedef struct __attribute__ ((packed)) uint4x8_t { uint4x2_t s0; uint4x2_t s1; uint4x2_t s2; uint4x2_t s3; } uint4x8_t;
typedef struct __attribute__ ((packed)) uint4x16_t { uint4x2_t s0; uint4x2_t s1; uint4x2_t s2; uint4x2_t s3; uint4x2_t s4; uint4x2_t s5; uint4x2_t s6; uint4x2_t s7; } uint4x16_t;
inline uchar2 cvt_uint4x2_to_uint8x2(uint4x2_t v) __attribute__((overloadable)) {
 const uchar v0 = v.s0 & 0x0F;
 const uchar v1 = (v.s0 & 0xF0) >> 4;
 return (uchar2)(v0, v1);
}
inline char2 cvt_uint4x2_to_int8x2(uint4x2_t v) __attribute__((overloadable)) {
 const char v0 = convert_char(v.s0 & 0x0F);
 const char v1 = convert_char((v.s0 & 0xF0) >> 4);
 return (char2)(v0, v1);
}
inline char2 cvt_int4x2_to_int8x2(int4x2_t v) __attribute__((overloadable)) {
 const char s_bit = (v.s0 & convert_char(0x08));
 const char mask = s_bit > 0 ? convert_char(0xF0) : convert_char(0x00);
 const char v0 = (v.s0 & convert_char(0x0F)) | mask;
 const char v1 = v.s0 >> 4;
 return (char2)(v0, v1);
}
inline uchar2 unpack_to_uchar(uint4x2_t v) __attribute__((overloadable)) {
 return cvt_uint4x2_to_uint8x2(v);
}
inline char2 unpack_to_char(int4x2_t v) __attribute__((overloadable)) {
 return cvt_int4x2_to_int8x2(v);
}
inline char2 unpack_to_char(uint4x2_t v) __attribute__((overloadable)) {
 return convert_char2(cvt_uint4x2_to_uint8x2(v));
}
inline char4 unpack_to_char(int4x4_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s1);
 return (char4)(v0.s0, v0.s1, v1.s0, v1.s1);
}
inline char4 unpack_to_char(uint4x4_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s1);
 return (char4)(v0.s0, v0.s1, v1.s0, v1.s1);
}
inline uchar4 unpack_to_uchar(uint4x4_t v) __attribute__((overloadable)) {
 uchar2 v0 = unpack_to_uchar(v.s0);
 uchar2 v1 = unpack_to_uchar(v.s1);
 return (uchar4)(v0.s0, v0.s1, v1.s0, v1.s1);
}
inline char4 unpack_transposed_to_char(int4x4_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s1);
 return (char4)(v0.s0, v1.s0, v0.s1, v1.s1);
}
inline char4 unpack_transposed_to_char(uint4x4_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s1);
 return (char4)(v0.s0, v1.s0, v0.s1, v1.s1);
}
inline uchar4 unpack_transposed_to_uchar(uint4x4_t v) __attribute__((overloadable)) {
 uchar2 v0 = unpack_to_uchar(v.s0);
 uchar2 v1 = unpack_to_uchar(v.s1);
 return (uchar4)(v0.s0, v1.s0, v0.s1, v1.s1);
}
inline uchar8 unpack_to_uchar(uint4x8_t v) __attribute__((overloadable)) {
 uchar2 v0 = unpack_to_uchar(v.s0);
 uchar2 v1 = unpack_to_uchar(v.s1);
 uchar2 v2 = unpack_to_uchar(v.s2);
 uchar2 v3 = unpack_to_uchar(v.s3);
 return (uchar8)(v0.s0, v0.s1, v1.s0, v1.s1, v2.s0, v2.s1, v3.s0, v3.s1);
}
inline char8 unpack_to_char(int4x8_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s1);
 char2 v2 = unpack_to_char(v.s2);
 char2 v3 = unpack_to_char(v.s3);
 return (char8)(v0.s0, v0.s1, v1.s0, v1.s1, v2.s0, v2.s1, v3.s0, v3.s1);
}
inline char8 unpack_to_char(uint4x8_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s1);
 char2 v2 = unpack_to_char(v.s2);
 char2 v3 = unpack_to_char(v.s3);
 return (char8)(v0.s0, v0.s1, v1.s0, v1.s1, v2.s0, v2.s1, v3.s0, v3.s1);
}
inline char8 unpack_transposed_to_char(int4x8_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s1);
 char2 v2 = unpack_to_char(v.s2);
 char2 v3 = unpack_to_char(v.s3);
 return (char8)(v0.s0, v1.s0, v2.s0, v3.s0, v0.s1, v1.s1, v2.s1, v3.s1);
}
inline char8 unpack_transposed_to_char(uint4x8_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s1);
 char2 v2 = unpack_to_char(v.s2);
 char2 v3 = unpack_to_char(v.s3);
 return (char8)(v0.s0, v1.s0, v2.s0, v3.s0, v0.s1, v1.s1, v2.s1, v3.s1);
}
inline uchar8 unpack_transposed_to_uchar(uint4x8_t v) __attribute__((overloadable)) {
 uchar2 v0 = unpack_to_uchar(v.s0);
 uchar2 v1 = unpack_to_uchar(v.s1);
 uchar2 v2 = unpack_to_uchar(v.s2);
 uchar2 v3 = unpack_to_uchar(v.s3);
 return (uchar8)(v0.s0, v1.s0, v2.s0, v3.s0, v0.s1, v1.s1, v2.s1, v3.s1);
}
inline float2 unpack_to_float(uint4x2_t v) __attribute__((overloadable)) {
 return convert_float2(cvt_uint4x2_to_uint8x2(v));
}
inline float2 unpack_to_float(int4x2_t v) __attribute__((overloadable)) {
 return convert_float2(cvt_int4x2_to_int8x2(v));
}
inline float4 unpack_to_float(uint4x4_t v) __attribute__((overloadable)) {
 float2 f0 = unpack_to_float(v.s0);
 float2 f1 = unpack_to_float(v.s1);
 return (float4)(f0.s0, f0.s1, f1.s0, f1.s1);
}
inline float4 unpack_to_float(int4x4_t v) __attribute__((overloadable)) {
 float2 f0 = unpack_to_float(v.s0);
 float2 f1 = unpack_to_float(v.s1);
 return (float4)(f0.s0, f0.s1, f1.s0, f1.s1);
}
inline float8 unpack_to_float(uint4x8_t v) __attribute__((overloadable)) {
 float2 f0 = unpack_to_float(v.s0);
 float2 f1 = unpack_to_float(v.s1);
 float2 f2 = unpack_to_float(v.s2);
 float2 f3 = unpack_to_float(v.s3);
 return (float8)(f0.s0, f0.s1, f1.s0, f1.s1, f2.s0, f2.s1, f3.s0, f3.s1);
}
inline float8 unpack_to_float(int4x8_t v) __attribute__((overloadable)) {
 float2 f0 = unpack_to_float(v.s0);
 float2 f1 = unpack_to_float(v.s1);
 float2 f2 = unpack_to_float(v.s2);
 float2 f3 = unpack_to_float(v.s3);
 return (float8)(f0.s0, f0.s1, f1.s0, f1.s1, f2.s0, f2.s1, f3.s0, f3.s1);
}
#if defined(cl_khr_fp16)
inline half2 unpack_to_half(uint4x2_t v) __attribute__((overloadable)) {
 return convert_half2(cvt_uint4x2_to_uint8x2(v));
}
inline half2 unpack_to_half(int4x2_t v) __attribute__((overloadable)) {
 return convert_half2(cvt_int4x2_to_int8x2(v));
}
inline half4 unpack_to_half(uint4x4_t v) __attribute__((overloadable)) {
 half2 f0 = unpack_to_half(v.s0);
 half2 f1 = unpack_to_half(v.s1);
 return (half4)(f0.s0, f0.s1, f1.s0, f1.s1);
}
inline half4 unpack_to_half_osv32_isv2(uint4x4_t v) __attribute__((overloadable)) {
 half2 f0 = unpack_to_half(v.s0);
 half2 f1 = unpack_to_half(v.s1);
 return (half4)(f0.s0, f0.s1, f1.s0, f1.s1);
}
inline half4 unpack_to_half(int4x4_t v) __attribute__((overloadable)) {
 half2 f0 = unpack_to_half(v.s0);
 half2 f1 = unpack_to_half(v.s1);
 return (half4)(f0.s0, f0.s1, f1.s0, f1.s1);
}

inline half4 unpack_to_half_osv32_isv2(int4x4_t v) __attribute__((overloadable)) {
 half2 f0 = unpack_to_half(v.s0);
 half2 f1 = unpack_to_half(v.s1);
 return (half4)(f0.s0, f0.s1, f1.s0, f1.s1);
}
inline half8 unpack_to_half(uint4x8_t v) __attribute__((overloadable)) {
 half2 f0 = unpack_to_half(v.s0);
 half2 f1 = unpack_to_half(v.s1);
 half2 f2 = unpack_to_half(v.s2);
 half2 f3 = unpack_to_half(v.s3);
 return (half8)(f0.s0, f0.s1, f1.s0, f1.s1, f2.s0, f2.s1, f3.s0, f3.s1);
}
inline half8 unpack_to_half_osv32_isv2(uint4x8_t v) __attribute__((overloadable)) {
 half2 f0 = unpack_to_half(v.s0);
 half2 f1 = unpack_to_half(v.s2);
 half2 f2 = unpack_to_half(v.s1);
 half2 f3 = unpack_to_half(v.s3);
 return (half8)(f0.s0, f0.s1, f1.s0, f1.s1, f2.s0, f2.s1, f3.s0, f3.s1);
}
inline half8 unpack_to_half(int4x8_t v) __attribute__((overloadable)) {
 half2 f0 = unpack_to_half(v.s0);
 half2 f1 = unpack_to_half(v.s1);
 half2 f2 = unpack_to_half(v.s2);
 half2 f3 = unpack_to_half(v.s3);
 return (half8)(f0.s0, f0.s1, f1.s0, f1.s1, f2.s0, f2.s1, f3.s0, f3.s1);
}
inline half8 unpack_to_half_osv32_isv2(int4x8_t v) __attribute__((overloadable)) {
 half2 f0 = unpack_to_half(v.s0);
 half2 f1 = unpack_to_half(v.s2);
 half2 f2 = unpack_to_half(v.s1);
 half2 f3 = unpack_to_half(v.s3);
 return (half8)(f0.s0, f0.s1, f1.s0, f1.s1, f2.s0, f2.s1, f3.s0, f3.s1);
}
inline char8 unpack_to_char_osv32_isv2(int4x8_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s2);
 char2 v2 = unpack_to_char(v.s1);
 char2 v3 = unpack_to_char(v.s3);
 return (char8)(v0.s0, v0.s1, v1.s0, v1.s1, v2.s0, v2.s1, v3.s0, v3.s1);
}
inline char8 unpack_to_char_osv32_isv2(uint4x8_t v) __attribute__((overloadable)) {
 char2 v0 = unpack_to_char(v.s0);
 char2 v1 = unpack_to_char(v.s2);
 char2 v2 = unpack_to_char(v.s1);
 char2 v3 = unpack_to_char(v.s3);
 return (char8)(v0.s0, v0.s1, v1.s0, v1.s1, v2.s0, v2.s1, v3.s0, v3.s1);
}
inline uchar8 unpack_to_uchar_osv32_isv2(uint4x8_t v) __attribute__((overloadable)) {
 uchar2 v0 = unpack_to_uchar(v.s0);
 uchar2 v1 = unpack_to_uchar(v.s2);
 uchar2 v2 = unpack_to_uchar(v.s1);
 uchar2 v3 = unpack_to_uchar(v.s3);
 return (uchar8)(v0.s0, v0.s1, v1.s0, v1.s1, v2.s0, v2.s1, v3.s0, v3.s1);
}
#endif
#define UNPACK_INT4x2(target_type, value) CAT(unpack_to_, target_type)(value)
#define UNPACK_INT4x2_OSV32_ISV2(target_type, value) CAT(CAT(unpack_to_, target_type), _osv32_isv2)(value)
#define UNPACK_INT4x4_OSV32_ISV2(target_type, value) CAT(CAT(unpack_to_, target_type), _osv32_isv2)(value)
#define UNPACK_TRANSPOSED_INT4x2(target_type, value) CAT(unpack_transposed_to_, target_type)(value)

#define BLOCK_READ_TYPE_size1 uchar
#define BLOCK_READ_TYPE_size2 ushort
#define BLOCK_READ_TYPE_size4 uint
#define BLOCK_READ_TYPE_size8 ulong
#define BLOCK_READ_TYPE(type_size) CAT(BLOCK_READ_TYPE_size, type_size)
#define BLOCK_READ_FUNC_size1 _sub_group_block_read_uc
#define BLOCK_READ_FUNC_size2 _sub_group_block_read_us
#define BLOCK_READ_FUNC_size4 _sub_group_block_read
#define BLOCK_READ_FUNC_size8 _sub_group_block_read_ul
#define BLOCK_READ_FUNC(type_size) CAT(BLOCK_READ_FUNC_size, type_size)
#define BLOCK_READN_FUNC_SIZE_DEF(type_size, vector_size) MAKE_VECTOR_TYPE(BLOCK_READ_FUNC(type_size), vector_size)
#define BLOCK_READN_FUNC_size1(vector_size) BLOCK_READN_FUNC_SIZE_DEF(1, vector_size)
#define BLOCK_READN_FUNC_size2(vector_size) BLOCK_READN_FUNC_SIZE_DEF(2, vector_size)
#define BLOCK_READN_FUNC_size4(vector_size) BLOCK_READN_FUNC_SIZE_DEF(4, vector_size)
#define BLOCK_READN_FUNC_size8(vector_size) BLOCK_READN_FUNC_SIZE_DEF(8, vector_size)
#define BLOCK_READN_FUNC(type_size, vector_size) CAT(BLOCK_READN_FUNC_size, type_size)(vector_size)
#define BLOCK_READN_RAW(type_size, vector_size, addr_space, ptr, offset) BLOCK_READN_FUNC(type_size, vector_size)((const addr_space BLOCK_READ_TYPE(type_size)*)(ptr) + (offset))
#define BLOCK_READN(type, vector_size, ptr, offset) AS_TYPE(MAKE_VECTOR_TYPE(type, vector_size), BLOCK_READN_RAW(TYPE_SIZE(type), vector_size, __global, ptr, offset))
#define BLOCK_READN_SLM(type, vector_size, ptr, offset) AS_TYPE(MAKE_VECTOR_TYPE(type, vector_size), BLOCK_READN_RAW(TYPE_SIZE(type), vector_size, __local, ptr, offset))
#define DT_INPUT_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT0_TYPE, 1, ptr, offset)
#define DT_INPUT_BLOCK_READ2(ptr, offset) BLOCK_READN(INPUT0_TYPE, 2, ptr, offset)
#define DT_INPUT_BLOCK_READ4(ptr, offset) BLOCK_READN(INPUT0_TYPE, 4, ptr, offset)
#define DT_INPUT_BLOCK_READ8(ptr, offset) BLOCK_READN(INPUT0_TYPE, 8, ptr, offset)
#define DT_INPUT_BLOCK_READ16(ptr, offset) BLOCK_READN(INPUT0_TYPE, 16, ptr, offset)
#define DT_BIAS_BLOCK_READ(ptr, offset) BLOCK_READN(BIAS_TYPE, 1, ptr, offset)
#define DT_BIAS_BLOCK_READ2(ptr, offset) BLOCK_READN(BIAS_TYPE, 2, ptr, offset)
#define DT_BIAS_BLOCK_READ4(ptr, offset) BLOCK_READN(BIAS_TYPE, 4, ptr, offset)
#define DT_BIAS_BLOCK_READ8(ptr, offset) BLOCK_READN(BIAS_TYPE, 8, ptr, offset)
#define DT_BIAS_BLOCK_READ16(ptr, offset) BLOCK_READN(BIAS_TYPE, 16, ptr, offset)
#define DT_FILTER_BLOCK_READ(ptr, offset) BLOCK_READN(FILTER_TYPE, 1, ptr, offset)
#define DT_FILTER_BLOCK_READ2(ptr, offset) BLOCK_READN(FILTER_TYPE, 2, ptr, offset)
#define DT_FILTER_BLOCK_READ4(ptr, offset) BLOCK_READN(FILTER_TYPE, 4, ptr, offset)
#define DT_FILTER_BLOCK_READ8(ptr, offset) BLOCK_READN(FILTER_TYPE, 8, ptr, offset)
#define DT_FILTER_BLOCK_READ16(ptr, offset) BLOCK_READN(FILTER_TYPE, 16, ptr, offset)
#define BLOCK_READ_IMPL_1 ret = ptr[idx];
#define BLOCK_READ_IMPL_2 ret.s0 = ptr[idx]; idx += get_max_sub_group_size(); ret.s1 = ptr[idx]; idx += get_max_sub_group_size();
#define BLOCK_READ_IMPL_4 BLOCK_READ_IMPL_2 ret.s2 = ptr[idx]; idx += get_max_sub_group_size(); ret.s3 = ptr[idx]; idx += get_max_sub_group_size();
#define BLOCK_READ_IMPL_8 BLOCK_READ_IMPL_4 ret.s4 = ptr[idx]; idx += get_max_sub_group_size(); ret.s5 = ptr[idx]; idx += get_max_sub_group_size(); ret.s6 = ptr[idx]; idx += get_max_sub_group_size(); ret.s7 = ptr[idx]; idx += get_max_sub_group_size();
#define BLOCK_READ_IMPL_16 BLOCK_READ_IMPL_8 ret.s8 = ptr[idx]; idx += get_max_sub_group_size(); ret.s9 = ptr[idx]; idx += get_max_sub_group_size(); ret.sa = ptr[idx]; idx += get_max_sub_group_size(); ret.sb = ptr[idx]; idx += get_max_sub_group_size(); ret.sc = ptr[idx]; idx += get_max_sub_group_size(); ret.sd = ptr[idx]; idx += get_max_sub_group_size(); ret.se = ptr[idx]; idx += get_max_sub_group_size(); ret.sf = ptr[idx]; idx += get_max_sub_group_size();
#define BLOCK_READ_IMPL(vec_size) CAT(BLOCK_READ_IMPL_, vec_size)
#define BLOCK_READ_FUNC_NAME(type_size, vec_size) MAKE_VECTOR_TYPE(BLOCK_READ_FUNC(type_size), vec_size)
#define DECLARE_BLOCK_READ_EMULATION(type_size, vec_size) inline MAKE_VECTOR_TYPE(BLOCK_READ_TYPE(type_size), vec_size) BLOCK_READ_FUNC_NAME(type_size, vec_size)(const __global BLOCK_READ_TYPE(type_size)* ptr) { uint idx = get_sub_group_local_id(); MAKE_VECTOR_TYPE(BLOCK_READ_TYPE(type_size), vec_size) ret; BLOCK_READ_IMPL(vec_size) return ret; }
#if defined(cl_intel_subgroups)
 #define _sub_group_block_read(ptr) intel_sub_group_block_read(ptr)
 #define _sub_group_block_read2(ptr) intel_sub_group_block_read2(ptr)
 #define _sub_group_block_read4(ptr) intel_sub_group_block_read4(ptr)
 #define _sub_group_block_read8(ptr) intel_sub_group_block_read8(ptr)
#elif (__OPENCL_C_VERSION__ >= 200)
 DECLARE_BLOCK_READ_EMULATION(4, 1)
 DECLARE_BLOCK_READ_EMULATION(4, 2)
 DECLARE_BLOCK_READ_EMULATION(4, 4)
 DECLARE_BLOCK_READ_EMULATION(4, 8)
#endif
#if defined(cl_intel_subgroups_short)
 #define _sub_group_block_read_us(ptr) intel_sub_group_block_read_us(ptr)
 #define _sub_group_block_read_us2(ptr) intel_sub_group_block_read_us2(ptr)
 #define _sub_group_block_read_us4(ptr) intel_sub_group_block_read_us4(ptr)
 #define _sub_group_block_read_us8(ptr) intel_sub_group_block_read_us8(ptr)
#elif (__OPENCL_C_VERSION__ >= 200)
 DECLARE_BLOCK_READ_EMULATION(2, 1)
 DECLARE_BLOCK_READ_EMULATION(2, 2)
 DECLARE_BLOCK_READ_EMULATION(2, 4)
 DECLARE_BLOCK_READ_EMULATION(2, 8)
#endif
#if defined(cl_intel_subgroups_char)
 #define _sub_group_block_read_uc(ptr) intel_sub_group_block_read_uc(ptr)
 #define _sub_group_block_read_uc2(ptr) intel_sub_group_block_read_uc2(ptr)
 #define _sub_group_block_read_uc4(ptr) intel_sub_group_block_read_uc4(ptr)
 #define _sub_group_block_read_uc8(ptr) intel_sub_group_block_read_uc8(ptr)
 #define _sub_group_block_read_uc16(ptr) intel_sub_group_block_read_uc16(ptr)
#elif (__OPENCL_C_VERSION__ >= 200)
 DECLARE_BLOCK_READ_EMULATION(1, 1)
 DECLARE_BLOCK_READ_EMULATION(1, 2)
 DECLARE_BLOCK_READ_EMULATION(1, 4)
 DECLARE_BLOCK_READ_EMULATION(1, 8)
 DECLARE_BLOCK_READ_EMULATION(1, 16)
#endif
#if defined(cl_intel_subgroups_long)
 #define _sub_group_block_read_ul(ptr) intel_sub_group_block_read_ul(ptr)
 #define _sub_group_block_read_ul2(ptr) intel_sub_group_block_read_ul2(ptr)
 #define _sub_group_block_read_ul4(ptr) intel_sub_group_block_read_ul4(ptr)
 #define _sub_group_block_read_ul8(ptr) intel_sub_group_block_read_ul8(ptr)
#elif (__OPENCL_C_VERSION__ >= 200)
 DECLARE_BLOCK_READ_EMULATION(8, 1)
 DECLARE_BLOCK_READ_EMULATION(8, 2)
 DECLARE_BLOCK_READ_EMULATION(8, 4)
 DECLARE_BLOCK_READ_EMULATION(8, 8)
#endif

#define BLOCK_WRITE_TYPE_size1 uchar
#define BLOCK_WRITE_TYPE_size2 ushort
#define BLOCK_WRITE_TYPE_size4 uint
#define BLOCK_WRITE_TYPE_size8 ulong
#define BLOCK_WRITE_TYPE(type_size) CAT(BLOCK_WRITE_TYPE_size, type_size)
#define BLOCK_WRITE_FUNC_size1 _sub_group_block_write_uc
#define BLOCK_WRITE_FUNC_size2 _sub_group_block_write_us
#define BLOCK_WRITE_FUNC_size4 _sub_group_block_write
#define BLOCK_WRITE_FUNC_size8 _sub_group_block_write_ul
#define BLOCK_WRITE_FUNC(type_size) CAT(BLOCK_WRITE_FUNC_size, type_size)
#define BLOCK_WRITEN_FUNC_SIZE_DEF(type_size, vector_size) MAKE_VECTOR_TYPE(BLOCK_WRITE_FUNC(type_size), vector_size)
#define BLOCK_WRITEN_FUNC_size1(vector_size) BLOCK_WRITEN_FUNC_SIZE_DEF(1, vector_size)
#define BLOCK_WRITEN_FUNC_size2(vector_size) BLOCK_WRITEN_FUNC_SIZE_DEF(2, vector_size)
#define BLOCK_WRITEN_FUNC_size4(vector_size) BLOCK_WRITEN_FUNC_SIZE_DEF(4, vector_size)
#define BLOCK_WRITEN_FUNC_size8(vector_size) BLOCK_WRITEN_FUNC_SIZE_DEF(8, vector_size)
#define BLOCK_WRITEN_FUNC(type_size, vector_size) CAT(BLOCK_WRITEN_FUNC_size, type_size)(vector_size)
#define BLOCK_WRITEN_RAW(type_size, vector_size, addr_space, ptr, offset, val) BLOCK_WRITEN_FUNC(type_size, vector_size)( (addr_space BLOCK_WRITE_TYPE(type_size)*)(ptr) + (offset), AS_TYPE(MAKE_VECTOR_TYPE(BLOCK_WRITE_TYPE(type_size), vector_size), val))
#define BLOCK_WRITEN(type, vector_size, ptr, offset, val) BLOCK_WRITEN_RAW(TYPE_SIZE(type), vector_size, __global, ptr, offset, val)
#define BLOCK_WRITEN_SLM(type, vector_size, ptr, offset, val) BLOCK_WRITEN_RAW(TYPE_SIZE(type), vector_size, __local, ptr, offset, val)
#define DT_OUTPUT_BLOCK_WRITE(ptr, offset, val) BLOCK_WRITEN(OUTPUT_TYPE, 1, ptr, offset, val)
#define DT_OUTPUT_BLOCK_WRITE2(ptr, offset, val) BLOCK_WRITEN(OUTPUT_TYPE, 2, ptr, offset, val)
#define DT_OUTPUT_BLOCK_WRITE4(ptr, offset, val) BLOCK_WRITEN(OUTPUT_TYPE, 4, ptr, offset, val)
#define DT_OUTPUT_BLOCK_WRITE8(ptr, offset, val) BLOCK_WRITEN(OUTPUT_TYPE, 8, ptr, offset, val)
#define DT_OUTPUT_BLOCK_WRITE16(ptr, offset, val) BLOCK_WRITEN(OUTPUT_TYPE, 16, ptr, offset, val)
#define BLOCK_WRITE_IMPL_1 out_ptr[idx] = v;
#define BLOCK_WRITE_IMPL_2 out_ptr[idx] = v.s0; idx += get_max_sub_group_size(); out_ptr[idx] = v.s1; idx += get_max_sub_group_size();
#define BLOCK_WRITE_IMPL_4 BLOCK_WRITE_IMPL_2 out_ptr[idx] = v.s2; idx += get_max_sub_group_size(); out_ptr[idx] = v.s3; idx += get_max_sub_group_size();
#define BLOCK_WRITE_IMPL_8 BLOCK_WRITE_IMPL_4 out_ptr[idx] = v.s4; idx += get_max_sub_group_size(); out_ptr[idx] = v.s5; idx += get_max_sub_group_size(); out_ptr[idx] = v.s6; idx += get_max_sub_group_size(); out_ptr[idx] = v.s7; idx += get_max_sub_group_size();
#define BLOCK_WRITE_IMPL_16 BLOCK_WRITE_IMPL_8 out_ptr[idx] = v.s8; idx += get_max_sub_group_size(); out_ptr[idx] = v.s9; idx += get_max_sub_group_size(); out_ptr[idx] = v.sa; idx += get_max_sub_group_size(); out_ptr[idx] = v.sb; idx += get_max_sub_group_size(); out_ptr[idx] = v.sc; idx += get_max_sub_group_size(); out_ptr[idx] = v.sd; idx += get_max_sub_group_size(); out_ptr[idx] = v.se; idx += get_max_sub_group_size(); out_ptr[idx] = v.sf; idx += get_max_sub_group_size();
#define BLOCK_WRITE_IMPL(vec_size) CAT(BLOCK_WRITE_IMPL_, vec_size)
#define BLOCK_WRITE_FUNC_NAME(type_size, vec_size) MAKE_VECTOR_TYPE(BLOCK_WRITE_FUNC(type_size), vec_size)
#define DECLARE_BLOCK_WRITE_EMULATION(type_size, vec_size) inline void BLOCK_WRITE_FUNC_NAME(type_size, vec_size)(__global BLOCK_WRITE_TYPE(type_size)* out_ptr, MAKE_VECTOR_TYPE(BLOCK_WRITE_TYPE(type_size), vec_size) v) { uint idx = get_sub_group_local_id(); BLOCK_WRITE_IMPL(vec_size) }
#if defined(cl_intel_subgroups)
 #define _sub_group_block_write(ptr, v) intel_sub_group_block_write(ptr, v)
 #define _sub_group_block_write2(ptr, v) intel_sub_group_block_write2(ptr, v)
 #define _sub_group_block_write4(ptr, v) intel_sub_group_block_write4(ptr, v)
 #define _sub_group_block_write8(ptr, v) intel_sub_group_block_write8(ptr, v)
#elif (__OPENCL_C_VERSION__ >= 200)
 DECLARE_BLOCK_WRITE_EMULATION(4, 1)
 DECLARE_BLOCK_WRITE_EMULATION(4, 2)
 DECLARE_BLOCK_WRITE_EMULATION(4, 4)
 DECLARE_BLOCK_WRITE_EMULATION(4, 8)
#endif
#if defined(cl_intel_subgroups_short)
 #define _sub_group_block_write_us(ptr, v) intel_sub_group_block_write_us(ptr, v)
 #define _sub_group_block_write_us2(ptr, v) intel_sub_group_block_write_us2(ptr, v)
 #define _sub_group_block_write_us4(ptr, v) intel_sub_group_block_write_us4(ptr, v)
 #define _sub_group_block_write_us8(ptr, v) intel_sub_group_block_write_us8(ptr, v)
#elif (__OPENCL_C_VERSION__ >= 200)
 DECLARE_BLOCK_WRITE_EMULATION(2, 1)
 DECLARE_BLOCK_WRITE_EMULATION(2, 2)
 DECLARE_BLOCK_WRITE_EMULATION(2, 4)
 DECLARE_BLOCK_WRITE_EMULATION(2, 8)
#endif
#if defined(cl_intel_subgroups_char)
 #define _sub_group_block_write_uc(ptr, v) intel_sub_group_block_write_uc(ptr, v)
 #define _sub_group_block_write_uc2(ptr, v) intel_sub_group_block_write_uc2(ptr, v)
 #define _sub_group_block_write_uc4(ptr, v) intel_sub_group_block_write_uc4(ptr, v)
 #define _sub_group_block_write_uc8(ptr, v) intel_sub_group_block_write_uc8(ptr, v)
 #define _sub_group_block_write_uc16(ptr, v) intel_sub_group_block_write_uc16(ptr, v)
#elif (__OPENCL_C_VERSION__ >= 200)
 DECLARE_BLOCK_WRITE_EMULATION(1, 1)
 DECLARE_BLOCK_WRITE_EMULATION(1, 2)
 DECLARE_BLOCK_WRITE_EMULATION(1, 4)
 DECLARE_BLOCK_WRITE_EMULATION(1, 8)
 DECLARE_BLOCK_WRITE_EMULATION(1, 16)
#endif
#if defined(cl_intel_subgroups_long)
 #define _sub_group_block_write_ul(ptr, v) intel_sub_group_block_write_ul(ptr, v)
 #define _sub_group_block_write_ul2(ptr, v) intel_sub_group_block_write_ul2(ptr, v)
 #define _sub_group_block_write_ul4(ptr, v) intel_sub_group_block_write_ul4(ptr, v)
 #define _sub_group_block_write_ul8(ptr, v) intel_sub_group_block_write_ul8(ptr, v)
#elif (__OPENCL_C_VERSION__ >= 200)
 DECLARE_BLOCK_WRITE_EMULATION(8, 1)
 DECLARE_BLOCK_WRITE_EMULATION(8, 2)
 DECLARE_BLOCK_WRITE_EMULATION(8, 4)
 DECLARE_BLOCK_WRITE_EMULATION(8, 8)
#endif

#ifdef cl_intel_subgroups
#define _sub_group_shuffle(v, c) intel_sub_group_shuffle(v, c)
#define _sub_group_shuffle_up(c, n, d) intel_sub_group_shuffle_up(c, n, d)
#define _sub_group_shuffle_down(c, n, d) intel_sub_group_shuffle_down(c, n, d)
#elif (__OPENCL_C_VERSION__ >= 200)
#define DECLARE_SUB_GROUP_SHUFFLE1(type, cast_type) inline type _sub_group_shuffle(type v, uint c) __attribute__((overloadable)) { return AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v), c)); }
#define DECLARE_SUB_GROUP_SHUFFLE2(type, cast_type) inline CAT(type, 2) _sub_group_shuffle(CAT(type, 2) v, uint c) __attribute__((overloadable)) { return (CAT(type, 2))( AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s0), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s1), c))); }
#define DECLARE_SUB_GROUP_SHUFFLE4(type, cast_type) inline CAT(type, 4) _sub_group_shuffle(CAT(type, 4) v, uint c) __attribute__((overloadable)) { return (CAT(type, 4))( AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s0), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s1), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s2), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s3), c))); }
#define DECLARE_SUB_GROUP_SHUFFLE8(type, cast_type) inline CAT(type, 8) _sub_group_shuffle(CAT(type, 8) v, uint c) __attribute__((overloadable)) { return (CAT(type, 8))( AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s0), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s1), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s2), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s3), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s4), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s5), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s6), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s7), c))); }
#define DECLARE_SUB_GROUP_SHUFFLE16(type, cast_type) inline CAT(type, 16) _sub_group_shuffle(CAT(type, 16) v, uint c) __attribute__((overloadable)) { return (CAT(type, 16))( AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s0), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s1), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s2), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s3), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s4), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s5), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s6), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s7), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s8), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.s9), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.sa), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.sb), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.sc), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.sd), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.se), c)), AS_TYPE(type, sub_group_broadcast(AS_TYPE(cast_type, v.sf), c))); }
#define DECLARE_SUB_GROUP_SHUFFLE(type) DECLARE_SUB_GROUP_SHUFFLE1(type, type) DECLARE_SUB_GROUP_SHUFFLE2(type, type) DECLARE_SUB_GROUP_SHUFFLE4(type, type) DECLARE_SUB_GROUP_SHUFFLE8(type, type) DECLARE_SUB_GROUP_SHUFFLE16(type, type)
#define DECLARE_SUB_GROUP_SHUFFLE_CASTED(type, cast_type) DECLARE_SUB_GROUP_SHUFFLE1(type, cast_type) DECLARE_SUB_GROUP_SHUFFLE2(type, cast_type) DECLARE_SUB_GROUP_SHUFFLE4(type, cast_type) DECLARE_SUB_GROUP_SHUFFLE8(type, cast_type) DECLARE_SUB_GROUP_SHUFFLE16(type, cast_type)
DECLARE_SUB_GROUP_SHUFFLE(int)
DECLARE_SUB_GROUP_SHUFFLE(uint)
DECLARE_SUB_GROUP_SHUFFLE(float)
#if defined(cl_khr_fp16)
 DECLARE_SUB_GROUP_SHUFFLE(half)
 DECLARE_SUB_GROUP_SHUFFLE_CASTED(short, half)
 DECLARE_SUB_GROUP_SHUFFLE_CASTED(ushort, half)
#endif
#endif

typedef struct half5 { half s0; half s1; half s2; half s3; half s4; } half5;
typedef struct half6 { half s0; half s1; half s2; half s3; half s4; half s5; } half6;
typedef struct half7 { half s0; half s1; half s2; half s3; half s4; half s5; half s6; } half7;
typedef struct half9 { half s0; half s1; half s2; half s3; half s4; half s5; half s6; half s7;
 half s8; } half9;
typedef struct half10 { half s0; half s1; half s2; half s3; half s4; half s5; half s6; half s7;
 half s8; half s9; } half10;
typedef struct half11 { half s0; half s1; half s2; half s3; half s4; half s5; half s6; half s7;
 half s8; half s9; half sa; } half11;
typedef struct half12 { half s0; half s1; half s2; half s3; half s4; half s5; half s6; half s7;
 half s8; half s9; half sa; half sb;} half12;
typedef struct half13 { half s0; half s1; half s2; half s3; half s4; half s5; half s6; half s7;
 half s8; half s9; half sa; half sb; half sc;} half13;
typedef struct half14 { half s0; half s1; half s2; half s3; half s4; half s5; half s6; half s7;
 half s8; half s9; half sa; half sb; half sc; half se;} half14;
typedef struct half15 { half s0; half s1; half s2; half s3; half s4; half s5; half s6; half s7;
 half s8; half s9; half sa; half sb; half sc; half se; half sf;} half15;
typedef struct half0 { half s0; } half0;
typedef struct float5 { float s0; float s1; float s2; float s3; float s4; } float5;
typedef struct float6 { float s0; float s1; float s2; float s3; float s4; float s5; } float6;
typedef struct float7 { float s0; float s1; float s2; float s3; float s4; float s5; float s6; } float7;
typedef struct float9 { float s0; float s1; float s2; float s3; float s4; float s5; float s6; float s7; float s8; } float9;
typedef struct float10 { float s0; float s1; float s2; float s3; float s4; float s5;
 float s6; float s7; float s8; float s9;} float10;
typedef struct float11 { float s0; float s1; float s2; float s3; float s4; float s5;
 float s6; float s7; float s8; float s9; float sa;} float11;
typedef struct float12 { float s0; float s1; float s2; float s3; float s4; float s5;
 float s6; float s7; float s8; float s9; float sa; float sb; } float12;
typedef struct float13 { float s0; float s1; float s2; float s3; float s4; float s5;
 float s6; float s7; float s8; float s9; float sa; float sb; float sc;} float13;
typedef struct float14 { float s0; float s1; float s2; float s3; float s4; float s5;
 float s6; float s7; float s8; float s9; float sa; float sb; float sc; float sd; } float14;
typedef struct float15 { float s0; float s1; float s2; float s3; float s4; float s5;
 float s6; float s7; float s8; float s9; float sa; float sb; float sc; float sd; float se; } float15;
typedef struct float0 { float s0; } float0;

//====================================================
// Kernel template: sdpa_micro_prefill 
// Kernel name: sdpa_micro_prefill_8660372428234100028_0_0__sa
#define KERNEL(name) __kernel void sdpa_micro_prefill_8660372428234100028_0_0__sa
#define KERNEL_ID sdpa_micro_prefill_8660372428234100028_0_0__sa
#define FUNC(name)  _##name##_sdpa_micro_prefill_8660372428234100028_0_0__sa
#define FUNC_CALL(name)  _##name##_sdpa_micro_prefill_8660372428234100028_0_0__sa
#define CONST_ARRAY_DECL(name) __constant size_t  _##name##_sdpa_micro_prefill_8660372428234100028_0_0__sa []
#define CONST_ARRAY_REF(name)  _##name##_sdpa_micro_prefill_8660372428234100028_0_0__sa
#define FP64_SUPPORTED 0
#define FP16_SUPPORTED 1
#define FP16_UNIT_USED 1
#define INT8_UNIT_USED 0
#define INT32_UNIT_USED 0
#define INT64_UNIT_USED 0
#define UINT8_UNIT_USED 0
#define UINT32_UNIT_USED 0
#define UNIT_TYPE half
#define UNIT_VAL_MAX HALF_MAX
#define UNIT_VAL_MIN -UNIT_VAL_MAX
#define UNIT_VAL_ONE 1.0h
#define UNIT_VAL_ZERO 0.0h
#define TO_UNIT_TYPE(v) convert_half(v)
#define TO_UNIT_TYPE_SAT(v) convert_half(v)
#define AS_UNIT_TYPE(v) as_half(v)
#define UNIT_MAX_FUNC fmax
#define UNIT_MIN_FUNC fmin
#define UNIT_ABS_FUNC fabs
#define UNIT_TYPE_SIZE 2
#define UNIT_IS_FP 1
#define NL_M as_float(0x0)/*0.000000e+00*/
#define NL_N as_float(0x0)/*0.000000e+00*/
#define ACTIVATION_FUNC_TYPE half
#define ACTIVATION_FUNC_VAL_MAX HALF_MAX
#define ACTIVATION_FUNC_VAL_MIN -ACTIVATION_FUNC_VAL_MAX
#define ACTIVATION_FUNC_VAL_ONE 1.0h
#define ACTIVATION_FUNC_VAL_ZERO 0.0h
#define TO_ACTIVATION_FUNC_TYPE(v) convert_half(v)
#define TO_ACTIVATION_FUNC_TYPE_SAT(v) convert_half(v)
#define AS_ACTIVATION_FUNC_TYPE(v) as_half(v)
#define ACTIVATION_FUNC_MAX_FUNC fmax
#define ACTIVATION_FUNC_MIN_FUNC fmin
#define ACTIVATION_FUNC_ABS_FUNC fabs
#define ACTIVATION_FUNC_TYPE_SIZE 2
#define ACTIVATION_FUNC_IS_FP 1
#define ACTIVATION_PARAMS NL_M, NL_N
#define ACTIVATION_FUNC(input, m, n) input
#define ACTIVATION(input, params) ACTIVATION_FUNC(input, params)
#define INPUT0_SIZE_X 128
#define INPUT0_SIZE_Y (shape_info[6] )
#define INPUT0_SIZE_Z 1
#define INPUT0_SIZE_W 1
#define INPUT0_SIZE_U 1
#define INPUT0_SIZE_V 1
#define INPUT0_FEATURE_NUM 28
#define INPUT0_BATCH_NUM (shape_info[0] )
#define INPUT0_PAD_BEFORE_SIZE_X 0
#define INPUT0_PAD_BEFORE_SIZE_Y 0
#define INPUT0_PAD_BEFORE_SIZE_Z 0
#define INPUT0_PAD_BEFORE_SIZE_W 0
#define INPUT0_PAD_BEFORE_SIZE_U 0
#define INPUT0_PAD_BEFORE_SIZE_V 0
#define INPUT0_PAD_BEFORE_FEATURE_NUM 0
#define INPUT0_PAD_BEFORE_BATCH_NUM 0
#define INPUT0_PAD_AFTER_SIZE_X 0
#define INPUT0_PAD_AFTER_SIZE_Y 0
#define INPUT0_PAD_AFTER_SIZE_Z 0
#define INPUT0_PAD_AFTER_SIZE_W 0
#define INPUT0_PAD_AFTER_SIZE_U 0
#define INPUT0_PAD_AFTER_SIZE_V 0
#define INPUT0_PAD_AFTER_FEATURE_NUM 0
#define INPUT0_PAD_AFTER_BATCH_NUM 0
#define INPUT0_X_PITCH 1
#define INPUT0_Y_PITCH 128
#define INPUT0_Z_PITCH (128*(shape_info[6]  + 0))
#define INPUT0_W_PITCH (128*(shape_info[6]  + 0)*1)
#define INPUT0_U_PITCH (128*(shape_info[6]  + 0)*1*1)
#define INPUT0_V_PITCH (128*(shape_info[6]  + 0)*1*1*1)
#define INPUT0_FEATURE_PITCH (128*(shape_info[6]  + 0)*1*1*1*1)
#define INPUT0_BATCH_PITCH (128*(shape_info[6]  + 0)*1*1*1*1*28)
#define INPUT0_GET_INDEX_SAFE(b, f, y, x) GET_DATA_INDEX_SAFE(INPUT0, b, f, y, x)
#define INPUT0_GET_INDEX(b, f, y, x) GET_DATA_INDEX(INPUT0, b, f, y, x)
#define INPUT0_GET_INDEX_RAW(b, f, y, x) GET_DATA_INDEX_RAW(INPUT0, b, f, y, x)
#define INPUT0_VIEW_OFFSET 0
#define INPUT0_LENGTH 0
#define INPUT0_DIMS 4
#define INPUT0_SIMPLE 1
#define INPUT0_GROUPED 0
#define INPUT0_LAYOUT_BFYX 1
#define INPUT0_TYPE half
#define INPUT0_VAL_MAX HALF_MAX
#define INPUT0_VAL_MIN -INPUT0_VAL_MAX
#define INPUT0_VAL_ONE 1.0h
#define INPUT0_VAL_ZERO 0.0h
#define TO_INPUT0_TYPE(v) convert_half(v)
#define TO_INPUT0_TYPE_SAT(v) convert_half(v)
#define AS_INPUT0_TYPE(v) as_half(v)
#define INPUT0_MAX_FUNC fmax
#define INPUT0_MIN_FUNC fmin
#define INPUT0_ABS_FUNC fabs
#define INPUT0_TYPE_SIZE 2
#define INPUT0_IS_FP 1
#define INPUT0_OFFSET ((INPUT0_X_PITCH*INPUT0_PAD_BEFORE_SIZE_X) + (INPUT0_Y_PITCH*INPUT0_PAD_BEFORE_SIZE_Y) + (INPUT0_Z_PITCH*INPUT0_PAD_BEFORE_SIZE_Z) + (INPUT0_W_PITCH*INPUT0_PAD_BEFORE_SIZE_W) + (INPUT0_FEATURE_PITCH*INPUT0_PAD_BEFORE_FEATURE_NUM) + (INPUT0_BATCH_PITCH*INPUT0_PAD_BEFORE_BATCH_NUM))
#define INPUT0_PAD_BEFORE (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT0_PAD_AFTER (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT1_SIZE_X 128
#define INPUT1_SIZE_Y (shape_info[14] )
#define INPUT1_SIZE_Z 1
#define INPUT1_SIZE_W 1
#define INPUT1_SIZE_U 1
#define INPUT1_SIZE_V 1
#define INPUT1_FEATURE_NUM 4
#define INPUT1_BATCH_NUM (shape_info[8] )
#define INPUT1_PAD_BEFORE_SIZE_X 0
#define INPUT1_PAD_BEFORE_SIZE_Y (shape_info[16])
#define INPUT1_PAD_BEFORE_SIZE_Z 0
#define INPUT1_PAD_BEFORE_SIZE_W 0
#define INPUT1_PAD_BEFORE_SIZE_U 0
#define INPUT1_PAD_BEFORE_SIZE_V 0
#define INPUT1_PAD_BEFORE_FEATURE_NUM 0
#define INPUT1_PAD_BEFORE_BATCH_NUM 0
#define INPUT1_PAD_AFTER_SIZE_X 0
#define INPUT1_PAD_AFTER_SIZE_Y (shape_info[17])
#define INPUT1_PAD_AFTER_SIZE_Z 0
#define INPUT1_PAD_AFTER_SIZE_W 0
#define INPUT1_PAD_AFTER_SIZE_U 0
#define INPUT1_PAD_AFTER_SIZE_V 0
#define INPUT1_PAD_AFTER_FEATURE_NUM 0
#define INPUT1_PAD_AFTER_BATCH_NUM 0
#define INPUT1_X_PITCH 1
#define INPUT1_Y_PITCH 128
#define INPUT1_Z_PITCH (128*(shape_info[14]  + (shape_info[16] + shape_info[17])))
#define INPUT1_W_PITCH (128*(shape_info[14]  + (shape_info[16] + shape_info[17]))*1)
#define INPUT1_U_PITCH (128*(shape_info[14]  + (shape_info[16] + shape_info[17]))*1*1)
#define INPUT1_V_PITCH (128*(shape_info[14]  + (shape_info[16] + shape_info[17]))*1*1*1)
#define INPUT1_FEATURE_PITCH (128*(shape_info[14]  + (shape_info[16] + shape_info[17]))*1*1*1*1)
#define INPUT1_BATCH_PITCH (128*(shape_info[14]  + (shape_info[16] + shape_info[17]))*1*1*1*1*4)
#define INPUT1_GET_INDEX_SAFE(b, f, y, x) GET_DATA_INDEX_SAFE(INPUT1, b, f, y, x)
#define INPUT1_GET_INDEX(b, f, y, x) GET_DATA_INDEX(INPUT1, b, f, y, x)
#define INPUT1_GET_INDEX_RAW(b, f, y, x) GET_DATA_INDEX_RAW(INPUT1, b, f, y, x)
#define INPUT1_VIEW_OFFSET 0
#define INPUT1_LENGTH 0
#define INPUT1_DIMS 4
#define INPUT1_SIMPLE 1
#define INPUT1_GROUPED 0
#define INPUT1_LAYOUT_BFYX 1
#define INPUT1_TYPE half
#define INPUT1_VAL_MAX HALF_MAX
#define INPUT1_VAL_MIN -INPUT1_VAL_MAX
#define INPUT1_VAL_ONE 1.0h
#define INPUT1_VAL_ZERO 0.0h
#define TO_INPUT1_TYPE(v) convert_half(v)
#define TO_INPUT1_TYPE_SAT(v) convert_half(v)
#define AS_INPUT1_TYPE(v) as_half(v)
#define INPUT1_MAX_FUNC fmax
#define INPUT1_MIN_FUNC fmin
#define INPUT1_ABS_FUNC fabs
#define INPUT1_TYPE_SIZE 2
#define INPUT1_IS_FP 1
#define INPUT1_OFFSET ((INPUT1_X_PITCH*INPUT1_PAD_BEFORE_SIZE_X) + (INPUT1_Y_PITCH*INPUT1_PAD_BEFORE_SIZE_Y) + (INPUT1_Z_PITCH*INPUT1_PAD_BEFORE_SIZE_Z) + (INPUT1_W_PITCH*INPUT1_PAD_BEFORE_SIZE_W) + (INPUT1_FEATURE_PITCH*INPUT1_PAD_BEFORE_FEATURE_NUM) + (INPUT1_BATCH_PITCH*INPUT1_PAD_BEFORE_BATCH_NUM))
#define INPUT1_PAD_BEFORE (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT1_PAD_AFTER (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT2_SIZE_X 128
#define INPUT2_SIZE_Y (shape_info[24] )
#define INPUT2_SIZE_Z 1
#define INPUT2_SIZE_W 1
#define INPUT2_SIZE_U 1
#define INPUT2_SIZE_V 1
#define INPUT2_FEATURE_NUM 4
#define INPUT2_BATCH_NUM (shape_info[18] )
#define INPUT2_PAD_BEFORE_SIZE_X 0
#define INPUT2_PAD_BEFORE_SIZE_Y (shape_info[26])
#define INPUT2_PAD_BEFORE_SIZE_Z 0
#define INPUT2_PAD_BEFORE_SIZE_W 0
#define INPUT2_PAD_BEFORE_SIZE_U 0
#define INPUT2_PAD_BEFORE_SIZE_V 0
#define INPUT2_PAD_BEFORE_FEATURE_NUM 0
#define INPUT2_PAD_BEFORE_BATCH_NUM 0
#define INPUT2_PAD_AFTER_SIZE_X 0
#define INPUT2_PAD_AFTER_SIZE_Y (shape_info[27])
#define INPUT2_PAD_AFTER_SIZE_Z 0
#define INPUT2_PAD_AFTER_SIZE_W 0
#define INPUT2_PAD_AFTER_SIZE_U 0
#define INPUT2_PAD_AFTER_SIZE_V 0
#define INPUT2_PAD_AFTER_FEATURE_NUM 0
#define INPUT2_PAD_AFTER_BATCH_NUM 0
#define INPUT2_X_PITCH 1
#define INPUT2_Y_PITCH 128
#define INPUT2_Z_PITCH (128*(shape_info[24]  + (shape_info[26] + shape_info[27])))
#define INPUT2_W_PITCH (128*(shape_info[24]  + (shape_info[26] + shape_info[27]))*1)
#define INPUT2_U_PITCH (128*(shape_info[24]  + (shape_info[26] + shape_info[27]))*1*1)
#define INPUT2_V_PITCH (128*(shape_info[24]  + (shape_info[26] + shape_info[27]))*1*1*1)
#define INPUT2_FEATURE_PITCH (128*(shape_info[24]  + (shape_info[26] + shape_info[27]))*1*1*1*1)
#define INPUT2_BATCH_PITCH (128*(shape_info[24]  + (shape_info[26] + shape_info[27]))*1*1*1*1*4)
#define INPUT2_GET_INDEX_SAFE(b, f, y, x) GET_DATA_INDEX_SAFE(INPUT2, b, f, y, x)
#define INPUT2_GET_INDEX(b, f, y, x) GET_DATA_INDEX(INPUT2, b, f, y, x)
#define INPUT2_GET_INDEX_RAW(b, f, y, x) GET_DATA_INDEX_RAW(INPUT2, b, f, y, x)
#define INPUT2_VIEW_OFFSET 0
#define INPUT2_LENGTH 0
#define INPUT2_DIMS 4
#define INPUT2_SIMPLE 1
#define INPUT2_GROUPED 0
#define INPUT2_LAYOUT_BFYX 1
#define INPUT2_TYPE half
#define INPUT2_VAL_MAX HALF_MAX
#define INPUT2_VAL_MIN -INPUT2_VAL_MAX
#define INPUT2_VAL_ONE 1.0h
#define INPUT2_VAL_ZERO 0.0h
#define TO_INPUT2_TYPE(v) convert_half(v)
#define TO_INPUT2_TYPE_SAT(v) convert_half(v)
#define AS_INPUT2_TYPE(v) as_half(v)
#define INPUT2_MAX_FUNC fmax
#define INPUT2_MIN_FUNC fmin
#define INPUT2_ABS_FUNC fabs
#define INPUT2_TYPE_SIZE 2
#define INPUT2_IS_FP 1
#define INPUT2_OFFSET ((INPUT2_X_PITCH*INPUT2_PAD_BEFORE_SIZE_X) + (INPUT2_Y_PITCH*INPUT2_PAD_BEFORE_SIZE_Y) + (INPUT2_Z_PITCH*INPUT2_PAD_BEFORE_SIZE_Z) + (INPUT2_W_PITCH*INPUT2_PAD_BEFORE_SIZE_W) + (INPUT2_FEATURE_PITCH*INPUT2_PAD_BEFORE_FEATURE_NUM) + (INPUT2_BATCH_PITCH*INPUT2_PAD_BEFORE_BATCH_NUM))
#define INPUT2_PAD_BEFORE (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT2_PAD_AFTER (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT3_SIZE_X (shape_info[35] )
#define INPUT3_SIZE_Y (shape_info[34] )
#define INPUT3_SIZE_Z 1
#define INPUT3_SIZE_W 1
#define INPUT3_SIZE_U 1
#define INPUT3_SIZE_V 1
#define INPUT3_FEATURE_NUM (shape_info[29] )
#define INPUT3_BATCH_NUM (shape_info[28] )
#define INPUT3_PAD_BEFORE_SIZE_X 0
#define INPUT3_PAD_BEFORE_SIZE_Y 0
#define INPUT3_PAD_BEFORE_SIZE_Z 0
#define INPUT3_PAD_BEFORE_SIZE_W 0
#define INPUT3_PAD_BEFORE_SIZE_U 0
#define INPUT3_PAD_BEFORE_SIZE_V 0
#define INPUT3_PAD_BEFORE_FEATURE_NUM 0
#define INPUT3_PAD_BEFORE_BATCH_NUM 0
#define INPUT3_PAD_AFTER_SIZE_X 0
#define INPUT3_PAD_AFTER_SIZE_Y 0
#define INPUT3_PAD_AFTER_SIZE_Z 0
#define INPUT3_PAD_AFTER_SIZE_W 0
#define INPUT3_PAD_AFTER_SIZE_U 0
#define INPUT3_PAD_AFTER_SIZE_V 0
#define INPUT3_PAD_AFTER_FEATURE_NUM 0
#define INPUT3_PAD_AFTER_BATCH_NUM 0
#define INPUT3_X_PITCH 1
#define INPUT3_Y_PITCH (shape_info[35]  + 0)
#define INPUT3_Z_PITCH ((shape_info[35]  + 0)*(shape_info[34]  + 0))
#define INPUT3_W_PITCH ((shape_info[35]  + 0)*(shape_info[34]  + 0)*1)
#define INPUT3_U_PITCH ((shape_info[35]  + 0)*(shape_info[34]  + 0)*1*1)
#define INPUT3_V_PITCH ((shape_info[35]  + 0)*(shape_info[34]  + 0)*1*1*1)
#define INPUT3_FEATURE_PITCH ((shape_info[35]  + 0)*(shape_info[34]  + 0)*1*1*1*1)
#define INPUT3_BATCH_PITCH ((shape_info[35]  + 0)*(shape_info[34]  + 0)*1*1*1*1*(shape_info[29]  + 0))
#define INPUT3_GET_INDEX_SAFE(b, f, y, x) GET_DATA_INDEX_SAFE(INPUT3, b, f, y, x)
#define INPUT3_GET_INDEX(b, f, y, x) GET_DATA_INDEX(INPUT3, b, f, y, x)
#define INPUT3_GET_INDEX_RAW(b, f, y, x) GET_DATA_INDEX_RAW(INPUT3, b, f, y, x)
#define INPUT3_VIEW_OFFSET 0
#define INPUT3_LENGTH 0
#define INPUT3_DIMS 4
#define INPUT3_SIMPLE 1
#define INPUT3_GROUPED 0
#define INPUT3_LAYOUT_BFYX 1
#define INPUT3_TYPE half
#define INPUT3_VAL_MAX HALF_MAX
#define INPUT3_VAL_MIN -INPUT3_VAL_MAX
#define INPUT3_VAL_ONE 1.0h
#define INPUT3_VAL_ZERO 0.0h
#define TO_INPUT3_TYPE(v) convert_half(v)
#define TO_INPUT3_TYPE_SAT(v) convert_half(v)
#define AS_INPUT3_TYPE(v) as_half(v)
#define INPUT3_MAX_FUNC fmax
#define INPUT3_MIN_FUNC fmin
#define INPUT3_ABS_FUNC fabs
#define INPUT3_TYPE_SIZE 2
#define INPUT3_IS_FP 1
#define INPUT3_OFFSET ((INPUT3_X_PITCH*INPUT3_PAD_BEFORE_SIZE_X) + (INPUT3_Y_PITCH*INPUT3_PAD_BEFORE_SIZE_Y) + (INPUT3_Z_PITCH*INPUT3_PAD_BEFORE_SIZE_Z) + (INPUT3_W_PITCH*INPUT3_PAD_BEFORE_SIZE_W) + (INPUT3_FEATURE_PITCH*INPUT3_PAD_BEFORE_FEATURE_NUM) + (INPUT3_BATCH_PITCH*INPUT3_PAD_BEFORE_BATCH_NUM))
#define INPUT3_PAD_BEFORE (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT3_PAD_AFTER (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define OUTPUT_SIZE_X 128
#define OUTPUT_SIZE_Y (shape_info[50] )
#define OUTPUT_SIZE_Z 1
#define OUTPUT_SIZE_W 1
#define OUTPUT_SIZE_U 1
#define OUTPUT_SIZE_V 1
#define OUTPUT_FEATURE_NUM 28
#define OUTPUT_BATCH_NUM (shape_info[44] )
#define OUTPUT_PAD_BEFORE_SIZE_X 0
#define OUTPUT_PAD_BEFORE_SIZE_Y 0
#define OUTPUT_PAD_BEFORE_SIZE_Z 0
#define OUTPUT_PAD_BEFORE_SIZE_W 0
#define OUTPUT_PAD_BEFORE_SIZE_U 0
#define OUTPUT_PAD_BEFORE_SIZE_V 0
#define OUTPUT_PAD_BEFORE_FEATURE_NUM 0
#define OUTPUT_PAD_BEFORE_BATCH_NUM 0
#define OUTPUT_PAD_AFTER_SIZE_X 0
#define OUTPUT_PAD_AFTER_SIZE_Y 0
#define OUTPUT_PAD_AFTER_SIZE_Z 0
#define OUTPUT_PAD_AFTER_SIZE_W 0
#define OUTPUT_PAD_AFTER_SIZE_U 0
#define OUTPUT_PAD_AFTER_SIZE_V 0
#define OUTPUT_PAD_AFTER_FEATURE_NUM 0
#define OUTPUT_PAD_AFTER_BATCH_NUM 0
#define OUTPUT_X_PITCH 1
#define OUTPUT_Y_PITCH 128
#define OUTPUT_Z_PITCH (128*(shape_info[50]  + 0))
#define OUTPUT_W_PITCH (128*(shape_info[50]  + 0)*1)
#define OUTPUT_U_PITCH (128*(shape_info[50]  + 0)*1*1)
#define OUTPUT_V_PITCH (128*(shape_info[50]  + 0)*1*1*1)
#define OUTPUT_FEATURE_PITCH (128*(shape_info[50]  + 0)*1*1*1*1)
#define OUTPUT_BATCH_PITCH (128*(shape_info[50]  + 0)*1*1*1*1*28)
#define OUTPUT_GET_INDEX_SAFE(b, f, y, x) GET_DATA_INDEX_SAFE(OUTPUT, b, f, y, x)
#define OUTPUT_GET_INDEX(b, f, y, x) GET_DATA_INDEX(OUTPUT, b, f, y, x)
#define OUTPUT_GET_INDEX_RAW(b, f, y, x) GET_DATA_INDEX_RAW(OUTPUT, b, f, y, x)
#define OUTPUT_VIEW_OFFSET 0
#define OUTPUT_LENGTH 0
#define OUTPUT_DIMS 4
#define OUTPUT_SIMPLE 1
#define OUTPUT_GROUPED 0
#define OUTPUT_LAYOUT_BFYX 1
#define OUTPUT_TYPE half
#define OUTPUT_VAL_MAX HALF_MAX
#define OUTPUT_VAL_MIN -OUTPUT_VAL_MAX
#define OUTPUT_VAL_ONE 1.0h
#define OUTPUT_VAL_ZERO 0.0h
#define TO_OUTPUT_TYPE(v) convert_half(v)
#define TO_OUTPUT_TYPE_SAT(v) convert_half(v)
#define AS_OUTPUT_TYPE(v) as_half(v)
#define OUTPUT_MAX_FUNC fmax
#define OUTPUT_MIN_FUNC fmin
#define OUTPUT_ABS_FUNC fabs
#define OUTPUT_TYPE_SIZE 2
#define OUTPUT_IS_FP 1
#define OUTPUT_OFFSET ((OUTPUT_X_PITCH*OUTPUT_PAD_BEFORE_SIZE_X) + (OUTPUT_Y_PITCH*OUTPUT_PAD_BEFORE_SIZE_Y) + (OUTPUT_Z_PITCH*OUTPUT_PAD_BEFORE_SIZE_Z) + (OUTPUT_W_PITCH*OUTPUT_PAD_BEFORE_SIZE_W) + (OUTPUT_FEATURE_PITCH*OUTPUT_PAD_BEFORE_FEATURE_NUM) + (OUTPUT_BATCH_PITCH*OUTPUT_PAD_BEFORE_BATCH_NUM))
#define OUTPUT_PAD_BEFORE (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define OUTPUT_PAD_AFTER (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define IS_DYNAMIC 1
#define OPTIONAL_SHAPE_INFO_ARG __global const int* shape_info,
#define OPTIONAL_SHAPE_INFO_TENSOR shape_info,
#define LayerID indirectsdpa:__module.model.layers.0.self_attn/aten::scaled_dot_product_attention/ScaledDotProductAttention
#define D_MAX 128
#define SUBGROUP_SIZE 8
#define INVERT_SCALE 0
#define SCALE_DATA_T half
#define WITH_ATTN_MASK 1
#define WITH_SCALE 0
#define Q_ALIGN 128
#define K_ALIGN 128
#define V_ALIGN 128
#define A_ALIGN 128
#define TRANSPOSE_K 0
#define REMAINDER_K 1
#define KV_GROUP_SIZE 7
#define BLOCK_Q 1
#define REMAINDER_Q 1
#define QRY_S0 INPUT0_BATCH_PITCH
#define QRY_S1 INPUT0_FEATURE_PITCH
#define QRY_S2 INPUT0_Y_PITCH
#define QRY_S3 INPUT0_X_PITCH
#define KEY_S0 INPUT1_BATCH_PITCH
#define KEY_S1 INPUT1_FEATURE_PITCH
#define KEY_S2 INPUT1_Y_PITCH
#define KEY_S3 INPUT1_X_PITCH
#define VAL_S0 INPUT2_BATCH_PITCH
#define VAL_S1 INPUT2_FEATURE_PITCH
#define VAL_S2 INPUT2_Y_PITCH
#define VAL_S3 INPUT2_X_PITCH
#define DST_S0 OUTPUT_BATCH_PITCH
#define DST_S1 OUTPUT_FEATURE_PITCH
#define DST_S2 OUTPUT_Y_PITCH
#define DST_S3 OUTPUT_X_PITCH
#define QRY_B0 1
#define QRY_SB0 1
#define QRY_B1 1
#define QRY_SB1 1
#define QRY_B2 1
#define QRY_SB2 1
#define QRY_B3 1
#define QRY_SB3 1
#define KEY_B0 1
#define KEY_SB0 1
#define KEY_B1 1
#define KEY_SB1 1
#define KEY_B2 1
#define KEY_SB2 1
#define KEY_B3 1
#define KEY_SB3 1
#define VAL_B0 1
#define VAL_SB0 1
#define VAL_B1 1
#define VAL_SB1 1
#define VAL_B2 1
#define VAL_SB2 1
#define VAL_B3 1
#define VAL_SB3 1
#define DST_B0 1
#define DST_SB0 1
#define DST_B1 1
#define DST_SB1 1
#define DST_B2 1
#define DST_SB2 1
#define DST_B3 1
#define DST_SB3 1

#ifndef MICRO_DECL_float_tile_16x16_blocked_8x8
#define MICRO_DECL_float_tile_16x16_blocked_8x8
typedef struct {
    float8 x[4];
} float_tile_16x16_blocked_8x8;
DECLARE_2D_TILE_OPS(float_tile_16x16_blocked_8x8,float,8,8,8,2,2)
#endif
#define ugemm_kq_sg_tile_m 16
#define ugemm_kq_sg_tile_n 16
#define ugemm_kq_wg_tile_m 256
#define ugemm_kq_wg_tile_n 32
#define ugemm_kq_sg_per_wg_m 16
#define ugemm_kq_sg_per_wg_n 2
#define ugemm_kq_sg_per_wg_k 1
#define ugemm_kq_slm_size 16384
#define ugemm_kq_barrier_count  1
#define ugemm_kq_systolic  1
typedef float_tile_16x16_blocked_8x8 ugemm_kq_c_type;
#define ugemm_kq_c_type_block0 8
#define ugemm_kq_c_type_nblock0 2
#define ugemm_kq_c_type_block1 8
#define ugemm_kq_c_type_nblock1 2
ugemm_kq_c_type ugemm_kq(const global half* a, int lda, const local half* b, int ldb, int m, int n, int k, int i0, int j0, int h0, int local_id_m, int local_id_n, const local char* slm) {
    ugemm_kq_c_type c;
    __asm__ volatile("{\n"
            ".implicit_PSEUDO_INPUT %0 offset=2560 size=256\n"
            ".implicit_PSEUDO_INPUT %1 offset=3072 size=256\n"
            ".implicit_PSEUDO_INPUT %2 offset=2816 size=256\n"
            ".implicit_PSEUDO_INPUT %3 offset=3328 size=256\n"
            ".decl COPY4 v_type=G type=uq num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY4 offset=264 size=8\n"
            ".decl COPY5 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY5 offset=276 size=4\n"
            ".decl COPY6 v_type=G type=ud num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY6 offset=256 size=4\n"
            ".decl COPY7 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY7 offset=272 size=4\n"
            ".decl COPY8 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY8 offset=284 size=4\n"
            ".decl COPY9 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY9 offset=280 size=4\n"
            ".decl COPY10 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY10 offset=288 size=4\n"
            ".decl COPY11 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY11 offset=296 size=4\n"
            ".decl COPY12 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY12 offset=292 size=4\n"
            ".decl COPY13 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY13 offset=300 size=4\n"
            ".decl COPY14 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY14 offset=308 size=4\n"
            ".decl COPY15 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY15 offset=304 size=4\n"
            ".decl COPY16 v_type=G type=ud num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY16 offset=260 size=4\n"
            "fence_sw\n"
            "mov (M1_NM, 1) COPY4(0,0)<1> %4(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY5(0,0)<1> %5(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY6(0,0)<1> %6(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY7(0,0)<1> %7(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY8(0,0)<1> %8(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY9(0,0)<1> %9(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY10(0,0)<1> %10(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY11(0,0)<1> %11(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY12(0,0)<1> %12(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY13(0,0)<1> %13(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY14(0,0)<1> %14(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY15(0,0)<1> %15(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY16(0,0)<1> %16(0,0)<1;1,0>\n"
            ".decl CLOBBER0 v_type=G type=ud num_elts=2\n"
            ".implicit_PSEUDO_INPUT CLOBBER0 offset=312 size=8\n"
            ".decl CLOBBER1 v_type=G type=ud num_elts=560\n"
            ".implicit_PSEUDO_INPUT CLOBBER1 offset=320 size=2240\n"
            ".decl CLOBBER2 v_type=G type=ud num_elts=448\n"
            ".implicit_PSEUDO_INPUT CLOBBER2 offset=3584 size=1792\n"
            "fence_sw\n"
            "add (M1,1) CLOBBER0(0,0)<1> CLOBBER0(0,0)<0;1,0> 0xcafefade:ud\n"
            "fence_sw\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY4(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY5(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY6(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY7(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY8(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY9(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY10(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY11(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY12(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY13(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY14(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY15(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY16(0,0)<1;1,0>\n"
            "mov (M1,2) CLOBBER0(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(6,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(8,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(10,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(12,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(14,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(16,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(18,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(20,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(22,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(24,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(26,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(28,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(30,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(32,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(34,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(36,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(38,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(40,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(42,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(44,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(46,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(48,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(50,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(52,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(54,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(56,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(58,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(60,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(62,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(64,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(66,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(68,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(6,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(8,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(10,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(12,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(14,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(16,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(18,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(20,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(22,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(24,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(26,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(28,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(30,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(32,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(34,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(36,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(38,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(40,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(42,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(44,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(46,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(48,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(50,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(52,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(54,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %0(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %0(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %0(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %0(6,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %2(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %2(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %2(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %2(6,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %1(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %1(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %1(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %1(6,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %3(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %3(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %3(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %3(6,0)<1> 0xAAAAAAAA:ud\n"
            ".decl DUMMY_DPAS_SRC v_type=G type=ud num_elts=64 alias=<CLOBBER1,0>\n"
            ".decl DUMMY_DPAS_DST v_type=G type=f num_elts=64 alias=<DUMMY_DPAS_SRC,0>\n"
            "dpas.bf.bf.8.1 (M1,8) DUMMY_DPAS_DST.0 V0.0 DUMMY_DPAS_SRC.0 DUMMY_DPAS_SRC(0,0)\n"
            "barrier\n"
            "fence_sw\n"
            "add (M1,1) CLOBBER0(0,0)<1> CLOBBER0(0,0)<0;1,0> 0xfadecafe:ud\n"
            "fence_sw\n"
        "}\n"
    : "=rw"(c.x[0]), "=rw"(c.x[1]), "=rw"(c.x[2]), "=rw"(c.x[3]) : "rw.u"(a), "rw.u"(lda), "rw.u"(b), "rw.u"(ldb), "rw.u"(m), "rw.u"(n), "rw.u"(k), "rw.u"(i0), "rw.u"(j0), "rw.u"(h0), "rw.u"(local_id_m), "rw.u"(local_id_n), "rw.u"(slm));
    // @_u_@0 610900802002c5094070000000000000610900802041417000000000ffffffff01090080000001000000003000000000610000802002e5090030000000000000610000802002050c00310000000000007a00008010c5240c100005018409000069000080208285088408000101000100690000802082a508a4080001010001006100008060060540240900000000000061000080600625404409000000000000690000802081450c84090001040004004001008060022509440c0006240900005b0000806086440944090501a4091000400000802426450c04400006c4080000400000802426650c24400006e4080000410400806006c1202409000184080000490000802206c50c2409000284080000610100806002e50cc40c000000000000610000802006c50cc020000000000000400100806006050804080006c40c0000410000806006812044090001a4080000490000802206850c44090002a4080000610100806002a50c840c000000000000610000802006850c80200000000000004e0100802202010044080002840c0000400000802002450844080002840c0000610000802002050a0020000000000000400000802002650864080002a40c00004001008020026508040a0002640800005b000080608604080408050664092000690000806086050d64090001010001006c0100806086050a040d00011f001f004e00008022020100040d0002440800004000008020024508040d000244080000610000802002050b00200000000000004004008020066508040a0002640800004001008020026508040b0002640800004000008020a1850c240c000110001000400000802426a50c24090006c4080000400000802426c50c44090006e4080000400200802422850c840c0002a40c0000620000802082a50ca40c005110001000620300802082c50cc40c0051100010007000008060820100640c005100010001200000810040000000000000d02d0000410000802002450ca408000184090000690100802082450c440c000103000300690000802081650984090001080008006900008020812509a40900010a000a004e0300802202010044080002440c0000400000802002450844080002440c0000610000802002050a00200000000000004001008020026508040a0002640800004000008020022509240900022408000040010080200245092409000264090000610001802002058145082200000000004e0100802202058204810002a4080000610000802002050a00200000000000004001008020022582040a000224810000690000802082450ca408000101000100680000802082650ca40800011f001f004e0200802202058304810002440c0000610000802002050a0020000000000000400300802002258324810002640c00004001008020022583040a0002248300004100008020824120a408000103000300490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c40200000000000004e0100802202058404810002440c0000610000802002050a0020000000000000400000802002258424810002640c00004001008020022584040a000224840000690000802082450ca408000102000200680000802082650ca40800011e001e004e0200802202058504810002440c0000610000802002050a0020000000000000400300802002258524810002640c00004001008020022585040a0002248500004100008020824120a408000105000500490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c40200000000000004e0100802202058604810002440c0000610000802002050a0020000000000000400000802002258624810002640c00004001008020022586040a0002248600004100008020824120a408000106000600490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c40200000000000004e0100802202058704810002440c0000610000802002050a0020000000000000400000802002258724810002640c00004001008020022587040a0002248700004100008020824120a408000107000700490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c40200000000000004e0100802202058804810002440c0000610000802002050a0020000000000000400000802002258824810002640c00004001008020022588040a00022488000061000080200205894409000000000000400100802082058a0489000500020002610000802006050d0408000000000000400100802082052e040d000500080008610000802002052f24090000000000004001008020820580042f000500020002610000806006650904090000000000006100048020410550000000000000000061000480204105520000000000000000610004802041055400000000000000006100048020410556000000000000000061000480204105580000000000000000610004802041055a0000000000000000610004802041055c0000000000000000610004802041055e00000000000000006100048020410560000000000000000061000480204105620000000000000000610004802041056400000000000000006100048020410566000000000000000061000480204105680000000000000000610004802041056a0000000000000000610004802041056c0000000000000000610004802041056e0000000000000000010900800000010000000030000000007000008060860100640900552000200020000081004000000000000050190000400000806086650964090005c1ffc1ff31460080000014700c8100ff0000300231470080000014720c8200ff0000300231480080000014740c8300ff0000300231490080000014760c8400ff00003002314a0080000014780c8500ff00003002314b00800000147a0c8600ff00003002314c00800000147c0c8700ff00003002314d00800000147e0c8800ff0000300231400080000044300c0d00ee0000380231410080000044380c2e00ee000038024036108020820581048100854000400040001081208225812481000501000100403710802082058204820085400040004000108120822582248200050100010040381080208205830483008540004000400010812082258324830005010001004039108020820584048400854000400040001081208225842484000501000100403a108020820585048500854000400040001081208225852485000501000100403b108020820586048600854000400040001081208225862486000501000100403c108020820587048700854000400040001081208225872487000501000100403d108020820588048800854000400040001081208225882488000501000100403000802082050d040d000500010001403100802082052e042e000500010001612603802002050e057046000000000061000380200205160571460000000000612703802002050f05724600000000006100038020020517057346000000000061280380200205100574460000000000610003802002051805754600000000006129038020020511057646000000000061000380200205190577460000000000612a0380200205120578460000000000610003802002051a0579460000000000612b038020020513057a460000000000610003802002051b057b460000000000612c038020020514057c460000000000610003802002051c057d460000000000612d038020020515057e460000000000610003802002051d057f460000000000610000802042450a0000000000002020319f0080000000000c0a0830000000000100008000000100000000e00000000031440080000000000c8908ee440e3800314f0080000000000c8a08ee44163800314e008000000c0b0c003ee200000000314e0080000000000c0a0830000000000100008000000100000000e00000000070000080608601006409006500000000010900800000010000000000000000002000008100400000000000008004000031460080000014700c8100ff0000300231470080000014720c8200ff0000300231480080000014740c8300ff0000300231490080000014760c8400ff00003002314a0080000014780c8500ff00003002314b00800000147a0c8600ff00003002314c00800000147c0c8700ff00003002314d00800000147e0c8800ff00003002314400800000840e0c2f00ee00003c000100008000420100000000200c000000314500800000841e0c8000ee00003c0031420080000044400c0d00ee0000380231430080000044480c2e00ee00003802400000806086650964090035e0ffe0ff403200802082050d040d000500010001403300802082052e042e0005000100010120008000000100000000000000000059440380a03a0750045001010430040e59400380a03a0758045801010430041659410380a03a0760046001010438040e59400380a03a076804680101043804160134008000000100000000000000000031a00080000044300c0d00ee0000380231910080000044380c2e00ee000038024036108020820581048100854000400040001081208225812481000501000100403710802082058204820085400040004000108120822582248200050100010040381080208205830483008540004000400010812082258324830005010001004039108020820584048400854000400040001081208225842484000501000100403a108020820585048500854000400040001081208225852485000501000100403b108020820586048600854000400040001081208225862486000501000100403c108020820587048700854000400040001081208225872487000501000100403d108020820588048800854000400040001081208225882488000501000100403000802082050d040d000500010001403100802082052e042e0005000100010100008000420100000000303000000059420380a03a0750045001010440041e59450380a03a0758045801010440042659430380a03a0760046001010448041e59440380a03a07680468010104480426612603802002050e0570460000000000613f0380200205160571460000000000612703802002050f05724600000000006100038020020517057346000000000061280380200205100574460000000000610003802002051805754600000000006129038020020511057646000000000061000380200205190577460000000000612a0380200205120578460000000000610003802002051a0579460000000000612b038020020513057a460000000000610003802002051b057b460000000000612c038020020514057c460000000000610003802002051c057d460000000000612d038020020515057e460000000000610003802002051d057f460000000000314d0080000000000c0a0830000000000100008000000100000000e00000000031a40080000000000c8908ee440e3800319f0080000000000c8a08ee44163800314e008000000c0b0c003ee200000000314e0080000000000c0a0830000000000100008000000100000000e000000000200000810040000000000000a0fbffff4000008060866509640900053f003f00610003801040050b0000000010325476610003801040850b0000000098badcfe6901038020810684840b1001020002006900038020810682040b100102000200400213802002068406844482048100006100038020022684248100000000000040011381208226842684440501000100400013802002068206824482048100006100038020022682248100000000000040011381208226822682440501000100400013802002068b06824482a4080000610203802002268b2682440000000000400113812082268b268b440501000100400013802002068d06844482a4080000610003802002268d2684440000000000400113812082268d268d440501000100690000802082450ca408000101000100680000802082650ca40800011f001f00400213802002068f06824482440c0000400203802002268f26824402640c0000400113812082268f268f440501000100400013802002069106844482440c0000400003802002269126844402640c0000400113812082269126914405010001004100008020824120a408000103000300490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069306824482440c0000400003802002269326824402640c000040011381208226932693440501000100400013802002069506844482440c0000400003802002269526844402640c000040011381208226952695440501000100690000802082450ca408000102000200680000802082650ca40800011e001e00400213802002069706824482440c0000400203802002269726824402640c000040011381208226972697440501000100400013802002069906844482440c0000400003802002269926844402640c0000400113812082269926994405010001004100008020824120a408000105000500490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069b06824482440c0000400003802002269b26824402640c0000400113812082269b269b440501000100400013802002069d06844482440c0000400003802002269d26844402640c0000400113812082269d269d4405010001004100008020824120a408000106000600490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069f06824482440c0000400003802002269f26824402640c0000400113812082269f269f44050100010040001380200206a106844482440c000040000380200226a126844402640c000040011381208226a126a14405010001004100008020824120a408000107000700490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c402000000000000040011380200206a306824482440c000040000380200226a326824402640c000040011381208226a326a344050100010040001380200206a506844482440c000040000380200226a526844402640c000040011381208226a526a5440501000100400000801486460c64090005e1ffe1ff680100801081450c440c0001010001004001008014a1450c440c000110001000610000801041050b00000000ffffffff6801008010011130040b0001440c00007000008060860100640900552100210020000081004000000000000010010000400000801486460c64090005e1ffe1ff680100801081450c440c0001010001004001008014a1450c440c000110001000610000801041050b00000000ffffffff6801008010011130040b0001440c0000720002801009040b050b0109050b050b72000280a00a840b850b020a850b850b3146448100001470248200fb000000023147448100001472248b00fb000000023148448100001474248f00fb000000023149448100001476249300fb00000002314a448100001478249700fb00000002314b44810000147a249b00fb00000002314c44810000147c249f00fb00000002012d0080000001000000000000000000314d44810000147e24a300fb00000002314400800000840e0c2f00ee00003c000100008000420100000000200c000000314500800000841e0c8000ee00003c0031420080000044400c0d00ee0000380231430080000044480c2e00ee00003802403200802082050d040d000500010001403300802082052e042e0005000100010120008000000100000000000000000059440380a03a0750045001010430040e59400380a03a0758045801010430041659410380a03a0760046001010438040e594f0380a03a07680468010104380416200000810040000000000000900200000134008000000100000000000000000031a00080000044300c0d00ee00003802013f008000000100000000000000000031910080000044380c2e00ee00003802700000806086010064090055310031002000008100400000000000001002000040361380208206820682448540004000400013812082268226824405010001004000138020820684068444854000400040001381208226842684440501000100403713802082068b068b448540004000400013812082268b268b440501000100400013802082068d068d448540004000400013812082268d268d440501000100403813802082068f068f448540004000400013812082268f268f440501000100400013802082069106914485400040004000138120822691269144050100010040391380208206930693448540004000400013812082269326934405010001004000138020820695069544854000400040001381208226952695440501000100403a1380208206970697448540004000400013812082269726974405010001004000138020820699069944854000400040001381208226992699440501000100403b13802082069b069b448540004000400013812082269b269b440501000100400013802082069d069d448540004000400013812082269d269d440501000100403c13802082069f069f448540004000400013812082269f269f44050100010040001380208206a106a144854000400040001381208226a126a1440501000100403d1380208206a306a344854000400040001381208226a326a344050100010040001380208206a506a544854000400040001381208226a526a5440501000100403000802082050d040d000500010001403100802082052e042e0005000100010100008000420100000000301400000059450380a03a0750045001010440041e0120008000000100000000000000000059440380a03a075804580101044004260121008000000100000000000000000059430380a03a0760046001010448041e012f008000000100000000000000000059420380a03a07680468010104480426700000806086010064090055210021002000008100400000000000006007000001310080000001000000000000000000612603802002050e057046000000000001300080000001000000000000000000613f0380200205160571460000000000612703802002050f05724600000000006100038020020517057346000000000061280380200205100574460000000000610003802002051805754600000000006129038020020511057646000000000061000380200205190577460000000000612a0380200205120578460000000000610003802002051a0579460000000000612b038020020513057a460000000000610003802002051b057b460000000000612c038020020514057c460000000000610003802002051c057d460000000000612d038020020515057e460000000000610003802002051d057f460000000000400000802086450c64090005e0ffe0ff6100038010400586000000001032547661000380104085860000000098badcfe400104801001058605865805440c0002700144805085058705865855f0fff0ff6c00048050850586058658050f000f00650103802002050e05864602050e4600650003802002050f05864602050f460065000380200205100586460205104600650003802002051105864602051146006500038020020512058646020512460065000380200205130586460205134600650003802002051405864602051446006500038020020515058646020515460065000380200205160587460205164600650003802002051705874602051746006500038020020518058746020518460065000380200205190587460205194600650003802002051a05874602051a4600650003802002051b05874602051b4600650003802002051c05874602051c4600650003802002051d05874602051d4600314f0080000000000c0a0830000000000100008000000100000000e00000000031f40080000000000c8908ee440e3800319f0080000000000c8a08ee441638003141008000000c0b0c003ee20000000031410080000000000c0a0830000000000100008000000100000000e000000000314400800000840e0c2f00ee00003c0070000080608601006409005531003100200000810040000000000000400000000135008000000100000000000000000031420080000044400c0d00ee0000380231430080000044480c2e00ee00003802403400802082052f042f00050002000220000081004000000000000030000000403200802082050d040d000500010001403300802082052e042e000500010001400000802086450c64090005e0ffe0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805440c00026c0104805085050b050b58050f000f006590038020020530040b0002053046006500038020020531240b0002053146006500038020020532440b0002053246006500038020020533640b0002053346006500038020020534840b0002053446006500038020020535a40b0002053546006500038020020536c40b0002053646006500038020020537e40b0002053746006500038020020538040b0002053846006500038020020539240b000205394600650003802002053a440b0002053a4600650003802002053b640b0002053b4600650003802002053c840b0002053c4600650003802002053da40b0002053d4600650003802002053ec40b0002053e4600650003802002053fe40b0002053f4600011f0080004201000000003030000000594a0380a03a0750045001010430040e594b0380a03a0758045801010430041601190080004201000000003008000000594c0380a03a0760046001010438040e01220080000001000000000000000000594d0380a03a0768046801010438041620000081004000000000000000020000010000800042010000000020003c0000314400800000840e0c2f00ee00003c00403400802082052f042f000500fe00fe400000802086450c64090005d0ffd0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805440c00026c0104805085050b050b58050f000f00013500800000010000000000000000006592038020020540040b0002054046006500038020020541240b0002054146006500038020020542440b0002054246006500038020020543640b0002054346006500038020020544840b0002054446006500038020020545a40b0002054546006500038020020546c40b0002054646006500038020020547e40b0002054746006523038020020548040b0002054846006500038020020549240b000205494600650003802002054a440b0002054a4600650003802002054b640b0002054b4600650003802002054c840b0002054c4600650003802002054da40b0002054d4600650003802002054ec40b0002054e4600650003802002054fe40b0002054f4600011f008000420100000000300004000059440380a03a0750045001010440040e594b0380a03a0758045801010440041601190080000001000000000000000000594c0380a03a0760046001010448040e594d0380a03a07680468010104480416200000800040000000000000300f0000610003801040050a0000000010325476610003801040850a0000000098badcfe6901038020810684840a1001020002006900038020810682040a100102000200400213802002068406844482048100006100038020022684248100000000000040011381208226842684440501000100400013802002068206824482048100006100038020022682248100000000000040011381208226822682440501000100400013802002068b06824482a4080000610203802002268b2682440000000000400113812082268b268b440501000100400013802002068d06844482a4080000610003802002268d2684440000000000400113812082268d268d440501000100690000802082450ca408000101000100680000802082650ca40800011f001f00400213802002068f06824482440c0000400203802002268f26824402640c0000400113812082268f268f440501000100400013802002069106844482440c0000400003802002269126844402640c0000400113812082269126914405010001004100008020824120a408000103000300490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069306824482440c0000400003802002269326824402640c000040011381208226932693440501000100400013802002069506844482440c0000400003802002269526844402640c000040011381208226952695440501000100690000802082450ca408000102000200680000802082650ca40800011e001e00400213802002069706824482440c0000400203802002269726824402640c000040011381208226972697440501000100400013802002069906844482440c0000400003802002269926844402640c0000400113812082269926994405010001004100008020824120a408000105000500490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069b06824482440c0000400003802002269b26824402640c0000400113812082269b269b440501000100400013802002069d06844482440c0000400003802002269d26844402640c0000400113812082269d269d4405010001004100008020824120a408000106000600490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069f06824482440c0000400003802002269f26824402640c0000400113812082269f269f44050100010040001380200206a106844482440c000040000380200226a126844402640c000040011381208226a126a14405010001004100008020824120a408000107000700490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c402000000000000040011380200206a306824482440c000040000380200226a326824402640c000040011381208226a326a344050100010040001380200206a506844482440c000040000380200226a526844402640c000040011381208226a526a5440501000100400000801486460c6409000501000100680100801081450c440c0001010001004001008014a1450c440c000110001000610000801041050a00000000ffffffff6801008010011130040a0001440c000070000080608601006409005501000100200000810040000000000000400a0000400000801486460c6409000501000100680100801081450c440c0001010001004001008014a1450c440c000110001000610000801041050a00000000ffffffff6801008010011130040a0001440c0000720002801009040a050a0109050a050a72000280a00a840a850a020a850a850a3146448100001470248200fb000000023147448100001472248b00fb000000023148448100001474248f00fb000000023149448100001476249300fb00000002314a448100001478249700fb00000002314b44810000147a249b00fb00000002314c44810000147c249f00fb00000002314d44810000147e24a300fb0000000231400080000044300c0d00ee0000380231410080000044380c2e00ee00003802700000806086010064090055110011002000008100400000000000001002000040361380208206820682448540004000400013812082268226824405010001004000138020820684068444854000400040001381208226842684440501000100403713802082068b068b448540004000400013812082268b268b440501000100400013802082068d068d448540004000400013812082268d268d440501000100403813802082068f068f448540004000400013812082268f268f440501000100400013802082069106914485400040004000138120822691269144050100010040391380208206930693448540004000400013812082269326934405010001004000138020820695069544854000400040001381208226952695440501000100403a1380208206970697448540004000400013812082269726974405010001004000138020820699069944854000400040001381208226992699440501000100403b13802082069b069b448540004000400013812082269b269b440501000100400013802082069d069d448540004000400013812082269d269d440501000100403c13802082069f069f448540004000400013812082269f269f44050100010040001380208206a106a144854000400040001381208226a126a1440501000100403d1380208206a306a344854000400040001381208226a326a344050100010040001380208206a506a544854000400040001381208226a526a5440501000100403000802082050d040d000500010001403100802082052e042e000500010001612603802002050e057046000000000061000380200205160571460000000000612703802002050f05724600000000006100038020020517057346000000000061280380200205100574460000000000610003802002051805754600000000006129038020020511057646000000000061000380200205190577460000000000612a0380200205120578460000000000610003802002051a0579460000000000612b038020020513057a460000000000610003802002051b057b460000000000612c038020020514057c460000000000610003802002051c057d460000000000612d038020020515057e460000000000610003802002051d057f460000000000610003801040050a0000000010325476611103801040850a0000000098badcfe400104801001050a050a580564090002700144805085050b050a5855f0fff0ff6c0004805085050a050a58050f000f00650103802002050e050a4602050e4600650003802002050f050a4602050f46006500038020020510050a4602051046006500038020020511050a4602051146006500038020020512050a4602051246006500038020020513050a4602051346006500038020020514050a4602051446006500038020020515050a4602051546006500038020020516050b4602051646006500038020020517050b4602051746006500038020020518050b4602051846006500038020020519050b460205194600650003802002051a050b4602051a4600650003802002051b050b4602051b4600650003802002051c050b4602051c4600650003802002051d050b4602051d4600610000802042450a0000000000002020319f0080000000000c0a0830000000000100008000000100000000e00000000031440080000000000c8908ee440e3800314f0080000000000c8a08ee44163800314e008000000c0b0c003ee200000000314e0080000000000c0a0830000000000100008000000100000000e000000000314400800000840e0c2f00ee00003c002000008100400000000000003000000031420080000044400c0d00ee0000380231430080000044480c2e00ee00003802403400802082052f042f00050002000220000081004000000000000030000000403200802082050d040d000500010001403300802082052e042e000500010001610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805640900026c0104805085050b050b58050f000f006590038020020530040b0002053046006500038020020531240b0002053146006500038020020532440b0002053246006500038020020533640b0002053346006500038020020534840b0002053446006500038020020535a40b0002053546006500038020020536c40b0002053646006500038020020537e40b0002053746006521038020020538040b0002053846006500038020020539240b000205394600650003802002053a440b0002053a4600650003802002053b640b0002053b4600650003802002053c840b0002053c4600650003802002053da40b0002053d4600650003802002053ec40b0002053e4600650003802002053fe40b0002053f4600011f0080004201000000003010000000594a0380a03a0750045001010430040e594b0380a03a0758045801010430041601190080000001000000000000000000594c0380a03a0760046001010438040e594d0380a03a07680468010104380416200000810040000000000000f0010000010000800042010000000020003c0000314400800000840e0c2f00ee00003c00403400802082052f042f000500fe00fe400000802086450c64090005f0fff0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805440c00026c0104805085050b050b58050f000f006592038020020540040b0002054046006500038020020541240b0002054146006500038020020542440b0002054246006500038020020543640b0002054346006500038020020544840b0002054446006500038020020545a40b0002054546006500038020020546c40b0002054646006500038020020547e40b0002054746006523038020020548040b0002054846006500038020020549240b000205494600650003802002054a440b0002054a4600650003802002054b640b0002054b4600650003802002054c840b0002054c4600650003802002054da40b0002054d4600650003802002054ec40b0002054e4600650003802002054fe40b0002054f4600011f008000420100000000300004000059440380a03a0750045001010440040e594b0380a03a0758045801010440041601190080000001000000000000000000594c0380a03a0760046001010448040e594d0380a03a076804680101044804162f000480004000000000000010000000200000800040000000000000d03100006900008060a545098409000103000300400100801406260944090002c40c0000410000802002650ca408000184090000690100802082650c640c000103000300690000802081450c84090001080008006900008020814509a40900010a000a004e0300802202010044080002640c0000400000802002450844080002640c0000610000802002050a00200000000000004001008020026508040a00026408000040000080200245094409000224080000400100802002650944090002440c00004000008054a1150d24090005000000006c0100805085150d140d00050f000f004000008054a1250d24090005010001006c0100805085250d240d00050f000f004000008054a1350d24090005020002006c0100805085350d340d00050f000f004000008054a1450d24090005030003006c0100805085450d440d00050f000f004000008054a1550d24090005040004006c0100805085550d540d00050f000f004000008054a1650d24090005050005006c0100805085650d640d00050f000f004000008054a1750d24090005060006006c0100805085750d740d00050f000f004000008054a1850d24090005070007006c0100805085850d840d00050f000f00610001802002058245082200000000004e0100802202058304820002a4080000610000802002050a00200000000000004001008020022583040a000224820000690000802082450ca408000101000100680000802082650ca40800011f001f004e0200802202058404820002440c0000610000802002050a0020000000000000400300802002258424820002640c00004001008020022584040a0002248400004100008020824120a408000103000300490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c40200000000000004e0100802202058504820002440c0000610000802002050a0020000000000000400000802002258524820002640c00004001008020022585040a000224850000690000802082450ca408000102000200680000802082650ca40800011e001e004e0200802202058604820002440c0000610000802002050a0020000000000000400300802002258624820002640c00004001008020022586040a0002248600004100008020824120a408000105000500490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c40200000000000004e0100802202058704820002440c0000610000802002050a0020000000000000400000802002258724820002640c00004001008020022587040a0002248700004100008020824120a408000106000600490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c40200000000000004e0100802202058804820002440c0000610000802002050a0020000000000000400000802002258824820002640c00004001008020022588040a0002248800004100008020824120a408000107000700490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c40200000000000004e0100802202058904820002440c0000610000802002050a0020000000000000400000802002258924820002640c00004001008020022589040a000224890000610000802002058a6409000000000000400100802082058b048a000500020002610000802006052e0408000000000000400100802082052f042e0005000800086100008020020580440900000000000040010080208205810480000500020002610000806006250904090000000000006100048020410550000000000000000061000480204105520000000000000000610004802041055400000000000000006100048020410556000000000000000061000480204105580000000000000000610004802041055a0000000000000000610004802041055c0000000000000000610004802041055e00000000000000006100048020410560000000000000000061000480204105620000000000000000610004802041056400000000000000006100048020410566000000000000000061000480204105680000000000000000610004802041056a0000000000000000610004802041056c0000000000000000610004802041056e00000000000000000109008000000100000000300000000070000080608601002409005520002000200000810040000000000000601b0000400000806086250924090005c1ffc1ff720002801009040a050a0109050a050a72000280a00a840a850a020a850a850a6100008010011130140d0000000000006100008010010131240d0000000000006100008010011131340d00000000000031464081000014700c8200ff000030026100008010011130440d00000000000031478081000014720c8300ff000030026100008010010131540d0000000000003148c081000014740c8400ff000030026100008010011131640d00000000000031494081000014760c8500ff000030026100008010011130740d000000000000314a8081000014780c8600ff000030026100008010010131840d000000000000314bc0810000147a0c8700ff00003002314c40810000147c0c8800ff00003002314d80810000147e0c8900ff0000300231400080000044300c2e00ee0000380231410080000044380c2f00ee000038024036108020820582048200854000400040001081208225822482000501000100403710802082058304830085400040004000108120822583248300050100010040381080208205840484008540004000400010812082258424840005010001004039108020820585048500854000400040001081208225852485000501000100403a108020820586048600854000400040001081208225862486000501000100403b108020820587048700854000400040001081208225872487000501000100403c108020820588048800854000400040001081208225882488000501000100403d108020820589048900854000400040001081208225892489000501000100403000802082052e042e000500010001403100802082052f042f000500010001612603802002050e057046000000000061000380200205160571460000000000612703802002050f05724600000000006100038020020517057346000000000061280380200205100574460000000000610003802002051805754600000000006129038020020511057646000000000061000380200205190577460000000000612a0380200205120578460000000000610003802002051a0579460000000000612b038020020513057a460000000000610003802002051b057b460000000000612c038020020514057c460000000000610003802002051c057d460000000000612d038020020515057e460000000000610003802002051d057f460000000000610000802042450a0000000000002020319f0080000000000c0a0830000000000100008000000100000000e00000000031440080000000000c8a08ee440e3800314f0080000000000c8b08ee44163800314e008000000c0b0c003ee200000000314e0080000000000c0a0830000000000100008000000100000000e000000000700000806086010024090065000000000109008000000100000000000000000020000081004000000000000040050000720002801009040b050b0109050b050b72000280a00a840b850b020a850b850b6100008010011131140d0000000000006100008010011130240d0000000000006100008010010131340d0000000000003146c081000014700c8200ff000030026100008010011131440d00000000000031474081000014720c8300ff000030026100008010011130540d00000000000031488081000014740c8400ff000030026100008010010131640d0000000000003149c081000014760c8500ff000030026100008010011131740d000000000000314a4081000014780c8600ff000030026100008010011130840d000000000000314b80810000147a0c8700ff00003002314cc0810000147c0c8800ff00003002012d0080000001000000000000000000314d40810000147e0c8900ff00003002314400800000840e0c8000ee00003c000100008000420100000000200c000000314500800000841e0c8100ee00003c0031420080000044400c2e00ee0000380231430080000044480c2f00ee00003802400000806086250924090035e0ffe0ff403200802082052e042e000500010001403300802082052f042f0005000100010120008000000100000000000000000059440380a03a0750045001010430040e59400380a03a0758045801010430041659410380a03a0760046001010438040e594f0380a03a076804680101043804160134008000000100000000000000000031a00080000044300c2e00ee00003802013f008000000100000000000000000031910080000044380c2f00ee000038024036108020820582048200854000400040001081208225822482000501000100403710802082058304830085400040004000108120822583248300050100010040381080208205840484008540004000400010812082258424840005010001004039108020820585048500854000400040001081208225852485000501000100403a108020820586048600854000400040001081208225862486000501000100403b108020820587048700854000400040001081208225872487000501000100403c108020820588048800854000400040001081208225882488000501000100403d108020820589048900854000400040001081208225892489000501000100403000802082052e042e000500010001403100802082052f042f0005000100010100008000420100000000303000000059420380a03a0750045001010440041e59450380a03a0758045801010440042659430380a03a0760046001010448041e594f0380a03a07680468010104480426612603802002050e057046000000000061000380200205160571460000000000612703802002050f05724600000000006100038020020517057346000000000061280380200205100574460000000000610003802002051805754600000000006129038020020511057646000000000061000380200205190577460000000000612a0380200205120578460000000000610003802002051a0579460000000000612b038020020513057a460000000000610003802002051b057b460000000000612c038020020514057c460000000000610003802002051c057d460000000000612d038020020515057e460000000000610003802002051d057f460000000000314d0080000000000c0a0830000000000100008000000100000000e00000000031a40080000000000c8a08ee440e3800319f0080000000000c8b08ee44163800314e008000000c0b0c003ee200000000314e0080000000000c0a0830000000000100008000000100000000e000000000200000810040000000000000e0faffff4000008060862509240900053f003f00610003801040050b0000000010325476610003801040850b0000000098badcfe6901038020810685840b1001020002006900038020810683040b100102000200400213802002068506854482048200006100038020022685248200000000000040011381208226852685440501000100400013802002068306834482048200006100038020022683248200000000000040011381208226832683440501000100400013802002068c06834482a4080000610203802002268c2683440000000000400113812082268c268c440501000100400013802002068e06854482a4080000610003802002268e2685440000000000400113812082268e268e440501000100690000802082450ca408000101000100680000802082650ca40800011f001f00400213802002069006834482440c0000400203802002269026834402640c000040011381208226902690440501000100400013802002069206854482440c0000400003802002269226854402640c0000400113812082269226924405010001004100008020824120a408000103000300490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069406834482440c0000400003802002269426834402640c000040011381208226942694440501000100400013802002069606854482440c0000400003802002269626854402640c000040011381208226962696440501000100690000802082450ca408000102000200680000802082650ca40800011e001e00400213802002069806834482440c0000400203802002269826834402640c000040011381208226982698440501000100400013802002069a06854482440c0000400003802002269a26854402640c0000400113812082269a269a4405010001004100008020824120a408000105000500490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069c06834482440c0000400003802002269c26834402640c0000400113812082269c269c440501000100400013802002069e06854482440c0000400003802002269e26854402640c0000400113812082269e269e4405010001004100008020824120a408000106000600490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c402000000000000040011380200206a006834482440c000040000380200226a026834402640c000040011381208226a026a044050100010040001380200206a206854482440c000040000380200226a226854402640c000040011381208226a226a24405010001004100008020824120a408000107000700490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c402000000000000040011380200206a406834482440c000040000380200226a426834402640c000040011381208226a426a444050100010040001380200206a606854482440c000040000380200226a626854402640c000040011381208226a626a6440501000100400000801486460c24090005e1ffe1ff680100801081450c440c0001010001004001008014a1450c440c000110001000610000801041050b00000000ffffffff680100801001950d040b0001440c000070000080608601002409005521002100200000810040000000000000c0010000400000801486460c24090005e1ffe1ff680100801081450c440c0001010001004001008014a1450c440c000110001000610000801041050b00000000ffffffff680100801001950d040b0001440c0000720002801009040b050b0109050b050b72000280a00a840b850b020a850b850b6102008010010131940d0000000000002200848100c0000020010000200100006100000010011131140d0000000000006100000010011130240d0000000000006100000010010131340d0000000000003146c40100001470248300fb000000026100000010011131440d0000000000003147440100001472248c00fb000000026100000010011130540d0000000000003148840100001474249000fb000000026100000010010131640d0000000000003149c40100001476249400fb000000026100000010011131740d000000000000314a440100001478249800fb000000026100000010011130840d000000000000314b84010000147a249c00fb00000002314cc4010000147c24a000fb00000002012d0080000001000000000000000000314d44010000147e24a400fb0000000225000480004000000000000010000000314400800000840e0c8000ee00003c000100008000420100000000200c000000314500800000841e0c8100ee00003c0031420080000044400c2e00ee0000380231430080000044480c2f00ee00003802403200802082052e042e000500010001403300802082052f042f0005000100010120008000000100000000000000000059440380a03a0750045001010430040e59400380a03a0758045801010430041659410380a03a0760046001010438040e594f0380a03a07680468010104380416200000810040000000000000900200000134008000000100000000000000000031a00080000044300c2e00ee00003802013f008000000100000000000000000031910080000044380c2f00ee00003802700000806086010024090055310031002000008100400000000000001002000040361380208206830683448540004000400013812082268326834405010001004000138020820685068544854000400040001381208226852685440501000100403713802082068c068c448540004000400013812082268c268c440501000100400013802082068e068e448540004000400013812082268e268e4405010001004038138020820690069044854000400040001381208226902690440501000100400013802082069206924485400040004000138120822692269244050100010040391380208206940694448540004000400013812082269426944405010001004000138020820696069644854000400040001381208226962696440501000100403a138020820698069844854000400040001381208226982698440501000100400013802082069a069a448540004000400013812082269a269a440501000100403b13802082069c069c448540004000400013812082269c269c440501000100400013802082069e069e448540004000400013812082269e269e440501000100403c1380208206a006a044854000400040001381208226a026a044050100010040001380208206a206a244854000400040001381208226a226a2440501000100403d1380208206a406a444854000400040001381208226a426a444050100010040001380208206a606a644854000400040001381208226a626a6440501000100403000802082052e042e000500010001403100802082052f042f0005000100010100008000420100000000301400000059450380a03a0750045001010440041e0120008000000100000000000000000059440380a03a075804580101044004260121008000000100000000000000000059430380a03a0760046001010448041e012f008000000100000000000000000059420380a03a07680468010104480426700000806086010024090055210021002000008100400000000000006007000001310080000001000000000000000000612603802002050e057046000000000001300080000001000000000000000000613f0380200205160571460000000000612703802002050f05724600000000006100038020020517057346000000000061280380200205100574460000000000610003802002051805754600000000006129038020020511057646000000000061000380200205190577460000000000612a0380200205120578460000000000610003802002051a0579460000000000612b038020020513057a460000000000610003802002051b057b460000000000612c038020020514057c460000000000610003802002051c057d460000000000612d038020020515057e460000000000610003802002051d057f460000000000400000802086450c24090005e0ffe0ff6100038010400587000000001032547661000380104085870000000098badcfe400104801001058705875805440c0002400104801081058805875805100010006c01058050850587058758050f000f00650103802002050e05874602050e4600650003802002050f05874602050f460065000380200205100587460205104600650003802002051105874602051146006500038020020512058746020512460065000380200205130587460205134600650003802002051405874602051446006500038020020515058746020515460065000380200205160588460205164600650003802002051705884602051746006500038020020518058846020518460065000380200205190588460205194600650003802002051a05884602051a4600650003802002051b05884602051b4600650003802002051c05884602051c4600650003802002051d05884602051d4600314f0080000000000c0a0830000000000100008000000100000000e00000000031f40080000000000c8a08ee440e3800319f0080000000000c8b08ee441638003141008000000c0b0c003ee20000000031410080000000000c0a0830000000000100008000000100000000e000000000314400800000840e0c8000ee00003c0070000080608601002409005531003100200000810040000000000000400000000135008000000100000000000000000031420080000044400c2e00ee0000380231430080000044480c2f00ee000038024034008020820580048000050002000220000081004000000000000030000000403200802082052e042e000500010001403300802082052f042f000500010001400000802086450c24090005e0ffe0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805440c00026c0104805085050b050b58050f000f006590038020020530040b0002053046006500038020020531240b0002053146006500038020020532440b0002053246006500038020020533640b0002053346006500038020020534840b0002053446006500038020020535a40b0002053546006500038020020536c40b0002053646006500038020020537e40b0002053746006500038020020538040b0002053846006500038020020539240b000205394600650003802002053a440b0002053a4600650003802002053b640b0002053b4600650003802002053c840b0002053c4600650003802002053da40b0002053d4600650003802002053ec40b0002053e4600650003802002053fe40b0002053f4600011f0080004201000000003030000000594a0380a03a0750045001010430040e594b0380a03a0758045801010430041601190080004201000000003008000000594c0380a03a0760046001010438040e01220080000001000000000000000000594d0380a03a0768046801010438041620000081004000000000000000020000010000800042010000000020003c0000314400800000840e0c8000ee00003c0040340080208205800480000500fe00fe400000802086450c24090005d0ffd0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805440c00026c0104805085050b050b58050f000f00013500800000010000000000000000006592038020020540040b0002054046006500038020020541240b0002054146006500038020020542440b0002054246006500038020020543640b0002054346006500038020020544840b0002054446006500038020020545a40b0002054546006500038020020546c40b0002054646006500038020020547e40b0002054746006523038020020548040b0002054846006500038020020549240b000205494600650003802002054a440b0002054a4600650003802002054b640b0002054b4600650003802002054c840b0002054c4600650003802002054da40b0002054d4600650003802002054ec40b0002054e4600650003802002054fe40b0002054f4600011f008000420100000000300004000059440380a03a0750045001010440040e594b0380a03a0758045801010440041601190080000001000000000000000000594c0380a03a0760046001010448040e594d0380a03a07680468010104480416200000800040000000000000e00f0000610003801040050a0000000010325476610003801040850a0000000098badcfe6901038020810685840a1001020002006900038020810683040a100102000200400213802002068506854482048200006100038020022685248200000000000040011381208226852685440501000100400013802002068306834482048200006100038020022683248200000000000040011381208226832683440501000100400013802002068c06834482a4080000610203802002268c2683440000000000400113812082268c268c440501000100400013802002068e06854482a4080000610003802002268e2685440000000000400113812082268e268e440501000100690000802082450ca408000101000100680000802082650ca40800011f001f00400213802002069006834482440c0000400203802002269026834402640c000040011381208226902690440501000100400013802002069206854482440c0000400003802002269226854402640c0000400113812082269226924405010001004100008020824120a408000103000300490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069406834482440c0000400003802002269426834402640c000040011381208226942694440501000100400013802002069606854482440c0000400003802002269626854402640c000040011381208226962696440501000100690000802082450ca408000102000200680000802082650ca40800011e001e00400213802002069806834482440c0000400203802002269826834402640c000040011381208226982698440501000100400013802002069a06854482440c0000400003802002269a26854402640c0000400113812082269a269a4405010001004100008020824120a408000105000500490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c4020000000000000400113802002069c06834482440c0000400003802002269c26834402640c0000400113812082269c269c440501000100400013802002069e06854482440c0000400003802002269e26854402640c0000400113812082269e269e4405010001004100008020824120a408000106000600490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c402000000000000040011380200206a006834482440c000040000380200226a026834402640c000040011381208226a026a044050100010040001380200206a206854482440c000040000380200226a226854402640c000040011381208226a226a24405010001004100008020824120a408000107000700490000802282450ca408000500000000610100802002650c440c000000000000610000802002450c402000000000000040011380200206a406834482440c000040000380200226a426834402640c000040011381208226a426a444050100010040001380200206a606854482440c000040000380200226a626854402640c000040011381208226a626a6440501000100400000801486460c2409000501000100680100801081450c440c0001010001004001008014a1450c440c000110001000610000801041050a00000000ffffffff680100801001950d040a0001440c000070000080608601002409005501000100200000810040000000000000f00a0000400000801486460c2409000501000100680100801081450c440c0001010001004001008014a1450c440c000110001000610000801041050a00000000ffffffff680100801001950d040a0001440c0000720002801009040a050a0109050a050a72000280a00a840a850a020a850a850a6102008010011130940d0000000000002200448100c0000010010000100100006100000010010131140d0000000000006100000010011131240d0000000000006100000010011130340d0000000000003146840100001470248300fb000000026100000010010131440d0000000000003147c40100001472248c00fb000000026100000010011131540d0000000000003148440100001474249000fb000000026100000010011130640d0000000000003149840100001476249400fb000000026100000010010131740d000000000000314ac40100001478249800fb000000026100000010011131840d000000000000314b44010000147a249c00fb00000002314c84010000147c24a000fb00000002314dc4010000147e24a400fb000000022500048000400000000000001000000031400080000044300c2e00ee0000380231410080000044380c2f00ee00003802700000806086010024090055110011002000008100400000000000001002000040361380208206830683448540004000400013812082268326834405010001004000138020820685068544854000400040001381208226852685440501000100403713802082068c068c448540004000400013812082268c268c440501000100400013802082068e068e448540004000400013812082268e268e4405010001004038138020820690069044854000400040001381208226902690440501000100400013802082069206924485400040004000138120822692269244050100010040391380208206940694448540004000400013812082269426944405010001004000138020820696069644854000400040001381208226962696440501000100403a138020820698069844854000400040001381208226982698440501000100400013802082069a069a448540004000400013812082269a269a440501000100403b13802082069c069c448540004000400013812082269c269c440501000100400013802082069e069e448540004000400013812082269e269e440501000100403c1380208206a006a044854000400040001381208226a026a044050100010040001380208206a206a244854000400040001381208226a226a2440501000100403d1380208206a406a444854000400040001381208226a426a444050100010040001380208206a606a644854000400040001381208226a626a6440501000100403000802082052e042e000500010001403100802082052f042f000500010001612603802002050e057046000000000061000380200205160571460000000000612703802002050f05724600000000006100038020020517057346000000000061280380200205100574460000000000610003802002051805754600000000006129038020020511057646000000000061000380200205190577460000000000612a0380200205120578460000000000610003802002051a0579460000000000612b038020020513057a460000000000610003802002051b057b460000000000612c038020020514057c460000000000610003802002051c057d460000000000612d038020020515057e460000000000610003802002051d057f460000000000610003801040050a0000000010325476611103801040850a0000000098badcfe400104801001050a050a580524090002400104801081050b050a5805100010006c0105805085050a050a58050f000f00650103802002050e050a4602050e4600650003802002050f050a4602050f46006500038020020510050a4602051046006500038020020511050a4602051146006500038020020512050a4602051246006500038020020513050a4602051346006500038020020514050a4602051446006500038020020515050a4602051546006500038020020516050b4602051646006500038020020517050b4602051746006500038020020518050b4602051846006500038020020519050b460205194600650003802002051a050b4602051a4600650003802002051b050b4602051b4600650003802002051c050b4602051c4600650003802002051d050b4602051d4600610000802042450a0000000000002020319f0080000000000c0a0830000000000100008000000100000000e00000000031440080000000000c8a08ee440e3800314f0080000000000c8b08ee44163800314e008000000c0b0c003ee200000000314e0080000000000c0a0830000000000100008000000100000000e000000000314400800000840e0c8000ee00003c002000008100400000000000003000000031420080000044400c2e00ee0000380231430080000044480c2f00ee000038024034008020820580048000050002000220000081004000000000000030000000403200802082052e042e000500010001403300802082052f042f000500010001610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805240900026c0104805085050b050b58050f000f006590038020020530040b0002053046006500038020020531240b0002053146006500038020020532440b0002053246006500038020020533640b0002053346006500038020020534840b0002053446006500038020020535a40b0002053546006500038020020536c40b0002053646006500038020020537e40b0002053746006521038020020538040b0002053846006500038020020539240b000205394600650003802002053a440b0002053a4600650003802002053b640b0002053b4600650003802002053c840b0002053c4600650003802002053da40b0002053d4600650003802002053ec40b0002053e4600650003802002053fe40b0002053f4600011f0080004201000000003010000000594a0380a03a0750045001010430040e594b0380a03a0758045801010430041601190080000001000000000000000000594c0380a03a0760046001010438040e594d0380a03a07680468010104380416200000810040000000000000f0010000010000800042010000000020003c0000314400800000840e0c8000ee00003c0040340080208205800480000500fe00fe400000802086450c24090005f0fff0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805440c00026c0104805085050b050b58050f000f006592038020020540040b0002054046006500038020020541240b0002054146006500038020020542440b0002054246006500038020020543640b0002054346006500038020020544840b0002054446006500038020020545a40b0002054546006500038020020546c40b0002054646006500038020020547e40b0002054746006523038020020548040b0002054846006500038020020549240b000205494600650003802002054a440b0002054a4600650003802002054b640b0002054b4600650003802002054c840b0002054c4600650003802002054da40b0002054d4600650003802002054ec40b0002054e4600650003802002054fe40b0002054f4600011f008000420100000000300004000059440380a03a0750045001010440040e594b0380a03a0758045801010440041601190080000001000000000000000000594c0380a03a0760046001010448040e594d0380a03a076804680101044804167000048020820100a40c0061000000007000848020820100c40c0061000000002e00048200c0000010000000100000002f0004800040000000000000100000006100008020020130e4090000000000006100008020020131040c0000000000006100008020024170c40900000000000001090080000001000000003000000000
    return c;
}
#ifndef MICRO_DECL_float_tile_16x8_blocked_8x8
#define MICRO_DECL_float_tile_16x8_blocked_8x8
typedef struct {
    float8 x[2];
} float_tile_16x8_blocked_8x8;
DECLARE_2D_TILE_OPS(float_tile_16x8_blocked_8x8,float,8,8,8,2,1)
#endif
#define ugemm_vs_sg_tile_m 16
#define ugemm_vs_sg_tile_n 8
#define ugemm_vs_wg_tile_m 128
#define ugemm_vs_wg_tile_n 32
#define ugemm_vs_sg_per_wg_m 8
#define ugemm_vs_sg_per_wg_n 4
#define ugemm_vs_sg_per_wg_k 1
#define ugemm_vs_slm_size 8192
#define ugemm_vs_barrier_count  1
#define ugemm_vs_systolic  1
typedef float_tile_16x8_blocked_8x8 ugemm_vs_c_type;
#define ugemm_vs_c_type_block0 8
#define ugemm_vs_c_type_nblock0 2
#define ugemm_vs_c_type_block1 8
#define ugemm_vs_c_type_nblock1 1
ugemm_vs_c_type ugemm_vs(const global half* a, int lda, const local half* b, int ldb, int m, int n, int k, int i0, int j0, int h0, int local_id_m, int local_id_n, const local char* slm) {
    ugemm_vs_c_type c;
    __asm__ volatile("{\n"
            ".implicit_PSEUDO_INPUT %0 offset=1280 size=256\n"
            ".implicit_PSEUDO_INPUT %1 offset=1536 size=256\n"
            ".decl COPY2 v_type=G type=uq num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY2 offset=264 size=8\n"
            ".decl COPY3 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY3 offset=272 size=4\n"
            ".decl COPY4 v_type=G type=ud num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY4 offset=256 size=4\n"
            ".decl COPY5 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY5 offset=276 size=4\n"
            ".decl COPY6 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY6 offset=280 size=4\n"
            ".decl COPY7 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY7 offset=284 size=4\n"
            ".decl COPY8 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY8 offset=288 size=4\n"
            ".decl COPY9 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY9 offset=292 size=4\n"
            ".decl COPY10 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY10 offset=296 size=4\n"
            ".decl COPY11 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY11 offset=300 size=4\n"
            ".decl COPY12 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY12 offset=304 size=4\n"
            ".decl COPY13 v_type=G type=d num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY13 offset=308 size=4\n"
            ".decl COPY14 v_type=G type=ud num_elts=1\n"
            ".implicit_PSEUDO_INPUT COPY14 offset=260 size=4\n"
            "fence_sw\n"
            "mov (M1_NM, 1) COPY2(0,0)<1> %2(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY3(0,0)<1> %3(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY4(0,0)<1> %4(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY5(0,0)<1> %5(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY6(0,0)<1> %6(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY7(0,0)<1> %7(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY8(0,0)<1> %8(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY9(0,0)<1> %9(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY10(0,0)<1> %10(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY11(0,0)<1> %11(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY12(0,0)<1> %12(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY13(0,0)<1> %13(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) COPY14(0,0)<1> %14(0,0)<1;1,0>\n"
            ".decl CLOBBER0 v_type=G type=ud num_elts=2\n"
            ".implicit_PSEUDO_INPUT CLOBBER0 offset=312 size=8\n"
            ".decl CLOBBER1 v_type=G type=ud num_elts=240\n"
            ".implicit_PSEUDO_INPUT CLOBBER1 offset=320 size=960\n"
            ".decl CLOBBER2 v_type=G type=ud num_elts=208\n"
            ".implicit_PSEUDO_INPUT CLOBBER2 offset=1792 size=832\n"
            "fence_sw\n"
            "add (M1,1) CLOBBER0(0,0)<1> CLOBBER0(0,0)<0;1,0> 0xcafefadf:ud\n"
            "fence_sw\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY2(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY3(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY4(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY5(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY6(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY7(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY8(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY9(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY10(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY11(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY12(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY13(0,0)<1;1,0>\n"
            "mov (M1_NM, 1) V0(0,0)<1> COPY14(0,0)<1;1,0>\n"
            "mov (M1,2) CLOBBER0(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(6,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(8,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(10,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(12,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(14,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(16,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(18,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(20,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(22,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(24,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(26,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER1(28,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(6,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(8,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(10,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(12,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(14,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(16,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(18,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(20,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(22,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) CLOBBER2(24,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %0(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %0(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %0(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %0(6,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %1(0,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %1(2,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %1(4,0)<1> 0xAAAAAAAA:ud\n"
            "mov (M1,16) %1(6,0)<1> 0xAAAAAAAA:ud\n"
            ".decl DUMMY_DPAS_SRC v_type=G type=ud num_elts=64 alias=<CLOBBER1,0>\n"
            ".decl DUMMY_DPAS_DST v_type=G type=f num_elts=64 alias=<DUMMY_DPAS_SRC,0>\n"
            "dpas.bf.bf.8.1 (M1,8) DUMMY_DPAS_DST.0 V0.0 DUMMY_DPAS_SRC.0 DUMMY_DPAS_SRC(0,0)\n"
            "barrier\n"
            "fence_sw\n"
            "add (M1,1) CLOBBER0(0,0)<1> CLOBBER0(0,0)<0;1,0> 0xfadecaff:ud\n"
            "fence_sw\n"
        "}\n"
    : "=rw"(c.x[0]), "=rw"(c.x[1]) : "rw.u"(a), "rw.u"(lda), "rw.u"(b), "rw.u"(ldb), "rw.u"(m), "rw.u"(n), "rw.u"(k), "rw.u"(i0), "rw.u"(j0), "rw.u"(h0), "rw.u"(local_id_m), "rw.u"(local_id_n), "rw.u"(slm));
    // @_u_@1 610900802002c5094070000000000000610900802041417000000000ffffffff01090080000001000000003000000000610000802002e5090030000000000000610000802002050c00310000000000007a00008010c5240c100005018409000069000080208285088408000101000100690000802082a508a408000101000100690000802081450c84090001040004004001008060022509440c0006240900005b0000806086440944090501a4090800690200806086c50c24090001010001006c0100806086050ac40c00011f001f004e00008022020100c40c0002440800004000008020024508c40c000244080000610000802002050b00200000000000004004008020066508040a0002640800004001008020026508040b000264080000410000806006412044090001a4080000490000802206450c44090002a4080000610100806002650c440c000000000000610000802006450c4020000000000000400100806006050804080006440c000041000080600681206409000184080000490000802206850c6409000284080000610100806002a50c840c000000000000610000802006850c80200000000000004e0100802202010044080002840c0000400000802002450844080002840c0000610000802002050a0020000000000000400000802002650864080002a40c00004001008020026508040a0002640800005b0000806086040804080506640910004000008020a1450c240c000110001000400000802426650c24090006c4080000400000802426850c44090006e4080000400200802422450c440c0002640c0000620000802082650c640c005110001000620300802082850c840c005108000800690000802082a50c8408000105000500690000802082050e84080001050005007005008060820100440c005110001000200000810040000000000000c01c00006800008020816509a409000501000100410000802002450c84080001a4090000690100802082450c440c00010300030069000080208265096409000108000800690000802081c50ca4090001070007004001008020026509c40c0002640900006900008020812509840900010a000a004e0000802202010044080002440c0000400000802002450844080002440c0000610000802002050a00200000000000004001008020026508040a00026408000040000080200225092409000224080000400100802002450924090002640900006500008010810100840900250100010040000081608605080408000580008000610001802002054045082200000000004e010080220205410440000284080000610000802002050a00200000000000004001008020022541040a000224400000690000802082c50c8408000101000100680000802082e50c840800011f001f004e0200802202054204400002c40c0000610000802002050a0020000000000000400300802002254224400002e40c00004001008020022542040a000224420000410000802082c1208408000103000300490000802282c50c8408000500000000610100802002e50cc40c000000000000610000802002c50cc0200000000000004e0100802202054304400002c40c0000610000802002050a0020000000000000400000802002254324400002e40c00004001008020022543040a000224430000690000802082c50c8408000102000200680000802082e50c840800011e001e004e0200802202054404400002c40c0000610000802002050a0020000000000000400300802002254424400002e40c00004001008020022544040a000224440000410000802082c1208408000105000500490000802282c50c8408000500000000610100802002e50cc40c000000000000610000802002c50cc0200000000000004e0100802202054504400002c40c0000610000802002050a0020000000000000400000802002254524400002e40c00004001008020022545040a000224450000410000802082c1208408000106000600490000802282c50c8408000500000000610100802002e50cc40c000000000000610000802002c50cc0200000000000004e0100802202054604400002c40c0000610000802002050a0020000000000000400000802002254624400002e40c00004001008020022546040a000224460000410000802082c1208408000107000700490000802282c50c8408000500000000610100802002e50cc40c000000000000610000802002c50cc0200000000000004e0100802202054704400002c40c0000610000802002050a0020000000000000400000802002254724400002e40c00004001008020022547040a0002244700006100008020020548440900000000000040010080208205490448000500010001610000802002050d2409000000000000610000802006051704080000000000006100008060066509040900000000000061000480204105280000000000000000610004802041052a0000000000000000610004802041052c0000000000000000610004802041052e000000000000000061000480204105300000000000000000610004802041053200000000000000006100048020410534000000000000000061000480204105360000000000000000010900800000010000000030000000007000008060860100640900552000200020000081004000000000000080100000400000806086650964090005c1ffc1ff3143008000000c380c4000ff00002c023144008000000c390c4100ff00002c023145008000000c3a0c4200ff00002c023146008000000c3b0c4300ff00002c023147008000000c3c0c4400ff00002c023148008000000c3d0c4500ff00002c023149008000000c3e0c4600ff00002c02314a008000000c3f0c4700ff00002c02314100800000240f0c1700ee00003402403310802002054004400082040e000040001081208225402440000501000100403410802002054104410082040e000040001081208225412441000501000100403510802002054204420082a40c000040001081208225422442000501000100403610802002054304430082a40c000040001081208225432443000501000100403710802002054404440082040e000040001081208225442444000501000100403810802002054504450082040e000040001081208225452445000501000100403910802002054604460082a40c000040001081208225462446000501000100403a10802002054704470082a40c0000400010812082254724470005010001004031008020820517041700050001000161230380100106180538460000000000610003801001061c853846000000000061240380100116180539460000000000610003801001161c85394600000000006125038010010619053a460000000000610003801001061d853a4600000000006126038010011619053b460000000000610003801001161d853b460000000000612703801001061a053c460000000000610003801001061e853c460000000000612803801001161a053d460000000000610003801001161e853d460000000000612903801001061b053e460000000000610003801001061f853e460000000000612a03801001161b053f460000000000610003801001161f853f460000000000610000802042450a0000000000002020319d0080000000000c0a0830000000000100008000000100000000e00000000031400080000000000c4808ee24183400314c0080000000000c4908ee241c3400314d008000000c0b0c003ee200000000314d0080000000000c0a0830000000000100008000000100000000e0000000007000008060860100640900650000000001090080000001000000000000000000200000810040000000000000300400003143008000000c380c4000ff00002c023144008000000c390c4100ff00002c023145008000000c3a0c4200ff00002c023146008000000c3b0c4300ff00002c023147008000000c3c0c4400ff00002c023148008000000c3d0c4500ff00002c023149008000000c3e0c4600ff00002c02314a008000000c3f0c4700ff00002c02013f008000000100000000000000000031400080000084180c0d00ee00003c0031420080000024130c1700ee00003402400000806086650964090035e0ffe0ff403000802082050d040d00050002000240320080208205170417000500010001012100800000010000000000000000005a400380a03a0728042801010418040f5a4f0380a03a0730043001010420040f013f008000000100000000000000000031a00080000084180c0d00ee00003c00319100800000240f0c1700ee00003402403310802002054004400082040e000040001081208225402440000501000100403410802002054104410082040e000040001081208225412441000501000100403510802002054204420082a40c000040001081208225422442000501000100403610802002054304430082a40c000040001081208225432443000501000100403710802002054404440082040e000040001081208225442444000501000100403810802002054504450082040e000040001081208225452445000501000100403910802002054604460082a40c000040001081208225462446000501000100403a10802002054704470082a40c000040001081208225472447000501000100403000802082050d040d000500fe00fe40310080208205170417000500010001012200800000010000000000000000005a400380a03a072804280101041804135a4f0380a03a073004300101042004130130008000000100000000000000000061230380100106180538460000000000613c03801001061c853846000000000061240380100116180539460000000000610003801001161c85394600000000006125038010010619053a460000000000610003801001061d853a4600000000006126038010011619053b460000000000610003801001161d853b460000000000612703801001061a053c460000000000610003801001061e853c460000000000612803801001161a053d460000000000610003801001161e853d460000000000612903801001061b053e460000000000610003801001061f853e460000000000612a03801001161b053f460000000000610003801001161f853f460000000000314d0080000000000c0a0830000000000100008000000100000000e00000000031a00080000000000c4808ee24183400319c0080000000000c4908ee241c3400314d008000000c0b0c003ee200000000314d0080000000000c0a0830000000000100008000000100000000e000000000200000810040000000000000f0fbffff4000008060866509640900053f003f0070010080608601006409005521002100200000810040000000000000700300005b0000801485440c64098505a4090800400101801481450c450c2205e0ffe0ff7001008050850100440c003500000000400000805085450c440c0005ffffffff61000480204105380000000000000000610004802041053a0000000000000000610004802041053c0000000000000000610004802041053e0000000000000000200000910040000000000000e00200007005008050850100440c003500000000400000805085450c440c0005ffffffff31e3008000000c380c4000ff00002c02403310802002054004400082040e000040001081208225402440000501000100200000910040000000000000800200007003008050850100440c003500000000400000805085450c440c0005ffffffff3144008000000c390c4100ff00002c02403410802002054104410082040e000040001081208225412441000501000100200000910040000000000000200200007003008050850100440c003500000000400000805085450c440c0005ffffffff3145008000000c3a0c4200ff00002c02403510802002054204420082a40c000040001081208225422442000501000100200000910040000000000000c00100007003008050850100440c003500000000400000805085450c440c0005ffffffff3146008000000c3b0c4300ff00002c02403610802002054304430082a40c000040001081208225432443000501000100200000910040000000000000600100007003008050850100440c003500000000400000805085450c440c0005ffffffff3147008000000c3c0c4400ff00002c02403710802002054404440082040e000040001081208225442444000501000100200000910040000000000000000100007003008050850100440c003500000000400000805085450c440c0005ffffffff3148008000000c3d0c4500ff00002c02403810802002054504450082040e000040001081208225452445000501000100200000910040000000000000a00000007003008050850100440c003500000000400000805085450c440c0005e7ffe7ff3149008000000c3e0c4600ff00002c02403910802002054604460082a40c00004000108120822546244600050100010020000091004000000000000040000000314a008000000c3f0c4700ff00002c02403a10802002054704470082a40c000040001081208225472447000501000100013f008000000100000000000000000031400080000084180c0d00ee00003c0031420080000024130c1700ee00003402403000802082050d040d00050002000240320080208205170417000500010001012100800000010000000000000000005a400380a03a0728042801010418040f5a4f0380a03a0730043001010420040f013f008000000100000000000000000031a00080000084180c0d00ee00003c00700000806086010064090055210021002000008100400000000000002000000031a100800000240f0c1700ee00003402403000802082050d040d000500fe00fe200000810040000000000000300000000131008000000100000000000000000040030080208205170417000500010001010000800042010000000030050000005a4e0380a03a072804280101041804135a4f0380a03a0730043001010420041320000081004000000000000010040000013e008000000100000000000000000061f30380100106180538460000000000613c03801001061c853846000000000061240380100116180539460000000000610003801001161c85394600000000006125038010010619053a460000000000610003801001061d853a4600000000006126038010011619053b460000000000610003801001161d853b460000000000612703801001061a053c460000000000610003801001061e853c460000000000612803801001161a053d460000000000610003801001161e853d460000000000612903801001061b053e460000000000610003801001061f853e460000000000612a03801001161b053f460000000000610003801001161f853f460000000000314d0080000000000c0a0830000000000100008000000100000000e00000000031a00080000000000c4808ee24183400319c0080000000000c4908ee241c3400314d008000000c0b0c003ee200000000314d0080000000000c0a0830000000000100008000000100000000e000000000013f008000000100000000000000000031400080000084180c0d00ee00003c00700000806086010064090055310031002000008100400000000000003000000001000080004201000000002000c0000031420080000024130c1700ee00003402403000802082050d040d000500020002200000810040000000000000300000000131008000000100000000000000000040320080208205170417000500010001400000802086a50c64090005e0ffe0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805a40c00026c0104805085050b050b58050f000f00659103802002050f050b4602050f46006500038020020510050b4602051046006500038020020511050b4602051146006500038020020512050b460205124600011900800042010000000030014000005a450380a03a0728042801010418040f012f00800000010000000000000000005a460380a03a0730043001010420040f200000810040000000000000100100000100008000420100000000206000000031400080000084180c0d00ee00003c00403000802082050d040d000500fe00fe400000802086a50c64090005d0ffd0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805a40c00026c0104805085050b050b58050f000f0001000080004201000000002000c000006592038020020513050b4602051346006500038020020514050b4602051446006500038020020515050b4602051546006500038020020516050b460205164600011900800042010000000030200000005a400380a03a072804280101041804135a460380a03a073004300101042004132000008000400000000000004007000070000080608601006409005501000100200000810040000000000000200700005b0000801485440c64098505a40908007001008050850100440c003500000000400000805085450c440c0005ffffffff61000480204105380000000000000000610004802041053a0000000000000000610004802041053c0000000000000000610004802041053e0000000000000000200000910040000000000000e00200007005008050850100440c003500000000400000805085450c440c0005ffffffff31e3008000000c380c4000ff00002c02403310802002054004400082040e000040001081208225402440000501000100200000910040000000000000800200007003008050850100440c003500000000400000805085450c440c0005ffffffff3144008000000c390c4100ff00002c02403410802002054104410082040e000040001081208225412441000501000100200000910040000000000000200200007003008050850100440c003500000000400000805085450c440c0005ffffffff3145008000000c3a0c4200ff00002c02403510802002054204420082a40c000040001081208225422442000501000100200000910040000000000000c00100007003008050850100440c003500000000400000805085450c440c0005ffffffff3146008000000c3b0c4300ff00002c02403610802002054304430082a40c000040001081208225432443000501000100200000910040000000000000600100007003008050850100440c003500000000400000805085450c440c0005ffffffff3147008000000c3c0c4400ff00002c02403710802002054404440082040e000040001081208225442444000501000100200000910040000000000000000100007003008050850100440c003500000000400000805085450c440c0005ffffffff3148008000000c3d0c4500ff00002c02403810802002054504450082040e000040001081208225452445000501000100200000910040000000000000a00000007003008050850100440c003500000000400000805085450c440c0005e7ffe7ff3149008000000c3e0c4600ff00002c02403910802002054604460082a40c00004000108120822546244600050100010020000091004000000000000040000000314a008000000c3f0c4700ff00002c02403a10802002054704470082a40c000040001081208225472447000501000100314100800000240f0c1700ee000034024031008020820517041700050001000161d30380100106180538460000000000610003801001061c853846000000000061240380100116180539460000000000610003801001161c853946000000000061f5038010010619053a460000000000610003801001061d853a4600000000006126038010011619053b460000000000610003801001161d853b460000000000612703801001061a053c460000000000610003801001061e853c460000000000612803801001161a053d460000000000610003801001161e853d460000000000612903801001061b053e460000000000610003801001061f853e460000000000612a03801001161b053f460000000000610003801001161f853f460000000000610000802042450a0000000000002020319d0080000000000c0a0830000000000100008000000100000000e00000000031400080000000000c4808ee24183400314c0080000000000c4908ee241c3400314d008000000c0b0c003ee200000000314d0080000000000c0a0830000000000100008000000100000000e00000000031400080000084180c0d00ee00003c00700000806086010064090055110011002000008100400000000000002000000031420080000024130c1700ee00003402403000802082050d040d0005000200022000008100400000000000002000000040320080208205170417000500010001610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805640900026c0104805085050b050b58050f000f00659103802002050f050b4602050f46006500038020020510050b4602051046006500038020020511050b4602051146006500038020020512050b460205124600011900800042010000000030010000005a450380a03a0728042801010418040f5a460380a03a0730043001010420040f200000810040000000000000000100000100008000420100000000206000000031400080000084180c0d00ee00003c00403000802082050d040d000500fe00fe400000802086a50c64090005f0fff0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805a40c00026c0104805085050b050b58050f000f006592038020020513050b4602051346006500038020020514050b4602051446006500038020020515050b4602051546006500038020020516050b460205164600011900800042010000000030200000005a400380a03a072804280101041804135a460380a03a073004300101042004132f000480004000000000000010000000200000800040000000000000302000006800008020816509a409000501000100410000802002450c84080001a4090000690100802082450c440c00010300030069000080208265096409000108000800690000802081c50ca4090001070007004001008020026509c40c0002640900006900008020812509840900010a000a004e0000802202010044080002440c0000400000802002450844080002440c0000610000802002050a00200000000000004001008020026508040a000264080000400000802002250924090002240800004001008020024509240900026409000065000080108101008409002501000100400000816086050804080005800080004000008014826609640c000501000100680100801081650964090001010001004001008014a165096409000108000800610000801041250e00000000ff00ff006801008010011130240e00016409000061000380104005280000000010325476690103802081064004281001020002004e0103802202050b0640440244080000610003802002050a01204600000000006102038020020640050b4600000000004002038020022640050a4602640800004e0203802202050b0640440284080000610003802002050a01204600000000006102038020020642050b4600000000004002038020022642050a460226404400690000802082c50c8408000101000100680000802082e50c840800011f001f004e0203802202050b06404402c40c0000610003802002050a01204600000000006102038020020644050b460000000000400003802002264426404402e40c00004001038020022644050a460226444400410000802082c1208408000103000300490000802282c50c8408000500000000610100802002e50cc40c000000000000610000802002c50cc0200000000000004e0103802202050b06404402c40c0000610003802002050a01204600000000006102038020020646050b460000000000400003802002264626404402e40c00004001038020022646050a460226464400690000802082c50c8408000102000200680000802082e50c840800011e001e004e0203802202050b06404402c40c0000610003802002050a01204600000000006102038020020648050b460000000000400003802002264826404402e40c00004001038020022648050a460226484400410000802082c1208408000105000500490000802282c50c8408000500000000610100802002e50cc40c000000000000610000802002c50cc0200000000000004e0103802202050b06404402c40c0000610003802002050a0120460000000000610203802002064a050b460000000000400003802002264a26404402e40c0000400103802002264a050a4602264a4400410000802082c1208408000106000600490000802282c50c8408000500000000610100802002e50cc40c000000000000610000802002c50cc0200000000000004e0103802202050b06404402c40c0000610003802002050a0120460000000000610203802002064c050b460000000000400003802002264c26404402e40c0000400103802002264c050a4602264c4400410000802082c1208408000107000700490000802282c50c8408000500000000610100802002e50cc40c000000000000610000802002c50cc0200000000000004e0103802202050b06404402c40c0000610003802002050a0120460000000000610203802002064e050b460000000000400003802002264e26404402e40c0000400103802002264e050a4602264e44006100008020020550440900000000000040010080208205510450000500010001610000802002050d2409000000000000610000802006051704080000000000006100008060066509040900000000000061000480204105280000000000000000610004802041052a0000000000000000610004802041052c0000000000000000610004802041052e0000000000000000610004802041053000000000000000006100048020410532000000000000000061000480204105340000000000000000610004802041053600000000000000000109008000000100000000300000000070000080608601006409005520002000200000810040000000000000c0110000400000806086650964090005c1ffc1ff720002801009040a050a0109050a050a72000280a00a840a850a020a850a850a3143438100000c38144000fb000000023144438100000c39144200fb000000023145438100000c3a144400fb000000023146438100000c3b144600fb000000023147438100000c3c144800fb000000023148438100000c3d144a00fb000000023149438100000c3e144c00fb00000002314a438100000c3f144e00fb00000002314100800000240f0c1700ee00003402403313802002064006404482040e000040001381208226402640440501000100403413802002064206424482a40c000040001381208226422642440501000100403513802002064406444482040e000040001381208226442644440501000100403613802002064606464482a40c000040001381208226462646440501000100403713802002064806484482040e000040001381208226482648440501000100403813802002064a064a4482a40c0000400013812082264a264a440501000100403913802002064c064c4482040e0000400013812082264c264c440501000100403a13802002064e064e4482a40c0000400013812082264e264e4405010001004031008020820517041700050001000161230380100106180538460000000000610003801001061c853846000000000061240380100116180539460000000000610003801001161c85394600000000006125038010010619053a460000000000610003801001061d853a4600000000006126038010011619053b460000000000610003801001161d853b460000000000612703801001061a053c460000000000610003801001061e853c460000000000612803801001161a053d460000000000610003801001161e853d460000000000612903801001061b053e460000000000610003801001061f853e460000000000612a03801001161b053f460000000000610003801001161f853f460000000000610000802042450a0000000000002020319d0080000000000c0a0830000000000100008000000100000000e00000000031400080000000000c5008ee24183400314c0080000000000c5108ee241c3400314d008000000c0b0c003ee200000000314d0080000000000c0a0830000000000100008000000100000000e000000000700000806086010064090065000000000109008000000100000000000000000020000081004000000000000050040000720002801009040b050b0109050b050b72000280a00a840b850b020a850b850b3143438100000c38144000fb000000023144438100000c39144200fb000000023145438100000c3a144400fb000000023146438100000c3b144600fb000000023147438100000c3c144800fb000000023148438100000c3d144a00fb000000023149438100000c3e144c00fb00000002314a438100000c3f144e00fb00000002013f008000000100000000000000000031400080000084180c0d00ee00003c0031420080000024130c1700ee00003402400000806086650964090035e0ffe0ff403000802082050d040d00050002000240320080208205170417000500010001012100800000010000000000000000005a400380a03a0728042801010418040f5a4f0380a03a0730043001010420040f013f008000000100000000000000000031a00080000084180c0d00ee00003c00319100800000240f0c1700ee00003402403313802002064006404482040e000040001381208226402640440501000100403413802002064206424482a40c000040001381208226422642440501000100403513802002064406444482040e000040001381208226442644440501000100403613802002064606464482a40c000040001381208226462646440501000100403713802002064806484482040e000040001381208226482648440501000100403813802002064a064a4482a40c0000400013812082264a264a440501000100403913802002064c064c4482040e0000400013812082264c264c440501000100403a13802002064e064e4482a40c0000400013812082264e264e440501000100403000802082050d040d000500fe00fe40310080208205170417000500010001012200800000010000000000000000005a400380a03a072804280101041804135a4f0380a03a073004300101042004130130008000000100000000000000000061230380100106180538460000000000613c03801001061c853846000000000061240380100116180539460000000000610003801001161c85394600000000006125038010010619053a460000000000610003801001061d853a4600000000006126038010011619053b460000000000610003801001161d853b460000000000612703801001061a053c460000000000610003801001061e853c460000000000612803801001161a053d460000000000610003801001161e853d460000000000612903801001061b053e460000000000610003801001061f853e460000000000612a03801001161b053f460000000000610003801001161f853f460000000000314d0080000000000c0a0830000000000100008000000100000000e00000000031a00080000000000c5008ee24183400319c0080000000000c5108ee241c3400314d008000000c0b0c003ee200000000314d0080000000000c0a0830000000000100008000000100000000e000000000200000810040000000000000d0fbffff4000008060866509640900053f003f0070010080608601006409005521002100200000810040000000000000700400005b0000801485440c64098505a4090800400101801481450c450c2205e0ffe0ff7001008050850100440c003500000000400000805085450c440c0005ffffffff61000480204105380000000000000000610004802041053a0000000000000000610004802041053c0000000000000000610004802041053e0000000000000000200000910040000000000000e00300007005008050850100440c003500000000400000805085450c440c0005ffffffff720002801009040b050b0109050b050b72000280a00a840b850b020a850b850b31f3438100000c38144000fb00000002403313802002064006404482040e000040001381208226402640440501000100200000910040000000000000600300007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040b050b0109050b050b72010280a00a840b850b020a850b850b3144438100000c39144200fb00000002403413802002064206424482a40c000040001381208226422642440501000100200000910040000000000000e00200007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040b050b0109050b050b72010280a00a840b850b020a850b850b3145438100000c3a144400fb00000002403513802002064406444482040e000040001381208226442644440501000100200000910040000000000000600200007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040b050b0109050b050b72010280a00a840b850b020a850b850b3146438100000c3b144600fb00000002403613802002064606464482a40c000040001381208226462646440501000100200000910040000000000000e00100007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040b050b0109050b050b72010280a00a840b850b020a850b850b3147438100000c3c144800fb00000002403713802002064806484482040e000040001381208226482648440501000100200000910040000000000000600100007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040b050b0109050b050b72010280a00a840b850b020a850b850b3148438100000c3d144a00fb00000002403813802002064a064a4482a40c0000400013812082264a264a440501000100200000910040000000000000e00000007004008050850100440c003500000000400000805085450c440c0005e7ffe7ff720502801009040b050b0109050b050b72010280a00a840b850b020a850b850b3149438100000c3e144c00fb00000002403913802002064c064c4482040e0000400013812082264c264c44050100010020000091004000000000000060000000720302801009040b050b0109050b050b72010280a00a840b850b020a850b850b314a438100000c3f144e00fb00000002403a13802002064e064e4482a40c0000400013812082264e264e440501000100013f008000000100000000000000000031400080000084180c0d00ee00003c0031420080000024130c1700ee00003402403000802082050d040d00050002000240320080208205170417000500010001012100800000010000000000000000005a400380a03a0728042801010418040f5a4f0380a03a0730043001010420040f013f008000000100000000000000000031a00080000084180c0d00ee00003c00700000806086010064090055210021002000008100400000000000002000000031a100800000240f0c1700ee00003402403000802082050d040d000500fe00fe200000810040000000000000300000000131008000000100000000000000000040030080208205170417000500010001010000800042010000000030050000005a4e0380a03a072804280101041804135a4f0380a03a0730043001010420041320000081004000000000000010040000013e008000000100000000000000000061f30380100106180538460000000000613c03801001061c853846000000000061240380100116180539460000000000610003801001161c85394600000000006125038010010619053a460000000000610003801001061d853a4600000000006126038010011619053b460000000000610003801001161d853b460000000000612703801001061a053c460000000000610003801001061e853c460000000000612803801001161a053d460000000000610003801001161e853d460000000000612903801001061b053e460000000000610003801001061f853e460000000000612a03801001161b053f460000000000610003801001161f853f460000000000314d0080000000000c0a0830000000000100008000000100000000e00000000031a00080000000000c5008ee24183400319c0080000000000c5108ee241c3400314d008000000c0b0c003ee200000000314d0080000000000c0a0830000000000100008000000100000000e000000000013f008000000100000000000000000031400080000084180c0d00ee00003c00700000806086010064090055310031002000008100400000000000003000000001000080004201000000002000c0000031420080000024130c1700ee00003402403000802082050d040d000500020002200000810040000000000000300000000131008000000100000000000000000040320080208205170417000500010001400000802086a50c64090005e0ffe0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805a40c00026c0104805085050b050b58050f000f00659103802002050f050b4602050f46006500038020020510050b4602051046006500038020020511050b4602051146006500038020020512050b460205124600011900800042010000000030014000005a450380a03a0728042801010418040f012f00800000010000000000000000005a460380a03a0730043001010420040f200000810040000000000000100100000100008000420100000000206000000031400080000084180c0d00ee00003c00403000802082050d040d000500fe00fe400000802086a50c64090005d0ffd0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805a40c00026c0104805085050b050b58050f000f0001000080004201000000002000c000006592038020020513050b4602051346006500038020020514050b4602051446006500038020020515050b4602051546006500038020020516050b460205164600011900800042010000000030200000005a400380a03a072804280101041804135a460380a03a073004300101042004132000008000400000000000004008000070000080608601006409005501000100200000810040000000000000200800005b0000801485440c64098505a40908007001008050850100440c003500000000400000805085450c440c0005ffffffff61000480204105380000000000000000610004802041053a0000000000000000610004802041053c0000000000000000610004802041053e0000000000000000200000910040000000000000e00300007005008050850100440c003500000000400000805085450c440c0005ffffffff720002801009040a050a0109050a050a72000280a00a840a850a020a850a850a31f3438100000c38144000fb00000002403313802002064006404482040e000040001381208226402640440501000100200000910040000000000000600300007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040a050a0109050a050a72010280a00a840a850a020a850a850a3144438100000c39144200fb00000002403413802002064206424482a40c000040001381208226422642440501000100200000910040000000000000e00200007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040a050a0109050a050a72010280a00a840a850a020a850a850a3145438100000c3a144400fb00000002403513802002064406444482040e000040001381208226442644440501000100200000910040000000000000600200007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040a050a0109050a050a72010280a00a840a850a020a850a850a3146438100000c3b144600fb00000002403613802002064606464482a40c000040001381208226462646440501000100200000910040000000000000e00100007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040a050a0109050a050a72010280a00a840a850a020a850a850a3147438100000c3c144800fb00000002403713802002064806484482040e000040001381208226482648440501000100200000910040000000000000600100007004008050850100440c003500000000400000805085450c440c0005ffffffff720502801009040a050a0109050a050a72010280a00a840a850a020a850a850a3148438100000c3d144a00fb00000002403813802002064a064a4482a40c0000400013812082264a264a440501000100200000910040000000000000e00000007004008050850100440c003500000000400000805085450c440c0005e7ffe7ff720502801009040a050a0109050a050a72010280a00a840a850a020a850a850a3149438100000c3e144c00fb00000002403913802002064c064c4482040e0000400013812082264c264c44050100010020000091004000000000000060000000720302801009040a050a0109050a050a72010280a00a840a850a020a850a850a314a438100000c3f144e00fb00000002403a13802002064e064e4482a40c0000400013812082264e264e440501000100314100800000240f0c1700ee000034024031008020820517041700050001000161d30380100106180538460000000000610003801001061c853846000000000061240380100116180539460000000000610003801001161c853946000000000061f5038010010619053a460000000000610003801001061d853a4600000000006126038010011619053b460000000000610003801001161d853b460000000000612703801001061a053c460000000000610003801001061e853c460000000000612803801001161a053d460000000000610003801001161e853d460000000000612903801001061b053e460000000000610003801001061f853e460000000000612a03801001161b053f460000000000610003801001161f853f460000000000610000802042450a0000000000002020319d0080000000000c0a0830000000000100008000000100000000e00000000031400080000000000c5008ee24183400314c0080000000000c5108ee241c3400314d008000000c0b0c003ee200000000314d0080000000000c0a0830000000000100008000000100000000e00000000031400080000084180c0d00ee00003c00700000806086010064090055110011002000008100400000000000002000000031420080000024130c1700ee00003402403000802082050d040d0005000200022000008100400000000000002000000040320080208205170417000500010001610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805640900026c0104805085050b050b58050f000f00659103802002050f050b4602050f46006500038020020510050b4602051046006500038020020511050b4602051146006500038020020512050b460205124600011900800042010000000030010000005a450380a03a0728042801010418040f5a460380a03a0730043001010420040f200000810040000000000000000100000100008000420100000000206000000031400080000084180c0d00ee00003c00403000802082050d040d000500fe00fe400000802086a50c64090005f0fff0ff610003801040050b0000000010325476610003801040850b0000000098badcfe400104801001050b050b5805a40c00026c0104805085050b050b58050f000f006592038020020513050b4602051346006500038020020514050b4602051446006500038020020515050b4602051546006500038020020516050b460205164600011900800042010000000030200000005a400380a03a072804280101041804135a460380a03a073004300101042004137000048020820100640c0061000000007000848020820100840c0061000000002e00048200c0000010000000100000002f0004800040000000000000100000006100008020020130e4090000000000006100008020020131040c0000000000006100008020024170c40900000000000001090080000001000000003000000000
    return c;
}

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define DIV_UP(x, y) (((x) + (y)-1) / (y))
#define sg_per_wg (ugemm_kq_sg_per_wg_m * ugemm_kq_sg_per_wg_n)
#define q_tile_sg_n DIV_UP(ugemm_kq_wg_tile_n, sg_per_wg)
typedef ugemm_kq_c_type s_tile_type;
typedef ugemm_vs_c_type a_tile_type;
DECLARE_2D_TILE(q_tile_type, uint, SUBGROUP_SIZE, D_MAX / 2, 1, 1, q_tile_sg_n)
#ifdef BLOCK_Q
DECLARE_2D_TILE_BLOCK_OPS(
 q_tile_type, uint, SUBGROUP_SIZE, D_MAX / 2, 1, 1, q_tile_sg_n)
#elif Q_ALIGN < 4
DECLARE_2D_TILE_LOAD_PACKED_HALF(
 q_tile_type, SUBGROUP_SIZE, D_MAX / 2, 1, 1, q_tile_sg_n)
#endif
#ifdef BLOCK_A
DECLARE_2D_TILE(a_tile_type_half, half, SUBGROUP_SIZE, ugemm_vs_sg_tile_m, 1, 1,
 ugemm_vs_sg_tile_n)
#else
DECLARE_2D_TILE(a_tile_type_half, half, SUBGROUP_SIZE, ugemm_vs_sg_tile_m, 8, 1,
 ugemm_vs_sg_tile_n / 8)
#endif
DECLARE_2D_TILE(s_tile_type_half2, uint, SUBGROUP_SIZE, ugemm_kq_c_type_block0,
 ugemm_kq_c_type_block1 / 2, ugemm_kq_c_type_nblock0,
 ugemm_kq_c_type_nblock1)
DECLARE_2D_TILE(
 s_sum_tile_type, float, SUBGROUP_SIZE, ugemm_kq_sg_tile_n, 1, 1, 1)
DECLARE_2D_TILE(
 a_scale_tile_type, float, SUBGROUP_SIZE, ugemm_vs_sg_tile_n, 1, 1, 1)
DECLARE_2D_TILE(mask_tile_type, half, SUBGROUP_SIZE, ugemm_kq_c_type_block0, ugemm_kq_c_type_block1, ugemm_kq_c_type_nblock0, ugemm_kq_c_type_nblock1)
DECLARE_2D_TILE(mask_tile_type_float, float, SUBGROUP_SIZE, ugemm_kq_c_type_block0, ugemm_kq_c_type_block1, ugemm_kq_c_type_nblock0, ugemm_kq_c_type_nblock1)
#ifdef BLOCK_A
DECLARE_2D_TILE_BLOCK_OPS(a_tile_type_half, half, SUBGROUP_SIZE,
 ugemm_vs_sg_tile_m, 1, 1, ugemm_vs_sg_tile_n)
#endif
#ifdef BLOCK_2D_A
DECLARE_2D_TILE_BLOCK2D_OPS(a_tile_type_half, half, SUBGROUP_SIZE,
 ugemm_vs_sg_tile_m, 8, 1, ugemm_vs_sg_tile_n / 8)
#endif
#ifdef BLOCK_A
DECLARE_2D_TILE_COPY_REBLOCK(a_tile_type, SUBGROUP_SIZE, ugemm_vs_c_type_block0,
 ugemm_vs_c_type_block1, ugemm_vs_c_type_nblock0,
 ugemm_vs_c_type_nblock1, a_tile_type_half, SUBGROUP_SIZE,
 ugemm_vs_sg_tile_m, 1, 1, ugemm_vs_sg_tile_n)
#else
DECLARE_2D_TILE_COPY_REBLOCK(a_tile_type, SUBGROUP_SIZE, ugemm_vs_c_type_block0,
 ugemm_vs_c_type_block1, ugemm_vs_c_type_nblock0,
 ugemm_vs_c_type_nblock1, a_tile_type_half, SUBGROUP_SIZE,
 ugemm_vs_sg_tile_m, 8, 1, ugemm_vs_sg_tile_n / 8)
#endif
DECLARE_2D_TILE_VREDUCE(s_tile_type, SUBGROUP_SIZE, ugemm_kq_c_type_block0,
 ugemm_kq_c_type_block1, ugemm_kq_c_type_nblock0,
 ugemm_kq_c_type_nblock1, s_sum_tile_type, SUBGROUP_SIZE,
 ugemm_kq_sg_tile_n, 1, 1, 1)
DECLARE_2D_TILE_HREDUCE(s_tile_type, SUBGROUP_SIZE, ugemm_kq_c_type_block0,
 ugemm_kq_c_type_block1, ugemm_kq_c_type_nblock0,
 ugemm_kq_c_type_nblock1, mask_tile_type_float, SUBGROUP_SIZE,
 ugemm_kq_sg_tile_m, 1, 1, 1)
DECLARE_2D_TILE_HREDUCE(a_tile_type, SUBGROUP_SIZE, ugemm_vs_c_type_block0,
 ugemm_vs_c_type_block1, ugemm_vs_c_type_nblock0,
 ugemm_vs_c_type_nblock1, a_scale_tile_type, SUBGROUP_SIZE,
 ugemm_vs_sg_tile_n, 1, 1, 1)
#if ugemm_kq_wg_tile_n == ugemm_vs_wg_tile_n && (ugemm_kq_sg_tile_n % ugemm_vs_sg_tile_n) == 0
DECLARE_2D_TILE_RSELECT(a_scale_tile_type, SUBGROUP_SIZE, ugemm_vs_sg_tile_n, 1,
 1, 1, s_sum_tile_type, SUBGROUP_SIZE, ugemm_kq_sg_tile_n, 1, 1, 1)
#endif
#if PREFETCH_REMAINDER
#define cooperative_prefetch_2d_maybe_rem cooperative_prefetch_2d_rem
#else
#define cooperative_prefetch_2d_maybe_rem( ptr, r, c, rmax, cmax, ld, sg_id, n_sg, sg_size, caching) cooperative_prefetch_2d(ptr, rmax, cmax, ld, sg_id, n_sg, sg_size, caching)
#endif
#if TRANSPOSE_K
#define cooperative_prefetch_2d_k( ptr, r, c, rmax, cmax, ld, sg_id, n_sg, sg_size, caching) cooperative_prefetch_2d_maybe_rem( ptr, c, r, cmax, rmax, ld, sg_id, n_sg, sg_size, caching)
#else
#define cooperative_prefetch_2d_k cooperative_prefetch_2d_maybe_rem
#endif
#if REMAINDER_Q
#define tile_load_block_rem_q tile_load_block
#define tile_store_block_rem_q tile_store_block
#else
#define tile_load_block_rem_q(t, ptr, n, ld, off_r, off_c) tile_load_block(t, ptr, ld, off_r, off_c)
#define tile_store_block_rem_q(t, ptr, n, ld, off_r, off_c) tile_store_block(t, ptr, ld, off_r, off_c)
#endif
#define binary_add(x, y) ((x) + (y))
__attribute__((intel_reqd_sub_group_size(SUBGROUP_SIZE)))
KERNEL(micro_sdpa)(OPTIONAL_SHAPE_INFO_ARG
 const global half *K, const global half *Q, const global half *V,
 global half *A,
#if WITH_ATTN_MASK
 const global half *msk,
#endif
#if WITH_SCALE
 global SCALE_DATA_T *scale_ptr,
#endif
 int d, int k, int q) {
 uint sg_ij = sub_group_broadcast(get_local_id(1), 0);
 uint b0 = get_group_id(1);
 uint b1 = get_group_id(2);
 uint wg_j0 = get_group_id(0) * ugemm_kq_wg_tile_n;
 uint ldk = TRANSPOSE_K ? KEY_S3 : KEY_S2;
 uint ldq = QRY_S2;
 uint ldv = VAL_S2;
 uint lda = DST_S2;
 uint sg_i_kq = sg_ij % ugemm_kq_sg_per_wg_m;
 uint sg_j_kq = sg_ij / ugemm_kq_sg_per_wg_m;
 uint sg_i_vs = sg_ij % ugemm_vs_sg_per_wg_m;
 uint sg_j_vs = sg_ij / ugemm_vs_sg_per_wg_m;
#define Q_slm_size (D_MAX * ugemm_kq_wg_tile_n * sizeof(half))
#define S_slm_size (ugemm_kq_wg_tile_m * ugemm_kq_wg_tile_n * sizeof(half))
#define S_sum_slm_size (ugemm_kq_wg_tile_n * ugemm_kq_sg_per_wg_m * sizeof(float))
#define S_max_slm_size (ugemm_kq_wg_tile_n * sizeof(float))
#define ugemm_slm_size MAX(ugemm_kq_slm_size, ugemm_vs_slm_size)
 local char slm[Q_slm_size + S_slm_size + S_sum_slm_size + S_max_slm_size
 + ugemm_slm_size];
 local half *Q_slm = (local half *)&slm[0];
 local half *S_slm = (local half *)&slm[Q_slm_size];
 local float *S_sum_slm = (local float *)&slm[Q_slm_size + S_slm_size];
 local float *S_max_slm
 = (local float *)&slm[Q_slm_size + S_slm_size + S_sum_slm_size];
 local uint *ugemm_slm = (local uint *)&slm[Q_slm_size + S_slm_size
 + S_sum_slm_size + S_max_slm_size];
 const bool need_sum_barrier = (ugemm_vs_barrier_count == 0);
 K += KEY_OFF(b1, (b0 / KV_GROUP_SIZE), 0, 0) + INPUT1_OFFSET;
 Q += QRY_OFF(b1, b0, 0, 0) + INPUT0_OFFSET;
 V += VAL_OFF(b1, (b0 / KV_GROUP_SIZE), 0, 0) + INPUT2_OFFSET;
 A += DST_OFF(b1, b0, 0, 0, 0);
 __builtin_assume_aligned(K, K_ALIGN);
 __builtin_assume_aligned(Q, Q_ALIGN);
 __builtin_assume_aligned(V, V_ALIGN);
 __builtin_assume_aligned(A, A_ALIGN);
 q_tile_type Q_tile;
 uint q0_copy = q_tile_sg_n * sg_ij;
#ifdef BLOCK_Q
 tile_load_block_rem_q(
 &Q_tile, (global uint *)Q, q, ldq >> 1, 0, wg_j0 + q0_copy);
#elif Q_ALIGN >= 4
 tile_load(&Q_tile, (global uint *)Q, (d + 1) >> 1, q, ldq >> 1, 0,
 wg_j0 + q0_copy);
#else
 tile_load_packed_half(&Q_tile, Q, d, q, ldq, 0, wg_j0 + q0_copy);
#endif
#if WITH_SCALE
 #if INVERT_SCALE
 float iscale = convert_float(*scale_ptr);
 float scale = native_recip(iscale);
 #else
 float scale = convert_float(*scale_ptr);
 float iscale = native_recip(scale);
 #endif
#else
 float iscale = sqrt(convert_float(INPUT1_SIZE_X));
 float scale = native_recip(iscale);
#endif
 scale *= 1.442695f;
#ifdef PREFETCH_K0
 cooperative_prefetch_2d_k(K, d, k, ugemm_kq_wg_tile_m, PREFETCH_D_MAX, ldk,
 sg_ij, sg_per_wg, SUBGROUP_SIZE, LSC_LDCC_L1C_L3C);
#endif
 const uint n_col_sg = DIV_UP(ugemm_kq_wg_tile_n, SUBGROUP_SIZE * sg_per_wg);
 const float neg_inf = -INFINITY;
#pragma unroll
 for (int q = 0; q < n_col_sg; q++)
 intel_sub_group_block_write(
 (local uint *)&S_max_slm[(q + sg_ij * n_col_sg)
 * SUBGROUP_SIZE],
 as_uint(neg_inf));
 a_tile_type A_tile;
 tile_fill(A_tile, 0.0f);
 tile_store_t_sys_src1(
 Q_tile, (local uint *)&Q_slm[0], D_MAX / 2, q0_copy, 0);
 s_sum_tile_type S_sum_tile;
 s_sum_tile_type S_max_tile, S_max_tile_old;
 tile_fill(S_sum_tile, 0.0f);
 tile_fill(S_max_tile, -INFINITY);
 barrier(CLK_LOCAL_MEM_FENCE);
 for (int k0 = 0; k0 < k; k0 += ugemm_kq_wg_tile_m) {
 bool first = (k0 == 0);
 bool last = (k0 + ugemm_kq_wg_tile_m >= k);
 uint sg_i0_kq = sg_i_kq * ugemm_kq_sg_tile_m;
 uint sg_j0_kq = sg_j_kq * ugemm_kq_sg_tile_n;
#if WITH_ATTN_MASK
 mask_tile_type mask_tile;
 tile_load_t(&mask_tile, msk, q, k, q, sg_j0_kq + wg_j0, k0 + sg_i0_kq);
#endif
#if REMAINDER_K
 mask_tile_type_float k_mask;
#pragma unroll
 for (int ii = 0; ii < ugemm_kq_sg_tile_m / SUBGROUP_SIZE; ii++)
 k_mask.x[0][ii] = (k0 + sg_i0_kq + ii * SUBGROUP_SIZE
 + get_sub_group_local_id()
 < k)
 ? nan(0u)
 : -INFINITY;
#endif
 s_tile_type S_tile
 = ugemm_kq(K, ldk, Q_slm, D_MAX, k, ugemm_kq_wg_tile_n, d, k0,
 0, 0, sg_i_kq, sg_j_kq, (local char *)ugemm_slm);
#if WITH_ATTN_MASK
#define unscale(x) ((x)*iscale)
 mask_tile_type_float mask_tile_float;
 tile_copy(mask_tile, mask_tile_float);
 tile_elementwise(mask_tile_float, unscale);
 tile_binary(S_tile, mask_tile_float, binary_add);
#endif
#if REMAINDER_K
 tile_hbroadcast_min(&S_tile, k_mask);
#endif
 tile_vreduce_max(S_tile, &S_max_tile);
 tile_atomic_max_full(
 S_max_tile, S_max_slm, ugemm_kq_wg_tile_n, sg_j0_kq, 0);
 intel_work_group_barrier_arrive(CLK_LOCAL_MEM_FENCE);
#ifdef PREFETCH_V
 cooperative_prefetch_2d_maybe_rem(V, d, k - k0, D_MAX,
 (ugemm_kq_wg_tile_m * PREFETCH_D_MAX) / D_MAX, ldv, sg_ij,
 sg_per_wg, SUBGROUP_SIZE, LSC_LDCC_L1C_L3C);
#endif
#ifndef ALT_MAX
 intel_work_group_barrier_wait(CLK_LOCAL_MEM_FENCE);
 tile_load_full(&S_max_tile, S_max_slm, ugemm_kq_wg_tile_n, sg_j0_kq, 0);
#endif
 tile_vbroadcast_sub(&S_tile, S_max_tile);
#define scaled_exp(x) native_vexp2(x *scale)
 tile_elementwise(S_tile, scaled_exp);
#ifdef ALT_MAX
 intel_work_group_barrier_wait(CLK_LOCAL_MEM_FENCE);
 s_sum_tile_type S_max_tile1;
 tile_copy(S_max_tile, S_max_tile1);
 tile_load_full(&S_max_tile, S_max_slm, ugemm_kq_wg_tile_n, sg_j0_kq, 0);
#define binary_exp_neg(x, y) native_vexp2(scale *((x) - (y)))
 tile_binary(S_max_tile1, S_max_tile, binary_exp_neg);
 tile_vbroadcast_mul(&S_tile, S_max_tile1);
#endif
 s_sum_tile_type S_sum_tile1;
 tile_fill(S_sum_tile1, 0.0f);
 tile_vreduce_add(S_tile, &S_sum_tile1);
 s_tile_type_half2 S_tile_half2;
 tile_copy_to_half2(S_tile, S_tile_half2);
 tile_store_t_sys_src2(S_tile_half2, (local uint *)S_slm,
 ugemm_vs_sg_tile_n, ugemm_kq_wg_tile_m / 2, sg_i0_kq / 2,
 sg_j0_kq);
 intel_work_group_barrier_arrive(CLK_LOCAL_MEM_FENCE);
 if (!first) {
#define binary_exp_sub(x, y) native_vexp2(scale *((x) - (y)))
#define binary_mul(x, y) ((x) * (y))
 tile_binary(S_max_tile_old, S_max_tile, binary_exp_sub);
 tile_binary(S_sum_tile, S_max_tile_old, binary_mul);
 a_scale_tile_type A_scale_tile;
#if ugemm_kq_wg_tile_n == ugemm_vs_wg_tile_n && ugemm_kq_sg_tile_n == ugemm_vs_sg_tile_n
 tile_copy(S_max_tile_old, A_scale_tile);
#elif ugemm_kq_wg_tile_n == ugemm_vs_wg_tile_n && (ugemm_kq_sg_tile_n % ugemm_vs_sg_tile_n) == 0
 tile_rselect(&A_scale_tile, S_max_tile_old,
 sg_j_vs % (ugemm_kq_sg_tile_n / ugemm_vs_sg_tile_n));
#else
#error unimplemented
#endif
 tile_hbroadcast_mul(&A_tile, A_scale_tile);
 }
 tile_binary(S_sum_tile, S_sum_tile1, binary_add);
 tile_copy(S_max_tile, S_max_tile_old);
 if (last) {
 tile_store_full(S_sum_tile, S_sum_slm, ugemm_kq_wg_tile_n, sg_j0_kq,
 sg_i_kq);
 }
#ifdef PREFETCH_K
 if (!last) {
#if TRANSPOSE_K
 const uint stride_k = ldk;
#else
 const uint stride_k = 1;
#endif
 cooperative_prefetch_2d_k(K + (k0 + ugemm_kq_wg_tile_m) * stride_k,
 k - k0 - ugemm_kq_wg_tile_m, d, ugemm_kq_wg_tile_m,
 PREFETCH_D_MAX, ldk, sg_ij, sg_per_wg, SUBGROUP_SIZE,
 LSC_LDCC_L1C_L3C);
 }
#endif
#if WITH_ATTN_MASK && defined(PREFETCH_MASK)
 if (!last) {
 cooperative_prefetch_2d(msk + k0 + ugemm_kq_wg_tile_m + sg_i0_kq + (sg_j0_kq + wg_j0) * q,
 ugemm_kq_sg_tile_m, ugemm_kq_sg_tile_n, 0, 0, 1, SUBGROUP_SIZE,
 LSC_LDCC_L1UC_L3C);
 }
#endif
 intel_work_group_barrier_wait(CLK_LOCAL_MEM_FENCE);
 if (last && need_sum_barrier)
 intel_work_group_barrier_arrive(CLK_LOCAL_MEM_FENCE);
 int k_chunk = min(k - k0, ugemm_kq_wg_tile_m);
 a_tile_type A_tile1 = ugemm_vs(V, ldv, S_slm, ugemm_kq_wg_tile_m, d,
 ugemm_kq_wg_tile_n, k_chunk, 0, 0, 0, sg_i_vs, sg_j_vs,
 (local char *)ugemm_slm);
 V += ldv * ugemm_kq_wg_tile_m;
 tile_binary(A_tile, A_tile1, binary_add);
 }
 if (need_sum_barrier) intel_work_group_barrier_wait(CLK_LOCAL_MEM_FENCE);
 a_scale_tile_type A_scale_tile, A_scale_tile_load;
 tile_fill(A_scale_tile, 0.0f);
#pragma unroll
 for (uint sg1 = 0; sg1 < ugemm_kq_sg_per_wg_m; sg1++) {
 tile_load_full(&A_scale_tile_load, S_sum_slm, ugemm_kq_wg_tile_n,
 ugemm_vs_sg_tile_n * sg_j_vs, sg1);
 tile_binary(A_scale_tile, A_scale_tile_load, binary_add);
 }
 tile_elementwise(A_scale_tile, native_vrecip);
 tile_hbroadcast_mul(&A_tile, A_scale_tile);
 a_tile_type_half A_tile_half;
 tile_copy_reblock(A_tile, &A_tile_half);
 uint sg_i0_vs = sg_i_vs * ugemm_vs_sg_tile_m;
 uint sg_j0_vs = sg_j_vs * ugemm_vs_sg_tile_n + wg_j0;
#ifdef BLOCK_2D_A
 tile_store_block2d(A_tile_half, A, d, q, lda, sg_i0_vs, sg_j0_vs);
#elif defined(BLOCK_A)
 tile_store_block_rem_q(A_tile_half, A, q, lda, sg_i0_vs, sg_j0_vs);
#else
 tile_store(A_tile_half, A, d, q, lda, sg_i0_vs, sg_j0_vs);
#endif
}
#ifdef MAX
#undef MAX
#endif
#ifdef DIV_UP
#undef DIV_UP
#endif
#ifdef sg_per_wg
#undef sg_per_wg
#endif
#ifdef q_tile_sg_n
#undef q_tile_sg_n
#endif
#ifdef cooperative_prefetch_2d_maybe_rem
#undef cooperative_prefetch_2d_maybe_rem
#endif
#ifdef cooperative_prefetch_2d_maybe_rem
#undef cooperative_prefetch_2d_maybe_rem
#endif
#ifdef cooperative_prefetch_2d_k
#undef cooperative_prefetch_2d_k
#endif
#ifdef cooperative_prefetch_2d_k
#undef cooperative_prefetch_2d_k
#endif
#ifdef tile_load_block_rem_q
#undef tile_load_block_rem_q
#endif
#ifdef tile_store_block_rem_q
#undef tile_store_block_rem_q
#endif
#ifdef tile_load_block_rem_q
#undef tile_load_block_rem_q
#endif
#ifdef tile_store_block_rem_q
#undef tile_store_block_rem_q
#endif
#ifdef binary_add
#undef binary_add
#endif
#ifdef Q_slm_size
#undef Q_slm_size
#endif
#ifdef S_slm_size
#undef S_slm_size
#endif
#ifdef S_sum_slm_size
#undef S_sum_slm_size
#endif
#ifdef S_max_slm_size
#undef S_max_slm_size
#endif
#ifdef ugemm_slm_size
#undef ugemm_slm_size
#endif
#ifdef unscale
#undef unscale
#endif
#ifdef scaled_exp
#undef scaled_exp
#endif
#ifdef binary_exp_neg
#undef binary_exp_neg
#endif
#ifdef binary_exp_sub
#undef binary_exp_sub
#endif
#ifdef binary_mul
#undef binary_mul
#endif
#undef KERNEL
#undef KERNEL_ID
#undef FUNC
#undef FUNC_CALL
#undef CONST_ARRAY_DECL
#undef CONST_ARRAY_REF
#ifdef FP64_SUPPORTED
#undef FP64_SUPPORTED
#endif
#ifdef FP16_SUPPORTED
#undef FP16_SUPPORTED
#endif
#ifdef FP16_UNIT_USED
#undef FP16_UNIT_USED
#endif
#ifdef INT8_UNIT_USED
#undef INT8_UNIT_USED
#endif
#ifdef INT32_UNIT_USED
#undef INT32_UNIT_USED
#endif
#ifdef INT64_UNIT_USED
#undef INT64_UNIT_USED
#endif
#ifdef UINT8_UNIT_USED
#undef UINT8_UNIT_USED
#endif
#ifdef UINT32_UNIT_USED
#undef UINT32_UNIT_USED
#endif
#ifdef UNIT_TYPE
#undef UNIT_TYPE
#endif
#ifdef UNIT_VAL_MAX
#undef UNIT_VAL_MAX
#endif
#ifdef UNIT_VAL_MIN
#undef UNIT_VAL_MIN
#endif
#ifdef UNIT_VAL_ONE
#undef UNIT_VAL_ONE
#endif
#ifdef UNIT_VAL_ZERO
#undef UNIT_VAL_ZERO
#endif
#ifdef TO_UNIT_TYPE
#undef TO_UNIT_TYPE
#endif
#ifdef TO_UNIT_TYPE_SAT
#undef TO_UNIT_TYPE_SAT
#endif
#ifdef AS_UNIT_TYPE
#undef AS_UNIT_TYPE
#endif
#ifdef UNIT_MAX_FUNC
#undef UNIT_MAX_FUNC
#endif
#ifdef UNIT_MIN_FUNC
#undef UNIT_MIN_FUNC
#endif
#ifdef UNIT_ABS_FUNC
#undef UNIT_ABS_FUNC
#endif
#ifdef UNIT_TYPE_SIZE
#undef UNIT_TYPE_SIZE
#endif
#ifdef UNIT_IS_FP
#undef UNIT_IS_FP
#endif
#ifdef NL_M
#undef NL_M
#endif
#ifdef NL_N
#undef NL_N
#endif
#ifdef ACTIVATION_FUNC_TYPE
#undef ACTIVATION_FUNC_TYPE
#endif
#ifdef ACTIVATION_FUNC_VAL_MAX
#undef ACTIVATION_FUNC_VAL_MAX
#endif
#ifdef ACTIVATION_FUNC_VAL_MIN
#undef ACTIVATION_FUNC_VAL_MIN
#endif
#ifdef ACTIVATION_FUNC_VAL_ONE
#undef ACTIVATION_FUNC_VAL_ONE
#endif
#ifdef ACTIVATION_FUNC_VAL_ZERO
#undef ACTIVATION_FUNC_VAL_ZERO
#endif
#ifdef TO_ACTIVATION_FUNC_TYPE
#undef TO_ACTIVATION_FUNC_TYPE
#endif
#ifdef TO_ACTIVATION_FUNC_TYPE_SAT
#undef TO_ACTIVATION_FUNC_TYPE_SAT
#endif
#ifdef AS_ACTIVATION_FUNC_TYPE
#undef AS_ACTIVATION_FUNC_TYPE
#endif
#ifdef ACTIVATION_FUNC_MAX_FUNC
#undef ACTIVATION_FUNC_MAX_FUNC
#endif
#ifdef ACTIVATION_FUNC_MIN_FUNC
#undef ACTIVATION_FUNC_MIN_FUNC
#endif
#ifdef ACTIVATION_FUNC_ABS_FUNC
#undef ACTIVATION_FUNC_ABS_FUNC
#endif
#ifdef ACTIVATION_FUNC_TYPE_SIZE
#undef ACTIVATION_FUNC_TYPE_SIZE
#endif
#ifdef ACTIVATION_FUNC_IS_FP
#undef ACTIVATION_FUNC_IS_FP
#endif
#ifdef ACTIVATION_PARAMS
#undef ACTIVATION_PARAMS
#endif
#ifdef ACTIVATION_FUNC
#undef ACTIVATION_FUNC
#endif
#ifdef ACTIVATION
#undef ACTIVATION
#endif
#ifdef INPUT0_SIZE_X
#undef INPUT0_SIZE_X
#endif
#ifdef INPUT0_SIZE_Y
#undef INPUT0_SIZE_Y
#endif
#ifdef INPUT0_SIZE_Z
#undef INPUT0_SIZE_Z
#endif
#ifdef INPUT0_SIZE_W
#undef INPUT0_SIZE_W
#endif
#ifdef INPUT0_SIZE_U
#undef INPUT0_SIZE_U
#endif
#ifdef INPUT0_SIZE_V
#undef INPUT0_SIZE_V
#endif
#ifdef INPUT0_FEATURE_NUM
#undef INPUT0_FEATURE_NUM
#endif
#ifdef INPUT0_BATCH_NUM
#undef INPUT0_BATCH_NUM
#endif
#ifdef INPUT0_PAD_BEFORE_SIZE_X
#undef INPUT0_PAD_BEFORE_SIZE_X
#endif
#ifdef INPUT0_PAD_BEFORE_SIZE_Y
#undef INPUT0_PAD_BEFORE_SIZE_Y
#endif
#ifdef INPUT0_PAD_BEFORE_SIZE_Z
#undef INPUT0_PAD_BEFORE_SIZE_Z
#endif
#ifdef INPUT0_PAD_BEFORE_SIZE_W
#undef INPUT0_PAD_BEFORE_SIZE_W
#endif
#ifdef INPUT0_PAD_BEFORE_SIZE_U
#undef INPUT0_PAD_BEFORE_SIZE_U
#endif
#ifdef INPUT0_PAD_BEFORE_SIZE_V
#undef INPUT0_PAD_BEFORE_SIZE_V
#endif
#ifdef INPUT0_PAD_BEFORE_FEATURE_NUM
#undef INPUT0_PAD_BEFORE_FEATURE_NUM
#endif
#ifdef INPUT0_PAD_BEFORE_BATCH_NUM
#undef INPUT0_PAD_BEFORE_BATCH_NUM
#endif
#ifdef INPUT0_PAD_AFTER_SIZE_X
#undef INPUT0_PAD_AFTER_SIZE_X
#endif
#ifdef INPUT0_PAD_AFTER_SIZE_Y
#undef INPUT0_PAD_AFTER_SIZE_Y
#endif
#ifdef INPUT0_PAD_AFTER_SIZE_Z
#undef INPUT0_PAD_AFTER_SIZE_Z
#endif
#ifdef INPUT0_PAD_AFTER_SIZE_W
#undef INPUT0_PAD_AFTER_SIZE_W
#endif
#ifdef INPUT0_PAD_AFTER_SIZE_U
#undef INPUT0_PAD_AFTER_SIZE_U
#endif
#ifdef INPUT0_PAD_AFTER_SIZE_V
#undef INPUT0_PAD_AFTER_SIZE_V
#endif
#ifdef INPUT0_PAD_AFTER_FEATURE_NUM
#undef INPUT0_PAD_AFTER_FEATURE_NUM
#endif
#ifdef INPUT0_PAD_AFTER_BATCH_NUM
#undef INPUT0_PAD_AFTER_BATCH_NUM
#endif
#ifdef INPUT0_X_PITCH
#undef INPUT0_X_PITCH
#endif
#ifdef INPUT0_Y_PITCH
#undef INPUT0_Y_PITCH
#endif
#ifdef INPUT0_Z_PITCH
#undef INPUT0_Z_PITCH
#endif
#ifdef INPUT0_W_PITCH
#undef INPUT0_W_PITCH
#endif
#ifdef INPUT0_U_PITCH
#undef INPUT0_U_PITCH
#endif
#ifdef INPUT0_V_PITCH
#undef INPUT0_V_PITCH
#endif
#ifdef INPUT0_FEATURE_PITCH
#undef INPUT0_FEATURE_PITCH
#endif
#ifdef INPUT0_BATCH_PITCH
#undef INPUT0_BATCH_PITCH
#endif
#ifdef INPUT0_GET_INDEX_SAFE
#undef INPUT0_GET_INDEX_SAFE
#endif
#ifdef INPUT0_GET_INDEX
#undef INPUT0_GET_INDEX
#endif
#ifdef INPUT0_GET_INDEX_RAW
#undef INPUT0_GET_INDEX_RAW
#endif
#ifdef INPUT0_VIEW_OFFSET
#undef INPUT0_VIEW_OFFSET
#endif
#ifdef INPUT0_LENGTH
#undef INPUT0_LENGTH
#endif
#ifdef INPUT0_DIMS
#undef INPUT0_DIMS
#endif
#ifdef INPUT0_SIMPLE
#undef INPUT0_SIMPLE
#endif
#ifdef INPUT0_GROUPED
#undef INPUT0_GROUPED
#endif
#ifdef INPUT0_LAYOUT_BFYX
#undef INPUT0_LAYOUT_BFYX
#endif
#ifdef INPUT0_TYPE
#undef INPUT0_TYPE
#endif
#ifdef INPUT0_VAL_MAX
#undef INPUT0_VAL_MAX
#endif
#ifdef INPUT0_VAL_MIN
#undef INPUT0_VAL_MIN
#endif
#ifdef INPUT0_VAL_ONE
#undef INPUT0_VAL_ONE
#endif
#ifdef INPUT0_VAL_ZERO
#undef INPUT0_VAL_ZERO
#endif
#ifdef TO_INPUT0_TYPE
#undef TO_INPUT0_TYPE
#endif
#ifdef TO_INPUT0_TYPE_SAT
#undef TO_INPUT0_TYPE_SAT
#endif
#ifdef AS_INPUT0_TYPE
#undef AS_INPUT0_TYPE
#endif
#ifdef INPUT0_MAX_FUNC
#undef INPUT0_MAX_FUNC
#endif
#ifdef INPUT0_MIN_FUNC
#undef INPUT0_MIN_FUNC
#endif
#ifdef INPUT0_ABS_FUNC
#undef INPUT0_ABS_FUNC
#endif
#ifdef INPUT0_TYPE_SIZE
#undef INPUT0_TYPE_SIZE
#endif
#ifdef INPUT0_IS_FP
#undef INPUT0_IS_FP
#endif
#ifdef INPUT0_OFFSET
#undef INPUT0_OFFSET
#endif
#ifdef INPUT0_PAD_BEFORE
#undef INPUT0_PAD_BEFORE
#endif
#ifdef INPUT0_PAD_AFTER
#undef INPUT0_PAD_AFTER
#endif
#ifdef INPUT1_SIZE_X
#undef INPUT1_SIZE_X
#endif
#ifdef INPUT1_SIZE_Y
#undef INPUT1_SIZE_Y
#endif
#ifdef INPUT1_SIZE_Z
#undef INPUT1_SIZE_Z
#endif
#ifdef INPUT1_SIZE_W
#undef INPUT1_SIZE_W
#endif
#ifdef INPUT1_SIZE_U
#undef INPUT1_SIZE_U
#endif
#ifdef INPUT1_SIZE_V
#undef INPUT1_SIZE_V
#endif
#ifdef INPUT1_FEATURE_NUM
#undef INPUT1_FEATURE_NUM
#endif
#ifdef INPUT1_BATCH_NUM
#undef INPUT1_BATCH_NUM
#endif
#ifdef INPUT1_PAD_BEFORE_SIZE_X
#undef INPUT1_PAD_BEFORE_SIZE_X
#endif
#ifdef INPUT1_PAD_BEFORE_SIZE_Y
#undef INPUT1_PAD_BEFORE_SIZE_Y
#endif
#ifdef INPUT1_PAD_BEFORE_SIZE_Z
#undef INPUT1_PAD_BEFORE_SIZE_Z
#endif
#ifdef INPUT1_PAD_BEFORE_SIZE_W
#undef INPUT1_PAD_BEFORE_SIZE_W
#endif
#ifdef INPUT1_PAD_BEFORE_SIZE_U
#undef INPUT1_PAD_BEFORE_SIZE_U
#endif
#ifdef INPUT1_PAD_BEFORE_SIZE_V
#undef INPUT1_PAD_BEFORE_SIZE_V
#endif
#ifdef INPUT1_PAD_BEFORE_FEATURE_NUM
#undef INPUT1_PAD_BEFORE_FEATURE_NUM
#endif
#ifdef INPUT1_PAD_BEFORE_BATCH_NUM
#undef INPUT1_PAD_BEFORE_BATCH_NUM
#endif
#ifdef INPUT1_PAD_AFTER_SIZE_X
#undef INPUT1_PAD_AFTER_SIZE_X
#endif
#ifdef INPUT1_PAD_AFTER_SIZE_Y
#undef INPUT1_PAD_AFTER_SIZE_Y
#endif
#ifdef INPUT1_PAD_AFTER_SIZE_Z
#undef INPUT1_PAD_AFTER_SIZE_Z
#endif
#ifdef INPUT1_PAD_AFTER_SIZE_W
#undef INPUT1_PAD_AFTER_SIZE_W
#endif
#ifdef INPUT1_PAD_AFTER_SIZE_U
#undef INPUT1_PAD_AFTER_SIZE_U
#endif
#ifdef INPUT1_PAD_AFTER_SIZE_V
#undef INPUT1_PAD_AFTER_SIZE_V
#endif
#ifdef INPUT1_PAD_AFTER_FEATURE_NUM
#undef INPUT1_PAD_AFTER_FEATURE_NUM
#endif
#ifdef INPUT1_PAD_AFTER_BATCH_NUM
#undef INPUT1_PAD_AFTER_BATCH_NUM
#endif
#ifdef INPUT1_X_PITCH
#undef INPUT1_X_PITCH
#endif
#ifdef INPUT1_Y_PITCH
#undef INPUT1_Y_PITCH
#endif
#ifdef INPUT1_Z_PITCH
#undef INPUT1_Z_PITCH
#endif
#ifdef INPUT1_W_PITCH
#undef INPUT1_W_PITCH
#endif
#ifdef INPUT1_U_PITCH
#undef INPUT1_U_PITCH
#endif
#ifdef INPUT1_V_PITCH
#undef INPUT1_V_PITCH
#endif
#ifdef INPUT1_FEATURE_PITCH
#undef INPUT1_FEATURE_PITCH
#endif
#ifdef INPUT1_BATCH_PITCH
#undef INPUT1_BATCH_PITCH
#endif
#ifdef INPUT1_GET_INDEX_SAFE
#undef INPUT1_GET_INDEX_SAFE
#endif
#ifdef INPUT1_GET_INDEX
#undef INPUT1_GET_INDEX
#endif
#ifdef INPUT1_GET_INDEX_RAW
#undef INPUT1_GET_INDEX_RAW
#endif
#ifdef INPUT1_VIEW_OFFSET
#undef INPUT1_VIEW_OFFSET
#endif
#ifdef INPUT1_LENGTH
#undef INPUT1_LENGTH
#endif
#ifdef INPUT1_DIMS
#undef INPUT1_DIMS
#endif
#ifdef INPUT1_SIMPLE
#undef INPUT1_SIMPLE
#endif
#ifdef INPUT1_GROUPED
#undef INPUT1_GROUPED
#endif
#ifdef INPUT1_LAYOUT_BFYX
#undef INPUT1_LAYOUT_BFYX
#endif
#ifdef INPUT1_TYPE
#undef INPUT1_TYPE
#endif
#ifdef INPUT1_VAL_MAX
#undef INPUT1_VAL_MAX
#endif
#ifdef INPUT1_VAL_MIN
#undef INPUT1_VAL_MIN
#endif
#ifdef INPUT1_VAL_ONE
#undef INPUT1_VAL_ONE
#endif
#ifdef INPUT1_VAL_ZERO
#undef INPUT1_VAL_ZERO
#endif
#ifdef TO_INPUT1_TYPE
#undef TO_INPUT1_TYPE
#endif
#ifdef TO_INPUT1_TYPE_SAT
#undef TO_INPUT1_TYPE_SAT
#endif
#ifdef AS_INPUT1_TYPE
#undef AS_INPUT1_TYPE
#endif
#ifdef INPUT1_MAX_FUNC
#undef INPUT1_MAX_FUNC
#endif
#ifdef INPUT1_MIN_FUNC
#undef INPUT1_MIN_FUNC
#endif
#ifdef INPUT1_ABS_FUNC
#undef INPUT1_ABS_FUNC
#endif
#ifdef INPUT1_TYPE_SIZE
#undef INPUT1_TYPE_SIZE
#endif
#ifdef INPUT1_IS_FP
#undef INPUT1_IS_FP
#endif
#ifdef INPUT1_OFFSET
#undef INPUT1_OFFSET
#endif
#ifdef INPUT1_PAD_BEFORE
#undef INPUT1_PAD_BEFORE
#endif
#ifdef INPUT1_PAD_AFTER
#undef INPUT1_PAD_AFTER
#endif
#ifdef INPUT2_SIZE_X
#undef INPUT2_SIZE_X
#endif
#ifdef INPUT2_SIZE_Y
#undef INPUT2_SIZE_Y
#endif
#ifdef INPUT2_SIZE_Z
#undef INPUT2_SIZE_Z
#endif
#ifdef INPUT2_SIZE_W
#undef INPUT2_SIZE_W
#endif
#ifdef INPUT2_SIZE_U
#undef INPUT2_SIZE_U
#endif
#ifdef INPUT2_SIZE_V
#undef INPUT2_SIZE_V
#endif
#ifdef INPUT2_FEATURE_NUM
#undef INPUT2_FEATURE_NUM
#endif
#ifdef INPUT2_BATCH_NUM
#undef INPUT2_BATCH_NUM
#endif
#ifdef INPUT2_PAD_BEFORE_SIZE_X
#undef INPUT2_PAD_BEFORE_SIZE_X
#endif
#ifdef INPUT2_PAD_BEFORE_SIZE_Y
#undef INPUT2_PAD_BEFORE_SIZE_Y
#endif
#ifdef INPUT2_PAD_BEFORE_SIZE_Z
#undef INPUT2_PAD_BEFORE_SIZE_Z
#endif
#ifdef INPUT2_PAD_BEFORE_SIZE_W
#undef INPUT2_PAD_BEFORE_SIZE_W
#endif
#ifdef INPUT2_PAD_BEFORE_SIZE_U
#undef INPUT2_PAD_BEFORE_SIZE_U
#endif
#ifdef INPUT2_PAD_BEFORE_SIZE_V
#undef INPUT2_PAD_BEFORE_SIZE_V
#endif
#ifdef INPUT2_PAD_BEFORE_FEATURE_NUM
#undef INPUT2_PAD_BEFORE_FEATURE_NUM
#endif
#ifdef INPUT2_PAD_BEFORE_BATCH_NUM
#undef INPUT2_PAD_BEFORE_BATCH_NUM
#endif
#ifdef INPUT2_PAD_AFTER_SIZE_X
#undef INPUT2_PAD_AFTER_SIZE_X
#endif
#ifdef INPUT2_PAD_AFTER_SIZE_Y
#undef INPUT2_PAD_AFTER_SIZE_Y
#endif
#ifdef INPUT2_PAD_AFTER_SIZE_Z
#undef INPUT2_PAD_AFTER_SIZE_Z
#endif
#ifdef INPUT2_PAD_AFTER_SIZE_W
#undef INPUT2_PAD_AFTER_SIZE_W
#endif
#ifdef INPUT2_PAD_AFTER_SIZE_U
#undef INPUT2_PAD_AFTER_SIZE_U
#endif
#ifdef INPUT2_PAD_AFTER_SIZE_V
#undef INPUT2_PAD_AFTER_SIZE_V
#endif
#ifdef INPUT2_PAD_AFTER_FEATURE_NUM
#undef INPUT2_PAD_AFTER_FEATURE_NUM
#endif
#ifdef INPUT2_PAD_AFTER_BATCH_NUM
#undef INPUT2_PAD_AFTER_BATCH_NUM
#endif
#ifdef INPUT2_X_PITCH
#undef INPUT2_X_PITCH
#endif
#ifdef INPUT2_Y_PITCH
#undef INPUT2_Y_PITCH
#endif
#ifdef INPUT2_Z_PITCH
#undef INPUT2_Z_PITCH
#endif
#ifdef INPUT2_W_PITCH
#undef INPUT2_W_PITCH
#endif
#ifdef INPUT2_U_PITCH
#undef INPUT2_U_PITCH
#endif
#ifdef INPUT2_V_PITCH
#undef INPUT2_V_PITCH
#endif
#ifdef INPUT2_FEATURE_PITCH
#undef INPUT2_FEATURE_PITCH
#endif
#ifdef INPUT2_BATCH_PITCH
#undef INPUT2_BATCH_PITCH
#endif
#ifdef INPUT2_GET_INDEX_SAFE
#undef INPUT2_GET_INDEX_SAFE
#endif
#ifdef INPUT2_GET_INDEX
#undef INPUT2_GET_INDEX
#endif
#ifdef INPUT2_GET_INDEX_RAW
#undef INPUT2_GET_INDEX_RAW
#endif
#ifdef INPUT2_VIEW_OFFSET
#undef INPUT2_VIEW_OFFSET
#endif
#ifdef INPUT2_LENGTH
#undef INPUT2_LENGTH
#endif
#ifdef INPUT2_DIMS
#undef INPUT2_DIMS
#endif
#ifdef INPUT2_SIMPLE
#undef INPUT2_SIMPLE
#endif
#ifdef INPUT2_GROUPED
#undef INPUT2_GROUPED
#endif
#ifdef INPUT2_LAYOUT_BFYX
#undef INPUT2_LAYOUT_BFYX
#endif
#ifdef INPUT2_TYPE
#undef INPUT2_TYPE
#endif
#ifdef INPUT2_VAL_MAX
#undef INPUT2_VAL_MAX
#endif
#ifdef INPUT2_VAL_MIN
#undef INPUT2_VAL_MIN
#endif
#ifdef INPUT2_VAL_ONE
#undef INPUT2_VAL_ONE
#endif
#ifdef INPUT2_VAL_ZERO
#undef INPUT2_VAL_ZERO
#endif
#ifdef TO_INPUT2_TYPE
#undef TO_INPUT2_TYPE
#endif
#ifdef TO_INPUT2_TYPE_SAT
#undef TO_INPUT2_TYPE_SAT
#endif
#ifdef AS_INPUT2_TYPE
#undef AS_INPUT2_TYPE
#endif
#ifdef INPUT2_MAX_FUNC
#undef INPUT2_MAX_FUNC
#endif
#ifdef INPUT2_MIN_FUNC
#undef INPUT2_MIN_FUNC
#endif
#ifdef INPUT2_ABS_FUNC
#undef INPUT2_ABS_FUNC
#endif
#ifdef INPUT2_TYPE_SIZE
#undef INPUT2_TYPE_SIZE
#endif
#ifdef INPUT2_IS_FP
#undef INPUT2_IS_FP
#endif
#ifdef INPUT2_OFFSET
#undef INPUT2_OFFSET
#endif
#ifdef INPUT2_PAD_BEFORE
#undef INPUT2_PAD_BEFORE
#endif
#ifdef INPUT2_PAD_AFTER
#undef INPUT2_PAD_AFTER
#endif
#ifdef INPUT3_SIZE_X
#undef INPUT3_SIZE_X
#endif
#ifdef INPUT3_SIZE_Y
#undef INPUT3_SIZE_Y
#endif
#ifdef INPUT3_SIZE_Z
#undef INPUT3_SIZE_Z
#endif
#ifdef INPUT3_SIZE_W
#undef INPUT3_SIZE_W
#endif
#ifdef INPUT3_SIZE_U
#undef INPUT3_SIZE_U
#endif
#ifdef INPUT3_SIZE_V
#undef INPUT3_SIZE_V
#endif
#ifdef INPUT3_FEATURE_NUM
#undef INPUT3_FEATURE_NUM
#endif
#ifdef INPUT3_BATCH_NUM
#undef INPUT3_BATCH_NUM
#endif
#ifdef INPUT3_PAD_BEFORE_SIZE_X
#undef INPUT3_PAD_BEFORE_SIZE_X
#endif
#ifdef INPUT3_PAD_BEFORE_SIZE_Y
#undef INPUT3_PAD_BEFORE_SIZE_Y
#endif
#ifdef INPUT3_PAD_BEFORE_SIZE_Z
#undef INPUT3_PAD_BEFORE_SIZE_Z
#endif
#ifdef INPUT3_PAD_BEFORE_SIZE_W
#undef INPUT3_PAD_BEFORE_SIZE_W
#endif
#ifdef INPUT3_PAD_BEFORE_SIZE_U
#undef INPUT3_PAD_BEFORE_SIZE_U
#endif
#ifdef INPUT3_PAD_BEFORE_SIZE_V
#undef INPUT3_PAD_BEFORE_SIZE_V
#endif
#ifdef INPUT3_PAD_BEFORE_FEATURE_NUM
#undef INPUT3_PAD_BEFORE_FEATURE_NUM
#endif
#ifdef INPUT3_PAD_BEFORE_BATCH_NUM
#undef INPUT3_PAD_BEFORE_BATCH_NUM
#endif
#ifdef INPUT3_PAD_AFTER_SIZE_X
#undef INPUT3_PAD_AFTER_SIZE_X
#endif
#ifdef INPUT3_PAD_AFTER_SIZE_Y
#undef INPUT3_PAD_AFTER_SIZE_Y
#endif
#ifdef INPUT3_PAD_AFTER_SIZE_Z
#undef INPUT3_PAD_AFTER_SIZE_Z
#endif
#ifdef INPUT3_PAD_AFTER_SIZE_W
#undef INPUT3_PAD_AFTER_SIZE_W
#endif
#ifdef INPUT3_PAD_AFTER_SIZE_U
#undef INPUT3_PAD_AFTER_SIZE_U
#endif
#ifdef INPUT3_PAD_AFTER_SIZE_V
#undef INPUT3_PAD_AFTER_SIZE_V
#endif
#ifdef INPUT3_PAD_AFTER_FEATURE_NUM
#undef INPUT3_PAD_AFTER_FEATURE_NUM
#endif
#ifdef INPUT3_PAD_AFTER_BATCH_NUM
#undef INPUT3_PAD_AFTER_BATCH_NUM
#endif
#ifdef INPUT3_X_PITCH
#undef INPUT3_X_PITCH
#endif
#ifdef INPUT3_Y_PITCH
#undef INPUT3_Y_PITCH
#endif
#ifdef INPUT3_Z_PITCH
#undef INPUT3_Z_PITCH
#endif
#ifdef INPUT3_W_PITCH
#undef INPUT3_W_PITCH
#endif
#ifdef INPUT3_U_PITCH
#undef INPUT3_U_PITCH
#endif
#ifdef INPUT3_V_PITCH
#undef INPUT3_V_PITCH
#endif
#ifdef INPUT3_FEATURE_PITCH
#undef INPUT3_FEATURE_PITCH
#endif
#ifdef INPUT3_BATCH_PITCH
#undef INPUT3_BATCH_PITCH
#endif
#ifdef INPUT3_GET_INDEX_SAFE
#undef INPUT3_GET_INDEX_SAFE
#endif
#ifdef INPUT3_GET_INDEX
#undef INPUT3_GET_INDEX
#endif
#ifdef INPUT3_GET_INDEX_RAW
#undef INPUT3_GET_INDEX_RAW
#endif
#ifdef INPUT3_VIEW_OFFSET
#undef INPUT3_VIEW_OFFSET
#endif
#ifdef INPUT3_LENGTH
#undef INPUT3_LENGTH
#endif
#ifdef INPUT3_DIMS
#undef INPUT3_DIMS
#endif
#ifdef INPUT3_SIMPLE
#undef INPUT3_SIMPLE
#endif
#ifdef INPUT3_GROUPED
#undef INPUT3_GROUPED
#endif
#ifdef INPUT3_LAYOUT_BFYX
#undef INPUT3_LAYOUT_BFYX
#endif
#ifdef INPUT3_TYPE
#undef INPUT3_TYPE
#endif
#ifdef INPUT3_VAL_MAX
#undef INPUT3_VAL_MAX
#endif
#ifdef INPUT3_VAL_MIN
#undef INPUT3_VAL_MIN
#endif
#ifdef INPUT3_VAL_ONE
#undef INPUT3_VAL_ONE
#endif
#ifdef INPUT3_VAL_ZERO
#undef INPUT3_VAL_ZERO
#endif
#ifdef TO_INPUT3_TYPE
#undef TO_INPUT3_TYPE
#endif
#ifdef TO_INPUT3_TYPE_SAT
#undef TO_INPUT3_TYPE_SAT
#endif
#ifdef AS_INPUT3_TYPE
#undef AS_INPUT3_TYPE
#endif
#ifdef INPUT3_MAX_FUNC
#undef INPUT3_MAX_FUNC
#endif
#ifdef INPUT3_MIN_FUNC
#undef INPUT3_MIN_FUNC
#endif
#ifdef INPUT3_ABS_FUNC
#undef INPUT3_ABS_FUNC
#endif
#ifdef INPUT3_TYPE_SIZE
#undef INPUT3_TYPE_SIZE
#endif
#ifdef INPUT3_IS_FP
#undef INPUT3_IS_FP
#endif
#ifdef INPUT3_OFFSET
#undef INPUT3_OFFSET
#endif
#ifdef INPUT3_PAD_BEFORE
#undef INPUT3_PAD_BEFORE
#endif
#ifdef INPUT3_PAD_AFTER
#undef INPUT3_PAD_AFTER
#endif
#ifdef OUTPUT_SIZE_X
#undef OUTPUT_SIZE_X
#endif
#ifdef OUTPUT_SIZE_Y
#undef OUTPUT_SIZE_Y
#endif
#ifdef OUTPUT_SIZE_Z
#undef OUTPUT_SIZE_Z
#endif
#ifdef OUTPUT_SIZE_W
#undef OUTPUT_SIZE_W
#endif
#ifdef OUTPUT_SIZE_U
#undef OUTPUT_SIZE_U
#endif
#ifdef OUTPUT_SIZE_V
#undef OUTPUT_SIZE_V
#endif
#ifdef OUTPUT_FEATURE_NUM
#undef OUTPUT_FEATURE_NUM
#endif
#ifdef OUTPUT_BATCH_NUM
#undef OUTPUT_BATCH_NUM
#endif
#ifdef OUTPUT_PAD_BEFORE_SIZE_X
#undef OUTPUT_PAD_BEFORE_SIZE_X
#endif
#ifdef OUTPUT_PAD_BEFORE_SIZE_Y
#undef OUTPUT_PAD_BEFORE_SIZE_Y
#endif
#ifdef OUTPUT_PAD_BEFORE_SIZE_Z
#undef OUTPUT_PAD_BEFORE_SIZE_Z
#endif
#ifdef OUTPUT_PAD_BEFORE_SIZE_W
#undef OUTPUT_PAD_BEFORE_SIZE_W
#endif
#ifdef OUTPUT_PAD_BEFORE_SIZE_U
#undef OUTPUT_PAD_BEFORE_SIZE_U
#endif
#ifdef OUTPUT_PAD_BEFORE_SIZE_V
#undef OUTPUT_PAD_BEFORE_SIZE_V
#endif
#ifdef OUTPUT_PAD_BEFORE_FEATURE_NUM
#undef OUTPUT_PAD_BEFORE_FEATURE_NUM
#endif
#ifdef OUTPUT_PAD_BEFORE_BATCH_NUM
#undef OUTPUT_PAD_BEFORE_BATCH_NUM
#endif
#ifdef OUTPUT_PAD_AFTER_SIZE_X
#undef OUTPUT_PAD_AFTER_SIZE_X
#endif
#ifdef OUTPUT_PAD_AFTER_SIZE_Y
#undef OUTPUT_PAD_AFTER_SIZE_Y
#endif
#ifdef OUTPUT_PAD_AFTER_SIZE_Z
#undef OUTPUT_PAD_AFTER_SIZE_Z
#endif
#ifdef OUTPUT_PAD_AFTER_SIZE_W
#undef OUTPUT_PAD_AFTER_SIZE_W
#endif
#ifdef OUTPUT_PAD_AFTER_SIZE_U
#undef OUTPUT_PAD_AFTER_SIZE_U
#endif
#ifdef OUTPUT_PAD_AFTER_SIZE_V
#undef OUTPUT_PAD_AFTER_SIZE_V
#endif
#ifdef OUTPUT_PAD_AFTER_FEATURE_NUM
#undef OUTPUT_PAD_AFTER_FEATURE_NUM
#endif
#ifdef OUTPUT_PAD_AFTER_BATCH_NUM
#undef OUTPUT_PAD_AFTER_BATCH_NUM
#endif
#ifdef OUTPUT_X_PITCH
#undef OUTPUT_X_PITCH
#endif
#ifdef OUTPUT_Y_PITCH
#undef OUTPUT_Y_PITCH
#endif
#ifdef OUTPUT_Z_PITCH
#undef OUTPUT_Z_PITCH
#endif
#ifdef OUTPUT_W_PITCH
#undef OUTPUT_W_PITCH
#endif
#ifdef OUTPUT_U_PITCH
#undef OUTPUT_U_PITCH
#endif
#ifdef OUTPUT_V_PITCH
#undef OUTPUT_V_PITCH
#endif
#ifdef OUTPUT_FEATURE_PITCH
#undef OUTPUT_FEATURE_PITCH
#endif
#ifdef OUTPUT_BATCH_PITCH
#undef OUTPUT_BATCH_PITCH
#endif
#ifdef OUTPUT_GET_INDEX_SAFE
#undef OUTPUT_GET_INDEX_SAFE
#endif
#ifdef OUTPUT_GET_INDEX
#undef OUTPUT_GET_INDEX
#endif
#ifdef OUTPUT_GET_INDEX_RAW
#undef OUTPUT_GET_INDEX_RAW
#endif
#ifdef OUTPUT_VIEW_OFFSET
#undef OUTPUT_VIEW_OFFSET
#endif
#ifdef OUTPUT_LENGTH
#undef OUTPUT_LENGTH
#endif
#ifdef OUTPUT_DIMS
#undef OUTPUT_DIMS
#endif
#ifdef OUTPUT_SIMPLE
#undef OUTPUT_SIMPLE
#endif
#ifdef OUTPUT_GROUPED
#undef OUTPUT_GROUPED
#endif
#ifdef OUTPUT_LAYOUT_BFYX
#undef OUTPUT_LAYOUT_BFYX
#endif
#ifdef OUTPUT_TYPE
#undef OUTPUT_TYPE
#endif
#ifdef OUTPUT_VAL_MAX
#undef OUTPUT_VAL_MAX
#endif
#ifdef OUTPUT_VAL_MIN
#undef OUTPUT_VAL_MIN
#endif
#ifdef OUTPUT_VAL_ONE
#undef OUTPUT_VAL_ONE
#endif
#ifdef OUTPUT_VAL_ZERO
#undef OUTPUT_VAL_ZERO
#endif
#ifdef TO_OUTPUT_TYPE
#undef TO_OUTPUT_TYPE
#endif
#ifdef TO_OUTPUT_TYPE_SAT
#undef TO_OUTPUT_TYPE_SAT
#endif
#ifdef AS_OUTPUT_TYPE
#undef AS_OUTPUT_TYPE
#endif
#ifdef OUTPUT_MAX_FUNC
#undef OUTPUT_MAX_FUNC
#endif
#ifdef OUTPUT_MIN_FUNC
#undef OUTPUT_MIN_FUNC
#endif
#ifdef OUTPUT_ABS_FUNC
#undef OUTPUT_ABS_FUNC
#endif
#ifdef OUTPUT_TYPE_SIZE
#undef OUTPUT_TYPE_SIZE
#endif
#ifdef OUTPUT_IS_FP
#undef OUTPUT_IS_FP
#endif
#ifdef OUTPUT_OFFSET
#undef OUTPUT_OFFSET
#endif
#ifdef OUTPUT_PAD_BEFORE
#undef OUTPUT_PAD_BEFORE
#endif
#ifdef OUTPUT_PAD_AFTER
#undef OUTPUT_PAD_AFTER
#endif
#ifdef IS_DYNAMIC
#undef IS_DYNAMIC
#endif
#ifdef OPTIONAL_SHAPE_INFO_ARG
#undef OPTIONAL_SHAPE_INFO_ARG
#endif
#ifdef OPTIONAL_SHAPE_INFO_TENSOR
#undef OPTIONAL_SHAPE_INFO_TENSOR
#endif
#ifdef LayerID
#undef LayerID
#endif
#ifdef D_MAX
#undef D_MAX
#endif
#ifdef SUBGROUP_SIZE
#undef SUBGROUP_SIZE
#endif
#ifdef INVERT_SCALE
#undef INVERT_SCALE
#endif
#ifdef SCALE_DATA_T
#undef SCALE_DATA_T
#endif
#ifdef WITH_ATTN_MASK
#undef WITH_ATTN_MASK
#endif
#ifdef WITH_SCALE
#undef WITH_SCALE
#endif
#ifdef Q_ALIGN
#undef Q_ALIGN
#endif
#ifdef K_ALIGN
#undef K_ALIGN
#endif
#ifdef V_ALIGN
#undef V_ALIGN
#endif
#ifdef A_ALIGN
#undef A_ALIGN
#endif
#ifdef TRANSPOSE_K
#undef TRANSPOSE_K
#endif
#ifdef REMAINDER_K
#undef REMAINDER_K
#endif
#ifdef KV_GROUP_SIZE
#undef KV_GROUP_SIZE
#endif
#ifdef BLOCK_Q
#undef BLOCK_Q
#endif
#ifdef REMAINDER_Q
#undef REMAINDER_Q
#endif
#ifdef QRY_S0
#undef QRY_S0
#endif
#ifdef QRY_S1
#undef QRY_S1
#endif
#ifdef QRY_S2
#undef QRY_S2
#endif
#ifdef QRY_S3
#undef QRY_S3
#endif
#ifdef KEY_S0
#undef KEY_S0
#endif
#ifdef KEY_S1
#undef KEY_S1
#endif
#ifdef KEY_S2
#undef KEY_S2
#endif
#ifdef KEY_S3
#undef KEY_S3
#endif
#ifdef VAL_S0
#undef VAL_S0
#endif
#ifdef VAL_S1
#undef VAL_S1
#endif
#ifdef VAL_S2
#undef VAL_S2
#endif
#ifdef VAL_S3
#undef VAL_S3
#endif
#ifdef DST_S0
#undef DST_S0
#endif
#ifdef DST_S1
#undef DST_S1
#endif
#ifdef DST_S2
#undef DST_S2
#endif
#ifdef DST_S3
#undef DST_S3
#endif
#ifdef QRY_B0
#undef QRY_B0
#endif
#ifdef QRY_SB0
#undef QRY_SB0
#endif
#ifdef QRY_B1
#undef QRY_B1
#endif
#ifdef QRY_SB1
#undef QRY_SB1
#endif
#ifdef QRY_B2
#undef QRY_B2
#endif
#ifdef QRY_SB2
#undef QRY_SB2
#endif
#ifdef QRY_B3
#undef QRY_B3
#endif
#ifdef QRY_SB3
#undef QRY_SB3
#endif
#ifdef KEY_B0
#undef KEY_B0
#endif
#ifdef KEY_SB0
#undef KEY_SB0
#endif
#ifdef KEY_B1
#undef KEY_B1
#endif
#ifdef KEY_SB1
#undef KEY_SB1
#endif
#ifdef KEY_B2
#undef KEY_B2
#endif
#ifdef KEY_SB2
#undef KEY_SB2
#endif
#ifdef KEY_B3
#undef KEY_B3
#endif
#ifdef KEY_SB3
#undef KEY_SB3
#endif
#ifdef VAL_B0
#undef VAL_B0
#endif
#ifdef VAL_SB0
#undef VAL_SB0
#endif
#ifdef VAL_B1
#undef VAL_B1
#endif
#ifdef VAL_SB1
#undef VAL_SB1
#endif
#ifdef VAL_B2
#undef VAL_B2
#endif
#ifdef VAL_SB2
#undef VAL_SB2
#endif
#ifdef VAL_B3
#undef VAL_B3
#endif
#ifdef VAL_SB3
#undef VAL_SB3
#endif
#ifdef DST_B0
#undef DST_B0
#endif
#ifdef DST_SB0
#undef DST_SB0
#endif
#ifdef DST_B1
#undef DST_B1
#endif
#ifdef DST_SB1
#undef DST_SB1
#endif
#ifdef DST_B2
#undef DST_B2
#endif
#ifdef DST_SB2
#undef DST_SB2
#endif
#ifdef DST_B3
#undef DST_B3
#endif
#ifdef DST_SB3
#undef DST_SB3
#endif

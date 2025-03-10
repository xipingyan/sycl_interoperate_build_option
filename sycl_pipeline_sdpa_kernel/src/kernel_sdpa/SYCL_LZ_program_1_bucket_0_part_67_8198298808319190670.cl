
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
// Kernel template: scatter_nd_update_ref 
// Kernel name: scatter_nd_update_ref_13885323453389162897_0_0__sa
#define KERNEL(name) __kernel void scatter_nd_update_ref_13885323453389162897_0_0__sa
#define KERNEL_ID scatter_nd_update_ref_13885323453389162897_0_0__sa
#define FUNC(name)  _##name##_scatter_nd_update_ref_13885323453389162897_0_0__sa
#define FUNC_CALL(name)  _##name##_scatter_nd_update_ref_13885323453389162897_0_0__sa
#define CONST_ARRAY_DECL(name) __constant size_t  _##name##_scatter_nd_update_ref_13885323453389162897_0_0__sa []
#define CONST_ARRAY_REF(name)  _##name##_scatter_nd_update_ref_13885323453389162897_0_0__sa
#define FP64_SUPPORTED 0
#define FP16_SUPPORTED 1
#define FP16_UNIT_USED 1
#define INT8_UNIT_USED 0
#define INT32_UNIT_USED 1
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
#define INPUT0_SIZE_X 1
#define INPUT0_SIZE_Y 1
#define INPUT0_SIZE_Z 1
#define INPUT0_SIZE_W 1
#define INPUT0_SIZE_U 1
#define INPUT0_SIZE_V 1
#define INPUT0_FEATURE_NUM 1
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
#define INPUT0_Y_PITCH 1
#define INPUT0_Z_PITCH (1*1)
#define INPUT0_W_PITCH (1*1*1)
#define INPUT0_U_PITCH (1*1*1*1)
#define INPUT0_V_PITCH (1*1*1*1*1)
#define INPUT0_FEATURE_PITCH (1*1*1*1*1*1)
#define INPUT0_BATCH_PITCH (1*1*1*1*1*1*1)
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
#define INPUT1_SIZE_X 1
#define INPUT1_SIZE_Y 1
#define INPUT1_SIZE_Z 1
#define INPUT1_SIZE_W 1
#define INPUT1_SIZE_U 1
#define INPUT1_SIZE_V 1
#define INPUT1_FEATURE_NUM 1
#define INPUT1_BATCH_NUM (shape_info[8] )
#define INPUT1_PAD_BEFORE_SIZE_X 0
#define INPUT1_PAD_BEFORE_SIZE_Y 0
#define INPUT1_PAD_BEFORE_SIZE_Z 0
#define INPUT1_PAD_BEFORE_SIZE_W 0
#define INPUT1_PAD_BEFORE_SIZE_U 0
#define INPUT1_PAD_BEFORE_SIZE_V 0
#define INPUT1_PAD_BEFORE_FEATURE_NUM 0
#define INPUT1_PAD_BEFORE_BATCH_NUM 0
#define INPUT1_PAD_AFTER_SIZE_X 0
#define INPUT1_PAD_AFTER_SIZE_Y 0
#define INPUT1_PAD_AFTER_SIZE_Z 0
#define INPUT1_PAD_AFTER_SIZE_W 0
#define INPUT1_PAD_AFTER_SIZE_U 0
#define INPUT1_PAD_AFTER_SIZE_V 0
#define INPUT1_PAD_AFTER_FEATURE_NUM 0
#define INPUT1_PAD_AFTER_BATCH_NUM 0
#define INPUT1_X_PITCH 1
#define INPUT1_Y_PITCH 1
#define INPUT1_Z_PITCH (1*1)
#define INPUT1_W_PITCH (1*1*1)
#define INPUT1_U_PITCH (1*1*1*1)
#define INPUT1_V_PITCH (1*1*1*1*1)
#define INPUT1_FEATURE_PITCH (1*1*1*1*1*1)
#define INPUT1_BATCH_PITCH (1*1*1*1*1*1*1)
#define INPUT1_GET_INDEX_SAFE(b, f, y, x) GET_DATA_INDEX_SAFE(INPUT1, b, f, y, x)
#define INPUT1_GET_INDEX(b, f, y, x) GET_DATA_INDEX(INPUT1, b, f, y, x)
#define INPUT1_GET_INDEX_RAW(b, f, y, x) GET_DATA_INDEX_RAW(INPUT1, b, f, y, x)
#define INPUT1_VIEW_OFFSET 0
#define INPUT1_LENGTH 0
#define INPUT1_DIMS 4
#define INPUT1_SIMPLE 1
#define INPUT1_GROUPED 0
#define INPUT1_LAYOUT_BFYX 1
#define INPUT1_TYPE int
#define INPUT1_VAL_MAX INT_MAX
#define INPUT1_VAL_MIN INT_MIN
#define INPUT1_VAL_ONE (int) 1
#define INPUT1_VAL_ZERO (int) 0
#define TO_INPUT1_TYPE(v) convert_int(v)
#define TO_INPUT1_TYPE_SAT(v) convert_int_sat(v)
#define AS_INPUT1_TYPE(v) as_int(v)
#define INPUT1_MAX_FUNC max
#define INPUT1_MIN_FUNC min
#define INPUT1_ABS_FUNC abs
#define INPUT1_TYPE_SIZE 4
#define INPUT1_IS_FP 0
#define INPUT1_OFFSET ((INPUT1_X_PITCH*INPUT1_PAD_BEFORE_SIZE_X) + (INPUT1_Y_PITCH*INPUT1_PAD_BEFORE_SIZE_Y) + (INPUT1_Z_PITCH*INPUT1_PAD_BEFORE_SIZE_Z) + (INPUT1_W_PITCH*INPUT1_PAD_BEFORE_SIZE_W) + (INPUT1_FEATURE_PITCH*INPUT1_PAD_BEFORE_FEATURE_NUM) + (INPUT1_BATCH_PITCH*INPUT1_PAD_BEFORE_BATCH_NUM))
#define INPUT1_PAD_BEFORE (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT1_PAD_AFTER (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT2_SIZE_X 1
#define INPUT2_SIZE_Y 1
#define INPUT2_SIZE_Z 1
#define INPUT2_SIZE_W 1
#define INPUT2_SIZE_U 1
#define INPUT2_SIZE_V 1
#define INPUT2_FEATURE_NUM 1
#define INPUT2_BATCH_NUM (shape_info[16] )
#define INPUT2_PAD_BEFORE_SIZE_X 0
#define INPUT2_PAD_BEFORE_SIZE_Y 0
#define INPUT2_PAD_BEFORE_SIZE_Z 0
#define INPUT2_PAD_BEFORE_SIZE_W 0
#define INPUT2_PAD_BEFORE_SIZE_U 0
#define INPUT2_PAD_BEFORE_SIZE_V 0
#define INPUT2_PAD_BEFORE_FEATURE_NUM 0
#define INPUT2_PAD_BEFORE_BATCH_NUM 0
#define INPUT2_PAD_AFTER_SIZE_X 0
#define INPUT2_PAD_AFTER_SIZE_Y 0
#define INPUT2_PAD_AFTER_SIZE_Z 0
#define INPUT2_PAD_AFTER_SIZE_W 0
#define INPUT2_PAD_AFTER_SIZE_U 0
#define INPUT2_PAD_AFTER_SIZE_V 0
#define INPUT2_PAD_AFTER_FEATURE_NUM 0
#define INPUT2_PAD_AFTER_BATCH_NUM 0
#define INPUT2_X_PITCH 1
#define INPUT2_Y_PITCH 1
#define INPUT2_Z_PITCH (1*1)
#define INPUT2_W_PITCH (1*1*1)
#define INPUT2_U_PITCH (1*1*1*1)
#define INPUT2_V_PITCH (1*1*1*1*1)
#define INPUT2_FEATURE_PITCH (1*1*1*1*1*1)
#define INPUT2_BATCH_PITCH (1*1*1*1*1*1*1)
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
#define OUTPUT_SIZE_X 1
#define OUTPUT_SIZE_Y 1
#define OUTPUT_SIZE_Z 1
#define OUTPUT_SIZE_W 1
#define OUTPUT_SIZE_U 1
#define OUTPUT_SIZE_V 1
#define OUTPUT_FEATURE_NUM 1
#define OUTPUT_BATCH_NUM (shape_info[24] )
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
#define OUTPUT_Y_PITCH 1
#define OUTPUT_Z_PITCH (1*1)
#define OUTPUT_W_PITCH (1*1*1)
#define OUTPUT_U_PITCH (1*1*1*1)
#define OUTPUT_V_PITCH (1*1*1*1*1)
#define OUTPUT_FEATURE_PITCH (1*1*1*1*1*1)
#define OUTPUT_BATCH_PITCH (1*1*1*1*1*1*1)
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


#define GET_UPDATES_INDEX(prefix, idx_order) CAT(prefix, _GET_INDEX)(idx_order)
#define GET_OUTPUT_INDEX(idx_order) OUTPUT_GET_INDEX(idx_order)
#if OUTPUT_DIMS == 4
 #define ORDER b,f,y,x
#elif OUTPUT_DIMS == 5
 #define ORDER b,f,z,y,x
#elif OUTPUT_DIMS == 6
 #define ORDER b,f,w,z,y,x
#endif
#if INPUT2_DIMS == 4
 #define UPD_ORDER upd_b,upd_f,upd_y,upd_x
#elif INPUT2_DIMS == 5
 #define UPD_ORDER upd_b,upd_f,upd_z,upd_y,upd_x
#elif INPUT2_DIMS == 6
 #define UPD_ORDER upd_b,upd_f,upd_w,upd_z,upd_y,upd_x
#endif
#if INPUT1_DIMS == 4
 #define IDX_ORDER idx_b,idx_f,idx_y,idx_x
#elif INPUT1_DIMS == 5
 #define IDX_ORDER idx_b,idx_f,idx_z,idx_y,idx_x
#elif INPUT1_DIMS == 6
 #define IDX_ORDER idx_b,idx_f,idx_w,idx_z,idx_y,idx_x
#endif
#define INDICES_MAX_DIM 6
KERNEL(scatter_nd_update_ref)(OPTIONAL_SHAPE_INFO_ARG
 const __global INPUT0_TYPE* data,
#ifdef IS_SECOND_ITER
 const __global INPUT1_TYPE* indices,
 const __global INPUT2_TYPE* updates,
#endif
 __global OUTPUT_TYPE* output
#if HAS_FUSED_OPS_DECLS
 , FUSED_OPS_DECLS
#endif
)
{
 const uint dim0 = get_global_id(0);
 const uint dim1 = get_global_id(1);
 const uint dim2 = get_global_id(2);
#ifndef IS_SECOND_ITER
 const uint x = dim0 % OUTPUT_SIZE_X;
 const uint y = dim0 / OUTPUT_SIZE_X;
 const uint z = dim1 % OUTPUT_SIZE_Z;
 const uint w = dim1 / OUTPUT_SIZE_Z;
 const uint f = dim2 % OUTPUT_FEATURE_NUM;
 const uint b = dim2 / OUTPUT_FEATURE_NUM;
 const uint input_idx = GET_UPDATES_INDEX(INPUT0, ORDER);
 const uint output_idx = GET_OUTPUT_INDEX(ORDER);
 INPUT0_TYPE val = data[input_idx];
 #if HAS_FUSED_OPS
 FUSED_OPS_FIRST_KERNEL;
 output[output_idx] = TO_OUTPUT_TYPE(FUSED_OPS_RESULT_FIRST_KERNEL);
 #else
 output[output_idx] = ACTIVATION(val, ACTIVATION_PARAMS);
 #endif
#else
 const uint dataND[] = {INPUT0_BLOCK_ND};
 const uint updatesND[] = {INPUT2_BLOCK_ND};
 const uint indicesND[] = {INPUT1_BLOCK_ND};
 const uint size_to_update = dataND[INDICES_LAST_DIM];
 #if INPUT1_DIMS == 4
 const uint indices_dim[INPUT1_DIMS] = {INPUT1_BATCH_NUM, INPUT1_FEATURE_NUM, INPUT1_SIZE_Y, INPUT1_SIZE_X};
 #elif INPUT1_DIMS == 5
 const uint indices_dim[INPUT1_DIMS] = {INPUT1_BATCH_NUM, INPUT1_FEATURE_NUM, INPUT1_SIZE_Z, INPUT1_SIZE_Y, INPUT1_SIZE_X};
 #elif INPUT1_DIMS == 6
 const uint indices_dim[INPUT1_DIMS] = {INPUT1_BATCH_NUM, INPUT1_FEATURE_NUM, INPUT1_SIZE_W, INPUT1_SIZE_Z, INPUT1_SIZE_Y, INPUT1_SIZE_X};
 #endif
 #if INPUT0_DIMS == 4
 const uint data_dim[INPUT0_DIMS] = {INPUT0_BATCH_NUM, INPUT0_FEATURE_NUM, INPUT0_SIZE_Y, INPUT0_SIZE_X};
 #elif INPUT0_DIMS == 5
 const uint data_dim[INPUT0_DIMS] = {INPUT0_BATCH_NUM, INPUT0_FEATURE_NUM, INPUT0_SIZE_Z, INPUT0_SIZE_Y, INPUT0_SIZE_X};
 #elif INPUT0_DIMS == 6
 const uint data_dim[INPUT0_DIMS] = {INPUT0_BATCH_NUM, INPUT0_FEATURE_NUM, INPUT0_SIZE_W, INPUT0_SIZE_Z, INPUT0_SIZE_Y, INPUT0_SIZE_X};
 #endif
 uint idx[INDICES_MAX_DIM] = {0};
 uint rmd_idx = dim2;
 for (int i = 0; i < INDICES_RANK - 1; ++i) {
 idx[i] = rmd_idx / indicesND[1 + i];
 rmd_idx %= indicesND[1 + i];
 }
 uint out[INDICES_MAX_DIM] = {0};
 for (int i = 0; i < indices_dim[INDICES_RANK - 1]; ++i) {
 idx[INDICES_RANK - 1] = i;
 const uint idx_b = idx[0];
 const uint idx_f = idx[1];
 #if INPUT1_DIMS == 4
 const uint idx_y = idx[2];
 const uint idx_x = idx[3];
 #elif INPUT1_DIMS == 5
 const uint idx_z = idx[2];
 const uint idx_y = idx[3];
 const uint idx_x = idx[4];
 #elif INPUT1_DIMS == 6
 const uint idx_w = idx[2];
 const uint idx_z = idx[3];
 const uint idx_y = idx[4];
 const uint idx_x = idx[5];
 #endif
 uint index = GET_UPDATES_INDEX(INPUT1, IDX_ORDER);
 out[i] = indices[index];
 if(out[i] >= data_dim[i])
 out[i] = data_dim[i] - 1;
 }
 for (int i = 0; i < size_to_update; ++i) {
 uint upd[INDICES_MAX_DIM] = {0};
 for (int j = 0; j < INDICES_RANK - 1; ++j) {
 upd[j] = idx[j];
 }
 uint data_rmd = i, updates_rmd = i;
 for (int j = indices_dim[INDICES_RANK - 1]; j < INPUT0_DIMS; ++j) {
 out[j] = data_rmd / dataND[j + 1];
 data_rmd %= dataND[j + 1];
 }
 for (int k = INDICES_RANK - 1; k < INPUT2_DIMS; ++k) {
 upd[k] = updates_rmd / updatesND[k + 1];
 updates_rmd %= updatesND[k + 1];
 }
 const uint upd_b = upd[0];
 const uint upd_f = upd[1];
 #if INPUT2_DIMS == 4
 const uint upd_y = upd[2];
 const uint upd_x = upd[3];
 #elif INPUT2_DIMS == 5
 const uint upd_z = upd[2];
 const uint upd_y = upd[3];
 const uint upd_x = upd[4];
 #elif INPUT2_DIMS == 6
 const uint upd_w = upd[2];
 const uint upd_z = upd[3];
 const uint upd_y = upd[4];
 const uint upd_x = upd[5];
 #endif
 uint upd_idx = GET_UPDATES_INDEX(INPUT2, UPD_ORDER);
 const uint b = out[0];
 const uint f = out[1];
 #if INPUT0_DIMS == 4
 const uint y = out[2];
 const uint x = out[3];
 #elif INPUT0_DIMS == 5
 const uint z = out[2];
 const uint y = out[3];
 const uint x = out[4];
 #elif INPUT0_DIMS == 6
 const uint w = out[2];
 const uint z = out[3];
 const uint y = out[4];
 const uint x = out[5];
 #endif
 uint out_idx = GET_OUTPUT_INDEX(ORDER);
 INPUT2_TYPE val = updates[upd_idx];
 #if HAS_FUSED_OPS
 FUSED_OPS_SECOND_KERNEL;
 output[out_idx] = TO_OUTPUT_TYPE(FUSED_OPS_RESULT_SECOND_KERNEL);
 #else
 output[out_idx] = ACTIVATION(val, ACTIVATION_PARAMS);
 #endif
 }
#endif
}
#ifdef GET_UPDATES_INDEX
#undef GET_UPDATES_INDEX
#endif
#ifdef GET_OUTPUT_INDEX
#undef GET_OUTPUT_INDEX
#endif
#ifdef ORDER
#undef ORDER
#endif
#ifdef UPD_ORDER
#undef UPD_ORDER
#endif
#ifdef IDX_ORDER
#undef IDX_ORDER
#endif
#ifdef INDICES_MAX_DIM
#undef INDICES_MAX_DIM
#endif
#ifdef GET_UPDATES_INDEX
#undef GET_UPDATES_INDEX
#endif
#ifdef GET_OUTPUT_INDEX
#undef GET_OUTPUT_INDEX
#endif
#ifdef ORDER
#undef ORDER
#endif
#ifdef ORDER
#undef ORDER
#endif
#ifdef ORDER
#undef ORDER
#endif
#ifdef UPD_ORDER
#undef UPD_ORDER
#endif
#ifdef UPD_ORDER
#undef UPD_ORDER
#endif
#ifdef UPD_ORDER
#undef UPD_ORDER
#endif
#ifdef IDX_ORDER
#undef IDX_ORDER
#endif
#ifdef IDX_ORDER
#undef IDX_ORDER
#endif
#ifdef IDX_ORDER
#undef IDX_ORDER
#endif
#ifdef INDICES_MAX_DIM
#undef INDICES_MAX_DIM
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

//====================================================
// Kernel template: scatter_nd_update_ref 
// Kernel name: scatter_nd_update_ref_13885323453389162897_1_0__sa
#define KERNEL(name) __kernel void scatter_nd_update_ref_13885323453389162897_1_0__sa
#define KERNEL_ID scatter_nd_update_ref_13885323453389162897_1_0__sa
#define FUNC(name)  _##name##_scatter_nd_update_ref_13885323453389162897_1_0__sa
#define FUNC_CALL(name)  _##name##_scatter_nd_update_ref_13885323453389162897_1_0__sa
#define CONST_ARRAY_DECL(name) __constant size_t  _##name##_scatter_nd_update_ref_13885323453389162897_1_0__sa []
#define CONST_ARRAY_REF(name)  _##name##_scatter_nd_update_ref_13885323453389162897_1_0__sa
#define FP64_SUPPORTED 0
#define FP16_SUPPORTED 1
#define FP16_UNIT_USED 1
#define INT8_UNIT_USED 0
#define INT32_UNIT_USED 1
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
#define INPUT0_SIZE_X 1
#define INPUT0_SIZE_Y 1
#define INPUT0_SIZE_Z 1
#define INPUT0_SIZE_W 1
#define INPUT0_SIZE_U 1
#define INPUT0_SIZE_V 1
#define INPUT0_FEATURE_NUM 1
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
#define INPUT0_Y_PITCH 1
#define INPUT0_Z_PITCH (1*1)
#define INPUT0_W_PITCH (1*1*1)
#define INPUT0_U_PITCH (1*1*1*1)
#define INPUT0_V_PITCH (1*1*1*1*1)
#define INPUT0_FEATURE_PITCH (1*1*1*1*1*1)
#define INPUT0_BATCH_PITCH (1*1*1*1*1*1*1)
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
#define INPUT1_SIZE_X 1
#define INPUT1_SIZE_Y 1
#define INPUT1_SIZE_Z 1
#define INPUT1_SIZE_W 1
#define INPUT1_SIZE_U 1
#define INPUT1_SIZE_V 1
#define INPUT1_FEATURE_NUM 1
#define INPUT1_BATCH_NUM (shape_info[8] )
#define INPUT1_PAD_BEFORE_SIZE_X 0
#define INPUT1_PAD_BEFORE_SIZE_Y 0
#define INPUT1_PAD_BEFORE_SIZE_Z 0
#define INPUT1_PAD_BEFORE_SIZE_W 0
#define INPUT1_PAD_BEFORE_SIZE_U 0
#define INPUT1_PAD_BEFORE_SIZE_V 0
#define INPUT1_PAD_BEFORE_FEATURE_NUM 0
#define INPUT1_PAD_BEFORE_BATCH_NUM 0
#define INPUT1_PAD_AFTER_SIZE_X 0
#define INPUT1_PAD_AFTER_SIZE_Y 0
#define INPUT1_PAD_AFTER_SIZE_Z 0
#define INPUT1_PAD_AFTER_SIZE_W 0
#define INPUT1_PAD_AFTER_SIZE_U 0
#define INPUT1_PAD_AFTER_SIZE_V 0
#define INPUT1_PAD_AFTER_FEATURE_NUM 0
#define INPUT1_PAD_AFTER_BATCH_NUM 0
#define INPUT1_X_PITCH 1
#define INPUT1_Y_PITCH 1
#define INPUT1_Z_PITCH (1*1)
#define INPUT1_W_PITCH (1*1*1)
#define INPUT1_U_PITCH (1*1*1*1)
#define INPUT1_V_PITCH (1*1*1*1*1)
#define INPUT1_FEATURE_PITCH (1*1*1*1*1*1)
#define INPUT1_BATCH_PITCH (1*1*1*1*1*1*1)
#define INPUT1_GET_INDEX_SAFE(b, f, y, x) GET_DATA_INDEX_SAFE(INPUT1, b, f, y, x)
#define INPUT1_GET_INDEX(b, f, y, x) GET_DATA_INDEX(INPUT1, b, f, y, x)
#define INPUT1_GET_INDEX_RAW(b, f, y, x) GET_DATA_INDEX_RAW(INPUT1, b, f, y, x)
#define INPUT1_VIEW_OFFSET 0
#define INPUT1_LENGTH 0
#define INPUT1_DIMS 4
#define INPUT1_SIMPLE 1
#define INPUT1_GROUPED 0
#define INPUT1_LAYOUT_BFYX 1
#define INPUT1_TYPE int
#define INPUT1_VAL_MAX INT_MAX
#define INPUT1_VAL_MIN INT_MIN
#define INPUT1_VAL_ONE (int) 1
#define INPUT1_VAL_ZERO (int) 0
#define TO_INPUT1_TYPE(v) convert_int(v)
#define TO_INPUT1_TYPE_SAT(v) convert_int_sat(v)
#define AS_INPUT1_TYPE(v) as_int(v)
#define INPUT1_MAX_FUNC max
#define INPUT1_MIN_FUNC min
#define INPUT1_ABS_FUNC abs
#define INPUT1_TYPE_SIZE 4
#define INPUT1_IS_FP 0
#define INPUT1_OFFSET ((INPUT1_X_PITCH*INPUT1_PAD_BEFORE_SIZE_X) + (INPUT1_Y_PITCH*INPUT1_PAD_BEFORE_SIZE_Y) + (INPUT1_Z_PITCH*INPUT1_PAD_BEFORE_SIZE_Z) + (INPUT1_W_PITCH*INPUT1_PAD_BEFORE_SIZE_W) + (INPUT1_FEATURE_PITCH*INPUT1_PAD_BEFORE_FEATURE_NUM) + (INPUT1_BATCH_PITCH*INPUT1_PAD_BEFORE_BATCH_NUM))
#define INPUT1_PAD_BEFORE (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT1_PAD_AFTER (size_t []){ 0,0,0,0,0,0,0,0,0, } 
#define INPUT2_SIZE_X 1
#define INPUT2_SIZE_Y 1
#define INPUT2_SIZE_Z 1
#define INPUT2_SIZE_W 1
#define INPUT2_SIZE_U 1
#define INPUT2_SIZE_V 1
#define INPUT2_FEATURE_NUM 1
#define INPUT2_BATCH_NUM (shape_info[16] )
#define INPUT2_PAD_BEFORE_SIZE_X 0
#define INPUT2_PAD_BEFORE_SIZE_Y 0
#define INPUT2_PAD_BEFORE_SIZE_Z 0
#define INPUT2_PAD_BEFORE_SIZE_W 0
#define INPUT2_PAD_BEFORE_SIZE_U 0
#define INPUT2_PAD_BEFORE_SIZE_V 0
#define INPUT2_PAD_BEFORE_FEATURE_NUM 0
#define INPUT2_PAD_BEFORE_BATCH_NUM 0
#define INPUT2_PAD_AFTER_SIZE_X 0
#define INPUT2_PAD_AFTER_SIZE_Y 0
#define INPUT2_PAD_AFTER_SIZE_Z 0
#define INPUT2_PAD_AFTER_SIZE_W 0
#define INPUT2_PAD_AFTER_SIZE_U 0
#define INPUT2_PAD_AFTER_SIZE_V 0
#define INPUT2_PAD_AFTER_FEATURE_NUM 0
#define INPUT2_PAD_AFTER_BATCH_NUM 0
#define INPUT2_X_PITCH 1
#define INPUT2_Y_PITCH 1
#define INPUT2_Z_PITCH (1*1)
#define INPUT2_W_PITCH (1*1*1)
#define INPUT2_U_PITCH (1*1*1*1)
#define INPUT2_V_PITCH (1*1*1*1*1)
#define INPUT2_FEATURE_PITCH (1*1*1*1*1*1)
#define INPUT2_BATCH_PITCH (1*1*1*1*1*1*1)
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
#define OUTPUT_SIZE_X 1
#define OUTPUT_SIZE_Y 1
#define OUTPUT_SIZE_Z 1
#define OUTPUT_SIZE_W 1
#define OUTPUT_SIZE_U 1
#define OUTPUT_SIZE_V 1
#define OUTPUT_FEATURE_NUM 1
#define OUTPUT_BATCH_NUM (shape_info[24] )
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
#define OUTPUT_Y_PITCH 1
#define OUTPUT_Z_PITCH (1*1)
#define OUTPUT_W_PITCH (1*1*1)
#define OUTPUT_U_PITCH (1*1*1*1)
#define OUTPUT_V_PITCH (1*1*1*1*1)
#define OUTPUT_FEATURE_PITCH (1*1*1*1*1*1)
#define OUTPUT_BATCH_PITCH (1*1*1*1*1*1*1)
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
#define IS_SECOND_ITER true
#define INPUT0_BLOCK_ND ((shape_info[0] )*(1*(1*(1*1)))),(1*(1*(1*1))),(1*(1*1)),(1*1),1,
#define INPUT1_BLOCK_ND ((shape_info[8] )*1),1,
#define INPUT2_BLOCK_ND ((shape_info[16] )*(1*(1*(1*1)))),(1*(1*(1*1))),(1*(1*1)),(1*1),1,
#define INDICES_RANK 2
#define INDICES_LAST_DIM 1


#define GET_UPDATES_INDEX(prefix, idx_order) CAT(prefix, _GET_INDEX)(idx_order)
#define GET_OUTPUT_INDEX(idx_order) OUTPUT_GET_INDEX(idx_order)
#if OUTPUT_DIMS == 4
 #define ORDER b,f,y,x
#elif OUTPUT_DIMS == 5
 #define ORDER b,f,z,y,x
#elif OUTPUT_DIMS == 6
 #define ORDER b,f,w,z,y,x
#endif
#if INPUT2_DIMS == 4
 #define UPD_ORDER upd_b,upd_f,upd_y,upd_x
#elif INPUT2_DIMS == 5
 #define UPD_ORDER upd_b,upd_f,upd_z,upd_y,upd_x
#elif INPUT2_DIMS == 6
 #define UPD_ORDER upd_b,upd_f,upd_w,upd_z,upd_y,upd_x
#endif
#if INPUT1_DIMS == 4
 #define IDX_ORDER idx_b,idx_f,idx_y,idx_x
#elif INPUT1_DIMS == 5
 #define IDX_ORDER idx_b,idx_f,idx_z,idx_y,idx_x
#elif INPUT1_DIMS == 6
 #define IDX_ORDER idx_b,idx_f,idx_w,idx_z,idx_y,idx_x
#endif
#define INDICES_MAX_DIM 6
KERNEL(scatter_nd_update_ref)(OPTIONAL_SHAPE_INFO_ARG
 const __global INPUT0_TYPE* data,
#ifdef IS_SECOND_ITER
 const __global INPUT1_TYPE* indices,
 const __global INPUT2_TYPE* updates,
#endif
 __global OUTPUT_TYPE* output
#if HAS_FUSED_OPS_DECLS
 , FUSED_OPS_DECLS
#endif
)
{
 const uint dim0 = get_global_id(0);
 const uint dim1 = get_global_id(1);
 const uint dim2 = get_global_id(2);
#ifndef IS_SECOND_ITER
 const uint x = dim0 % OUTPUT_SIZE_X;
 const uint y = dim0 / OUTPUT_SIZE_X;
 const uint z = dim1 % OUTPUT_SIZE_Z;
 const uint w = dim1 / OUTPUT_SIZE_Z;
 const uint f = dim2 % OUTPUT_FEATURE_NUM;
 const uint b = dim2 / OUTPUT_FEATURE_NUM;
 const uint input_idx = GET_UPDATES_INDEX(INPUT0, ORDER);
 const uint output_idx = GET_OUTPUT_INDEX(ORDER);
 INPUT0_TYPE val = data[input_idx];
 #if HAS_FUSED_OPS
 FUSED_OPS_FIRST_KERNEL;
 output[output_idx] = TO_OUTPUT_TYPE(FUSED_OPS_RESULT_FIRST_KERNEL);
 #else
 output[output_idx] = ACTIVATION(val, ACTIVATION_PARAMS);
 #endif
#else
 const uint dataND[] = {INPUT0_BLOCK_ND};
 const uint updatesND[] = {INPUT2_BLOCK_ND};
 const uint indicesND[] = {INPUT1_BLOCK_ND};
 const uint size_to_update = dataND[INDICES_LAST_DIM];
 #if INPUT1_DIMS == 4
 const uint indices_dim[INPUT1_DIMS] = {INPUT1_BATCH_NUM, INPUT1_FEATURE_NUM, INPUT1_SIZE_Y, INPUT1_SIZE_X};
 #elif INPUT1_DIMS == 5
 const uint indices_dim[INPUT1_DIMS] = {INPUT1_BATCH_NUM, INPUT1_FEATURE_NUM, INPUT1_SIZE_Z, INPUT1_SIZE_Y, INPUT1_SIZE_X};
 #elif INPUT1_DIMS == 6
 const uint indices_dim[INPUT1_DIMS] = {INPUT1_BATCH_NUM, INPUT1_FEATURE_NUM, INPUT1_SIZE_W, INPUT1_SIZE_Z, INPUT1_SIZE_Y, INPUT1_SIZE_X};
 #endif
 #if INPUT0_DIMS == 4
 const uint data_dim[INPUT0_DIMS] = {INPUT0_BATCH_NUM, INPUT0_FEATURE_NUM, INPUT0_SIZE_Y, INPUT0_SIZE_X};
 #elif INPUT0_DIMS == 5
 const uint data_dim[INPUT0_DIMS] = {INPUT0_BATCH_NUM, INPUT0_FEATURE_NUM, INPUT0_SIZE_Z, INPUT0_SIZE_Y, INPUT0_SIZE_X};
 #elif INPUT0_DIMS == 6
 const uint data_dim[INPUT0_DIMS] = {INPUT0_BATCH_NUM, INPUT0_FEATURE_NUM, INPUT0_SIZE_W, INPUT0_SIZE_Z, INPUT0_SIZE_Y, INPUT0_SIZE_X};
 #endif
 uint idx[INDICES_MAX_DIM] = {0};
 uint rmd_idx = dim2;
 for (int i = 0; i < INDICES_RANK - 1; ++i) {
 idx[i] = rmd_idx / indicesND[1 + i];
 rmd_idx %= indicesND[1 + i];
 }
 uint out[INDICES_MAX_DIM] = {0};
 for (int i = 0; i < indices_dim[INDICES_RANK - 1]; ++i) {
 idx[INDICES_RANK - 1] = i;
 const uint idx_b = idx[0];
 const uint idx_f = idx[1];
 #if INPUT1_DIMS == 4
 const uint idx_y = idx[2];
 const uint idx_x = idx[3];
 #elif INPUT1_DIMS == 5
 const uint idx_z = idx[2];
 const uint idx_y = idx[3];
 const uint idx_x = idx[4];
 #elif INPUT1_DIMS == 6
 const uint idx_w = idx[2];
 const uint idx_z = idx[3];
 const uint idx_y = idx[4];
 const uint idx_x = idx[5];
 #endif
 uint index = GET_UPDATES_INDEX(INPUT1, IDX_ORDER);
 out[i] = indices[index];
 if(out[i] >= data_dim[i])
 out[i] = data_dim[i] - 1;
 }
 for (int i = 0; i < size_to_update; ++i) {
 uint upd[INDICES_MAX_DIM] = {0};
 for (int j = 0; j < INDICES_RANK - 1; ++j) {
 upd[j] = idx[j];
 }
 uint data_rmd = i, updates_rmd = i;
 for (int j = indices_dim[INDICES_RANK - 1]; j < INPUT0_DIMS; ++j) {
 out[j] = data_rmd / dataND[j + 1];
 data_rmd %= dataND[j + 1];
 }
 for (int k = INDICES_RANK - 1; k < INPUT2_DIMS; ++k) {
 upd[k] = updates_rmd / updatesND[k + 1];
 updates_rmd %= updatesND[k + 1];
 }
 const uint upd_b = upd[0];
 const uint upd_f = upd[1];
 #if INPUT2_DIMS == 4
 const uint upd_y = upd[2];
 const uint upd_x = upd[3];
 #elif INPUT2_DIMS == 5
 const uint upd_z = upd[2];
 const uint upd_y = upd[3];
 const uint upd_x = upd[4];
 #elif INPUT2_DIMS == 6
 const uint upd_w = upd[2];
 const uint upd_z = upd[3];
 const uint upd_y = upd[4];
 const uint upd_x = upd[5];
 #endif
 uint upd_idx = GET_UPDATES_INDEX(INPUT2, UPD_ORDER);
 const uint b = out[0];
 const uint f = out[1];
 #if INPUT0_DIMS == 4
 const uint y = out[2];
 const uint x = out[3];
 #elif INPUT0_DIMS == 5
 const uint z = out[2];
 const uint y = out[3];
 const uint x = out[4];
 #elif INPUT0_DIMS == 6
 const uint w = out[2];
 const uint z = out[3];
 const uint y = out[4];
 const uint x = out[5];
 #endif
 uint out_idx = GET_OUTPUT_INDEX(ORDER);
 INPUT2_TYPE val = updates[upd_idx];
 #if HAS_FUSED_OPS
 FUSED_OPS_SECOND_KERNEL;
 output[out_idx] = TO_OUTPUT_TYPE(FUSED_OPS_RESULT_SECOND_KERNEL);
 #else
 output[out_idx] = ACTIVATION(val, ACTIVATION_PARAMS);
 #endif
 }
#endif
}
#ifdef GET_UPDATES_INDEX
#undef GET_UPDATES_INDEX
#endif
#ifdef GET_OUTPUT_INDEX
#undef GET_OUTPUT_INDEX
#endif
#ifdef ORDER
#undef ORDER
#endif
#ifdef UPD_ORDER
#undef UPD_ORDER
#endif
#ifdef IDX_ORDER
#undef IDX_ORDER
#endif
#ifdef INDICES_MAX_DIM
#undef INDICES_MAX_DIM
#endif
#ifdef GET_UPDATES_INDEX
#undef GET_UPDATES_INDEX
#endif
#ifdef GET_OUTPUT_INDEX
#undef GET_OUTPUT_INDEX
#endif
#ifdef ORDER
#undef ORDER
#endif
#ifdef ORDER
#undef ORDER
#endif
#ifdef ORDER
#undef ORDER
#endif
#ifdef UPD_ORDER
#undef UPD_ORDER
#endif
#ifdef UPD_ORDER
#undef UPD_ORDER
#endif
#ifdef UPD_ORDER
#undef UPD_ORDER
#endif
#ifdef IDX_ORDER
#undef IDX_ORDER
#endif
#ifdef IDX_ORDER
#undef IDX_ORDER
#endif
#ifdef IDX_ORDER
#undef IDX_ORDER
#endif
#ifdef INDICES_MAX_DIM
#undef INDICES_MAX_DIM
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
#ifdef IS_SECOND_ITER
#undef IS_SECOND_ITER
#endif
#ifdef INPUT0_BLOCK_ND
#undef INPUT0_BLOCK_ND
#endif
#ifdef INPUT1_BLOCK_ND
#undef INPUT1_BLOCK_ND
#endif
#ifdef INPUT2_BLOCK_ND
#undef INPUT2_BLOCK_ND
#endif
#ifdef INDICES_RANK
#undef INDICES_RANK
#endif
#ifdef INDICES_LAST_DIM
#undef INDICES_LAST_DIM
#endif

//====================================================
// Kernel template: sdpa_opt_single_token 
// Kernel name: sdpa_opt_single_token_8660372428234100028_0_0__sa
#define KERNEL(name) __kernel void sdpa_opt_single_token_8660372428234100028_0_0__sa
#define KERNEL_ID sdpa_opt_single_token_8660372428234100028_0_0__sa
#define FUNC(name)  _##name##_sdpa_opt_single_token_8660372428234100028_0_0__sa
#define FUNC_CALL(name)  _##name##_sdpa_opt_single_token_8660372428234100028_0_0__sa
#define CONST_ARRAY_DECL(name) __constant size_t  _##name##_sdpa_opt_single_token_8660372428234100028_0_0__sa []
#define CONST_ARRAY_REF(name)  _##name##_sdpa_opt_single_token_8660372428234100028_0_0__sa
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
#define BROADCAST_GROUP_SIZE 7
#define DO_BROADCAST_KEY_VALUE f /= 7;
#define IS_CAUSAL 0
#define HAS_ATTN_MASK_INPUT 1
#define HAS_SCALE_INPUT 0
#define IS_KV_COMPRESSED 0
#define INPUT0_DIMS_ORDER b,f,w,z,y,x
#define INPUT1_DIMS_ORDER b,f,w,z,y,x
#define INPUT2_DIMS_ORDER b,f,w,z,y,x
#define TARGET_SEQ_LEN (shape_info[6] )
#define NUM_HEADS 28
#define NUM_KV_HEADS -1
#define SOURCE_SEQ_LEN (shape_info[14] )
#define SOFTMAX_ACCUMULATOR_TYPE float
#define SOFTMAX_ACCUMULATOR_VAL_MAX FLT_MAX
#define SOFTMAX_ACCUMULATOR_VAL_MIN -SOFTMAX_ACCUMULATOR_VAL_MAX
#define SOFTMAX_ACCUMULATOR_VAL_ONE 1.0f
#define SOFTMAX_ACCUMULATOR_VAL_ZERO 0.0f
#define TO_SOFTMAX_ACCUMULATOR_TYPE(v) convert_float(v)
#define TO_SOFTMAX_ACCUMULATOR_TYPE_SAT(v) convert_float(v)
#define AS_SOFTMAX_ACCUMULATOR_TYPE(v) as_float(v)
#define SOFTMAX_ACCUMULATOR_MAX_FUNC fmax
#define SOFTMAX_ACCUMULATOR_MIN_FUNC fmin
#define SOFTMAX_ACCUMULATOR_ABS_FUNC fabs
#define SOFTMAX_ACCUMULATOR_TYPE_SIZE 4
#define SOFTMAX_ACCUMULATOR_IS_FP 1
#define SUBGROUP_SIZE 16
#define HEAD_SIZE 128
#define SEQ_LEN_PARTITION_SIZE 256
#define TARGET_SEQ_LEN_BLOCK_SIZE 1
#define SDPA_STAGE_0 1
#define SG_SCALE_FACTOR 2
#define STATIC_SCALE_VALUE_INV as_float(0x413504f3)/*1.131371e+01*/
#define STATIC_SCALE_VALUE as_float(0x3db504f2)/*8.838834e-02*/


inline uint FUNC(get_input0_index_nt)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#if INPUT0_SIMPLE
 return GET_DATA_INDEX_6D(INPUT0, b, f, w, z, y, x);
#else
#if INPUT0_DIMS == 4
 return INPUT0_GET_INDEX(b, f, y, x);
#elif INPUT0_DIMS == 5
 return INPUT0_GET_INDEX(b, f, z, y, x);
#elif INPUT0_DIMS == 6
 return INPUT0_GET_INDEX(b, f, w, z, y, x);
#else
# error sdpa_opt.cl : Unsupported input 0 format
#endif
#endif
}
inline uint FUNC(get_input0_index)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef INPUT0_DIMS_ORDER
 return FUNC_CALL(get_input0_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT0_DIMS_ORDER);
#else
 return FUNC_CALL(get_input0_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR b, f, w, z, y, x);
#endif
}
inline uint FUNC(get_input1_index_nt)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef DO_BROADCAST_KEY_VALUE
 DO_BROADCAST_KEY_VALUE;
#endif
#if INPUT1_SIMPLE
 return GET_DATA_INDEX_6D(INPUT1, b, f, w, z, y, x);
#else
#if INPUT1_DIMS == 4
 return INPUT1_GET_INDEX(b, f, y, x);
#elif INPUT1_DIMS == 5
 return INPUT1_GET_INDEX(b, f, z, y, x);
#elif INPUT1_DIMS == 6
 return INPUT1_GET_INDEX(b, f, w, z, y, x);
#else
# error sdpa_opt.cl : Unsupported input 1 format
#endif
#endif
}
inline uint FUNC(get_input1_index)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef INPUT1_DIMS_ORDER
 return FUNC_CALL(get_input1_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT1_DIMS_ORDER);
#else
 return FUNC_CALL(get_input1_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR b, f, w, z, y, x);
#endif
}
inline uint FUNC(get_input2_index_nt)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef DO_BROADCAST_KEY_VALUE
 DO_BROADCAST_KEY_VALUE;
#endif
#if INPUT2_SIMPLE
 return GET_DATA_INDEX_6D_SAFE(INPUT2, b, f, w, z, y, x);
#else
#if INPUT2_DIMS == 4
 return INPUT2_GET_INDEX(b, f, y, x);
#elif INPUT2_DIMS == 5
 return INPUT2_GET_INDEX(b, f, z, y, x);
#elif INPUT2_DIMS == 6
 return INPUT2_GET_INDEX(b, f, w, z, y, x);
#else
# error sdpa_opt.cl : Unsupported input 1 format
#endif
#endif
}
inline uint FUNC(get_input2_index)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef INPUT2_DIMS_ORDER
 return FUNC_CALL(get_input2_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT2_DIMS_ORDER);
#else
 return FUNC_CALL(get_input2_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR b, f, w, z, y, x);
#endif
}
#ifdef BEAM_TABLE_TYPE
inline uint FUNC(get_bt_index_nt)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#if BEAM_TABLE_SIMPLE
 return GET_DATA_INDEX_6D_SAFE(BEAM_TABLE, b, f, w, z, y, x);
#else
# error sdpa_opt.cl : Unsupported beam table format
#endif
}
inline uint FUNC(get_bt_index_key)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
 return FUNC_CALL(get_bt_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT1_DIMS_ORDER);
}
inline uint FUNC(get_bt_index_value)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
 return FUNC_CALL(get_bt_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT2_DIMS_ORDER);
}
#endif
#define OUTPUT_BLOCK_READ(ptr, offset) BLOCK_READN(OUTPUT_TYPE, 1, ptr, offset)
#define OUTPUT_BLOCK_WRITE(ptr, offset, val) BLOCK_WRITEN(OUTPUT_TYPE, 1, ptr, offset, val)
#define VALUE_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT2_TYPE, 1, ptr, offset)
#define SUBGROUPS_PER_WG (HEAD_SIZE * SG_SCALE_FACTOR / SUBGROUP_SIZE)
#if IS_KV_COMPRESSED
#if COMPRESSED_PER_HEAD
 #define GET_COMPRESSION_INDEX(INPUT, b, f, y, x) GET_DATA_INDEX(INPUT, (b), (f), (y), (0));
#else
 #define GET_COMPRESSION_INDEX(INPUT, b, f, y, x) GET_DATA_INDEX(INPUT, (b), (0), (y), (0));
#endif
#endif
#ifdef SDPA_STAGE_0
#if TARGET_SEQ_LEN_BLOCK_SIZE == 1
REQD_SUB_GROUP_SIZE(SUBGROUP_SIZE)
__attribute__((reqd_work_group_size(1, 1, HEAD_SIZE * SG_SCALE_FACTOR)))
KERNEL(sdpa_opt)(
 OPTIONAL_SHAPE_INFO_ARG
 const __global INPUT0_TYPE* query_input,
 const __global INPUT1_TYPE* key_input,
 const __global INPUT2_TYPE* value_input,
#if HAS_ATTN_MASK_INPUT
 const __global INPUT3_TYPE* attn_mask,
#endif
#if HAS_SCALE_INPUT
 const __global INPUT4_TYPE* scale,
#endif
 __global OUTPUT_TYPE* output,
#if IS_KV_COMPRESSED
 const __global KEY_COMPRESSION_SCALE_TYPE* key_scale,
 const __global VALUE_COMPRESSION_SCALE_TYPE* val_scale,
#endif
#ifdef BEAM_TABLE_TYPE
 const __global BEAM_TABLE_TYPE* beam_table,
#endif
 __global SOFTMAX_ACCUMULATOR_TYPE* exp_sums,
 __global SOFTMAX_ACCUMULATOR_TYPE* max_logits,
 __global OUTPUT_TYPE* tmp_out
)
{
 const uint batch_idx = get_global_id(0);
 const uint b0_idx = batch_idx / NUM_HEADS;
 const uint b1_idx = batch_idx % NUM_HEADS;
 const uint target_seq_idx = get_global_id(1);
 const uint lid = get_local_id(2);
#if SG_SCALE_FACTOR == 2
 const uint head_size_idx = lid % HEAD_SIZE;
#elif SG_SCALE_FACTOR == 1
 const uint head_size_idx = lid;
#else
 #error "sdpa_opt.cl: Unsupported scale factor"
#endif
#if SUBGROUPS_PER_WG > SUBGROUP_SIZE
 #error "sdpa_opt.cl: Number of subgroups per work group should be less than subgroup_size
#endif
 const uint sgid = get_sub_group_id();
 const uint sglid = get_sub_group_local_id();
 const uint partition_idx = get_group_id(2);
 const uint num_of_partitions = get_num_groups(2);
 const uint wi_num_per_partition = get_local_size(2);
 const uint start_partition_idx = partition_idx * SEQ_LEN_PARTITION_SIZE;
 const uint partition_seq_len =
 ((partition_idx + 1) < num_of_partitions) ? (SEQ_LEN_PARTITION_SIZE)
 : (SOURCE_SEQ_LEN - partition_idx * SEQ_LEN_PARTITION_SIZE);
 __local INPUT0_TYPE query_local[HEAD_SIZE * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE qk_local[SEQ_LEN_PARTITION_SIZE * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE qk_max_vals[SUBGROUPS_PER_WG * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE qk_sum_vals[SUBGROUPS_PER_WG * TARGET_SEQ_LEN_BLOCK_SIZE];
 {
 SOFTMAX_ACCUMULATOR_TYPE qk_max[TARGET_SEQ_LEN_BLOCK_SIZE] = {SOFTMAX_ACCUMULATOR_VAL_MIN};
 for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
 qk_max[i] = SOFTMAX_ACCUMULATOR_VAL_MIN;
 }
 {
#if HAS_SCALE_INPUT
 const OUTPUT_TYPE scale_val = *scale;
#else
 const OUTPUT_TYPE scale_val = OUTPUT_VAL_ONE / sqrt(TO_OUTPUT_TYPE(HEAD_SIZE));
#endif
 {
 #define QUERY_STEP_LOCAL SUBGROUP_SIZE * SUBGROUPS_PER_WG
 uint query_local_offset = sgid * SUBGROUP_SIZE + sglid;
 const uint seq_idx_end = 1;
#ifdef INPUT0_DIMS_ORDER
 uint query_offset = FUNC_CALL(get_input0_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, target_seq_idx, (sgid * SUBGROUP_SIZE));
 uint query_offset_next_seq = FUNC_CALL(get_input0_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, target_seq_idx + 1, (sgid * SUBGROUP_SIZE));
 const uint query_pitch = query_offset_next_seq - query_offset;
#else
 uint query_offset = INPUT0_GET_INDEX(b0_idx, b1_idx, target_seq_idx, (sgid * SUBGROUP_SIZE));
 const uint query_pitch = QUERY_STEP_LOCAL;
#endif
#if SG_SCALE_FACTOR == 2
 if (sgid < HEAD_SIZE / SUBGROUP_SIZE) {
#else
 {
#endif
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 #define QUERY_BLOCK_SIZE 1
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, QUERY_BLOCK_SIZE, query_input, query_offset);
 query_local[query_local_offset] = val * scale_val;
 query_local_offset += QUERY_STEP_LOCAL;
 query_offset += query_pitch;
 }
 }
 #undef QUERY_BLOCK_SIZE
 #undef QUERY_STEP
 barrier(CLK_LOCAL_MEM_FENCE);
 }
 for (uint seq_len = sgid; seq_len < partition_seq_len; seq_len += (HEAD_SIZE / SUBGROUP_SIZE) * SG_SCALE_FACTOR) {
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_key)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + seq_len, 0)];
#else
 const uint b_idx = b0_idx;
#endif
#ifdef INPUT1_DIMS_ORDER
 uint key_offset = FUNC_CALL(get_input1_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + seq_len, 0);
#else
 uint key_offset = INPUT1_GET_INDEX(b_idx, b1_idx, start_partition_idx + seq_len, 0);
#endif
 SOFTMAX_ACCUMULATOR_TYPE acc[TARGET_SEQ_LEN_BLOCK_SIZE] = {SOFTMAX_ACCUMULATOR_VAL_ZERO};
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(KEY_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + seq_len, 0);
 KEY_COMPRESSION_SCALE_TYPE comp_scale = key_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 KEY_COMPRESSION_SCALE_TYPE comp_zp = key_scale[comp_offset + 1];
#endif
#endif
 uint head_idx_index = 0;
 #define KEY_BLOCK_SIZE 8
 for (; head_idx_index + (KEY_BLOCK_SIZE * SUBGROUP_SIZE) <= HEAD_SIZE; head_idx_index += SUBGROUP_SIZE * KEY_BLOCK_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, KEY_BLOCK_SIZE, ptr, offset);
 #define KEY_BLOCK MAKE_VECTOR_TYPE(INPUT1_TYPE, KEY_BLOCK_SIZE)
 #define KEY_BLOCK_UNCOMPRESSED MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, KEY_BLOCK_SIZE)
 #define TO_KEY_BLOCK_UNCOMPRESSED_TYPE(val) CAT(convert_, KEY_BLOCK_UNCOMPRESSED)(val)
 #define QUERY_BLOCK MAKE_VECTOR_TYPE(INPUT0_TYPE, KEY_BLOCK_SIZE)
 KEY_BLOCK key_vec_packed = KEY_BLOCK_READ(key_input, key_offset + head_idx_index);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed) - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed)) * comp_scale;
#else
 KEY_BLOCK key_vals = key_vec_packed;
#endif
 uint query_offset = head_idx_index + sglid;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 QUERY_BLOCK query_vals_reg;
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 query_vals_reg[i] = query_local[query_offset + i * SUBGROUP_SIZE];
 }
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 acc[seq_idx] = mad(TO_SOFTMAX_ACCUMULATOR_TYPE(query_vals_reg[i]), TO_SOFTMAX_ACCUMULATOR_TYPE(key_vals[i]), acc[seq_idx]);
 }
 query_offset += HEAD_SIZE;
 }
 }
 #define KEY_BLOCK_SIZE 4
 for (; head_idx_index + (KEY_BLOCK_SIZE * SUBGROUP_SIZE) <= HEAD_SIZE; head_idx_index += SUBGROUP_SIZE * KEY_BLOCK_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, KEY_BLOCK_SIZE, ptr, offset);
 #define KEY_BLOCK MAKE_VECTOR_TYPE(INPUT1_TYPE, KEY_BLOCK_SIZE)
 #define KEY_BLOCK_UNCOMPRESSED MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, KEY_BLOCK_SIZE)
 #define TO_KEY_BLOCK_UNCOMPRESSED_TYPE(val) CAT(convert_, KEY_BLOCK_UNCOMPRESSED)(val)
 #define QUERY_BLOCK MAKE_VECTOR_TYPE(INPUT0_TYPE, KEY_BLOCK_SIZE)
 KEY_BLOCK key_vec_packed = KEY_BLOCK_READ(key_input, key_offset + head_idx_index);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed) - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed)) * comp_scale;
#else
 KEY_BLOCK key_vals = key_vec_packed;
#endif
 uint query_offset = head_idx_index + sglid;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 QUERY_BLOCK query_vals_reg;
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 query_vals_reg[i] = query_local[query_offset + i * SUBGROUP_SIZE];
 }
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 acc[seq_idx] = mad(TO_SOFTMAX_ACCUMULATOR_TYPE(query_vals_reg[i]), TO_SOFTMAX_ACCUMULATOR_TYPE(key_vals[i]), acc[seq_idx]);
 }
 query_offset += HEAD_SIZE;
 }
 }
 #define KEY_BLOCK_SIZE 2
 for (; head_idx_index + (KEY_BLOCK_SIZE * SUBGROUP_SIZE) <= HEAD_SIZE; head_idx_index += SUBGROUP_SIZE * KEY_BLOCK_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, KEY_BLOCK_SIZE, ptr, offset);
 #define KEY_BLOCK MAKE_VECTOR_TYPE(INPUT1_TYPE, KEY_BLOCK_SIZE)
 #define KEY_BLOCK_UNCOMPRESSED MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, KEY_BLOCK_SIZE)
 #define TO_KEY_BLOCK_UNCOMPRESSED_TYPE(val) CAT(convert_, KEY_BLOCK_UNCOMPRESSED)(val)
 #define QUERY_BLOCK MAKE_VECTOR_TYPE(INPUT0_TYPE, KEY_BLOCK_SIZE)
 KEY_BLOCK key_vec_packed = KEY_BLOCK_READ(key_input, key_offset + head_idx_index);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed) - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed)) * comp_scale;
#else
 KEY_BLOCK key_vals = key_vec_packed;
#endif
 uint query_offset = head_idx_index + sglid;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 QUERY_BLOCK query_vals_reg;
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 query_vals_reg[i] = query_local[query_offset + i * SUBGROUP_SIZE];
 }
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 acc[seq_idx] = mad(TO_SOFTMAX_ACCUMULATOR_TYPE(query_vals_reg[i]), TO_SOFTMAX_ACCUMULATOR_TYPE(key_vals[i]), acc[seq_idx]);
 }
 query_offset += HEAD_SIZE;
 }
 }
 #define KEY_BLOCK_SIZE 1
 for (; head_idx_index + (KEY_BLOCK_SIZE * SUBGROUP_SIZE) <= HEAD_SIZE; head_idx_index += SUBGROUP_SIZE * KEY_BLOCK_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, KEY_BLOCK_SIZE, ptr, offset);
 #define KEY_BLOCK MAKE_VECTOR_TYPE(INPUT1_TYPE, KEY_BLOCK_SIZE)
 #define KEY_BLOCK_UNCOMPRESSED MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, KEY_BLOCK_SIZE)
 #define TO_KEY_BLOCK_UNCOMPRESSED_TYPE(val) CAT(convert_, KEY_BLOCK_UNCOMPRESSED)(val)
 #define QUERY_BLOCK MAKE_VECTOR_TYPE(INPUT0_TYPE, KEY_BLOCK_SIZE)
 KEY_BLOCK key_vec_packed = KEY_BLOCK_READ(key_input, key_offset + head_idx_index);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed) - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed)) * comp_scale;
#else
 KEY_BLOCK key_vals = key_vec_packed;
#endif
 uint query_offset = head_idx_index + sglid;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 QUERY_BLOCK query_vals_reg;
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 query_vals_reg = query_local[query_offset + i * SUBGROUP_SIZE];
 }
 acc[seq_idx] = mad(TO_SOFTMAX_ACCUMULATOR_TYPE(query_vals_reg), TO_SOFTMAX_ACCUMULATOR_TYPE(key_vals), acc[seq_idx]);
 query_offset += HEAD_SIZE;
 }
 }
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc[seq_idx] = sub_group_reduce_add(acc[seq_idx]);
 qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len] = acc[seq_idx];
 }
 }
 {
 barrier(CLK_LOCAL_MEM_FENCE);
 SOFTMAX_ACCUMULATOR_TYPE qk_val[TARGET_SEQ_LEN_BLOCK_SIZE];
 const uint seq_idx_end = 1;
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 for (uint seq_len = sgid * SUBGROUP_SIZE + sglid; seq_len < partition_seq_len; seq_len += (HEAD_SIZE * SG_SCALE_FACTOR)) {
 qk_val[seq_idx] = qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len];
#if IS_CAUSAL
 if (start_partition_idx + seq_len > target_seq_idx + seq_idx)
 qk_val[seq_idx] += INPUT0_VAL_MIN;
#elif !IS_CAUSAL && HAS_ATTN_MASK_INPUT
 const uint attn_mask_offset = INPUT3_GET_INDEX_SAFE(b0_idx, b1_idx, target_seq_idx + seq_idx, start_partition_idx + seq_len);
 qk_val[seq_idx] += attn_mask[attn_mask_offset];
#endif
 qk_max[seq_idx] = SOFTMAX_ACCUMULATOR_MAX_FUNC(qk_max[seq_idx], TO_SOFTMAX_ACCUMULATOR_TYPE(qk_val[seq_idx]));
 qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len] = qk_val[seq_idx];
 }
 }
 }
 }
 {
 const uint seq_idx_end = 1;
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 qk_max[seq_idx] = sub_group_reduce_max(qk_max[seq_idx]);
 }
 if (sglid == 0) {
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 qk_max_vals[seq_idx * SUBGROUPS_PER_WG + sgid] = qk_max[seq_idx];
 }
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_max[seq_idx] = SOFTMAX_ACCUMULATOR_VAL_MIN;
 if (sglid < SUBGROUPS_PER_WG)
 qk_max[seq_idx] = qk_max_vals[seq_idx * SUBGROUPS_PER_WG + sglid];
 qk_max[seq_idx] = sub_group_reduce_max(qk_max[seq_idx]);
 }
 SOFTMAX_ACCUMULATOR_TYPE exp_sum[TARGET_SEQ_LEN_BLOCK_SIZE] = {SOFTMAX_ACCUMULATOR_VAL_ZERO};
 const uint qk_num_per_wi = CEIL_DIV(partition_seq_len, SUBGROUPS_PER_WG * SUBGROUP_SIZE);
 for (uint qk_idx = 0; qk_idx < qk_num_per_wi; qk_idx++) {
 const uint local_data_idx = qk_idx * (SUBGROUPS_PER_WG * SUBGROUP_SIZE) + lid;
 if (local_data_idx < partition_seq_len) {
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE qk_new = native_exp(TO_SOFTMAX_ACCUMULATOR_TYPE(qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + local_data_idx]) - qk_max[seq_idx]);
 qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + local_data_idx] = TO_OUTPUT_TYPE(qk_new);
 exp_sum[seq_idx] += qk_new;
 }
 }
 }
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 exp_sum[seq_idx] = sub_group_reduce_add(exp_sum[seq_idx]);
 if (sglid == 0)
 qk_sum_vals[seq_idx * SUBGROUPS_PER_WG + sgid] = exp_sum[seq_idx];
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 exp_sum[seq_idx] = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 if (sglid < SUBGROUPS_PER_WG)
 exp_sum[seq_idx] = qk_sum_vals[seq_idx * SUBGROUPS_PER_WG + sglid];
 exp_sum[seq_idx] = sub_group_reduce_add(exp_sum[seq_idx]);
 }
 for (uint qk_idx = 0; qk_idx < qk_num_per_wi; qk_idx++) {
 const uint local_data_idx = qk_idx * (SUBGROUPS_PER_WG * SUBGROUP_SIZE) + lid;
 if (local_data_idx < partition_seq_len) {
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE qk_new = TO_SOFTMAX_ACCUMULATOR_TYPE(qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + local_data_idx]) / exp_sum[seq_idx];
 qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + local_data_idx] = TO_OUTPUT_TYPE(qk_new);
 }
 }
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 {
 if (num_of_partitions > 1 && lid == 0) {
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 (seq_idx + target_seq_idx) * (num_of_partitions) +
 partition_idx;
 exp_sums[exp_sums_offset] = exp_sum[seq_idx];
 const uint max_logits_offset = exp_sums_offset;
 max_logits[max_logits_offset] = qk_max[seq_idx];
 }
 }
 }
 }
 }
 {
 OUTPUT_TYPE acc[TARGET_SEQ_LEN_BLOCK_SIZE] = {OUTPUT_VAL_ZERO};
#ifndef BEAM_TABLE_TYPE
#ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, 0, 0);
 uint value_offset_next_seq = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, 1, 0);
 const uint value_pitch = value_offset_next_seq - value_offset;
#else
 const uint value_pitch = HEAD_SIZE;
#endif
#endif
#if SG_SCALE_FACTOR > 1
 const uint seq_len_start = (sgid / (HEAD_SIZE / SUBGROUP_SIZE)) * (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR / SUBGROUP_SIZE);
 const uint seq_len_end = min(seq_len_start + (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR / SUBGROUP_SIZE), partition_seq_len / SUBGROUP_SIZE);
#else
 const uint seq_len_start = 0;
 const uint seq_len_end = partition_seq_len / SUBGROUP_SIZE;
#endif
 for (uint seq_len = seq_len_start; seq_len < seq_len_end; seq_len++) {
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, sgid * SUBGROUP_SIZE)];
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, sgid * SUBGROUP_SIZE);
#else
 const uint b_idx = b0_idx;
#ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE), head_size_idx);
#else
 uint value_offset = INPUT2_GET_INDEX(b_idx, b1_idx, start_partition_idx + (seq_len * SUBGROUP_SIZE), head_size_idx);
#endif
#endif
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 OUTPUT_TYPE qk_val[TARGET_SEQ_LEN_BLOCK_SIZE];
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len * SUBGROUP_SIZE + sglid];
 }
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, sub_group_broadcast(value_offset, i));
#else
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - sub_group_broadcast(comp_zp, i)) * sub_group_broadcast(comp_scale, i);
#elif IS_KV_COMPRESSED
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * sub_group_broadcast(comp_scale, i));
#else
 INPUT2_TYPE value_val = value_packed;
#endif
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc[seq_idx] = mad(sub_group_broadcast(qk_val[seq_idx], i), value_val, acc[seq_idx]);
 }
#ifndef BEAM_TABLE_TYPE
 value_offset += value_pitch;
#endif
 }
 }
#if SG_SCALE_FACTOR > 1
 if (sgid >= HEAD_SIZE / SUBGROUP_SIZE) {
#endif
 for (uint seq_len = (partition_seq_len / SUBGROUP_SIZE) * SUBGROUP_SIZE; seq_len < partition_seq_len; seq_len++) {
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + seq_len, head_size_idx)];
#else
 const uint b_idx = b0_idx;
#endif
#ifdef INPUT2_DIMS_ORDER
 const uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + seq_len, head_size_idx);
#else
 const uint value_offset = INPUT2_GET_INDEX(b_idx, b1_idx, start_partition_idx + seq_len, head_size_idx);
#endif
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + seq_len, 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 OUTPUT_TYPE qk_val[TARGET_SEQ_LEN_BLOCK_SIZE];
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len];
 }
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 const VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 const VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * comp_scale);
#else
 const INPUT2_TYPE value_val = value_packed;
#endif
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc[seq_idx] = mad(qk_val[seq_idx], value_val, acc[seq_idx]);
 }
 }
#if SG_SCALE_FACTOR > 1
 }
#endif
#if SG_SCALE_FACTOR > 1
 if ((partition_seq_len > (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR)) || (partition_seq_len % SUBGROUP_SIZE != 0)) {
 if (sgid >= HEAD_SIZE / SUBGROUP_SIZE) {
 query_local[head_size_idx] = acc[0];
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 if (sgid < HEAD_SIZE / SUBGROUP_SIZE) {
 acc[0] += query_local[head_size_idx];
 }
 }
#endif
#if SG_SCALE_FACTOR > 1
 if (sgid < HEAD_SIZE / SUBGROUP_SIZE) {
#endif
 if (num_of_partitions > 1) {
 const uint seq_idx_end = 1;
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 const uint tmp_out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 (target_seq_idx + seq_idx) * (num_of_partitions * HEAD_SIZE) +
 partition_idx * (HEAD_SIZE) +
 head_size_idx;
 tmp_out[tmp_out_offset] = acc[seq_idx];
 }
 } else {
 const uint seq_idx_end = 1;
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 const uint output_offset = OUTPUT_GET_INDEX(b0_idx, b1_idx, target_seq_idx + seq_idx, head_size_idx);
 output[output_offset] = acc[seq_idx];
 }
 }
#if SG_SCALE_FACTOR > 1
 }
#endif
 }
}
#else
#if IS_PAGED_ATTENTION
 #define SOURCE_SEQ_LEN (subsequence_begins[gws_seq_indexes_correspondence[((uint)get_global_id(1))] + 1] - subsequence_begins[gws_seq_indexes_correspondence[((uint)get_global_id(1))]])
 #define TARGET_SEQ_LEN (subsequence_begins[gws_seq_indexes_correspondence[((uint)get_global_id(1))] + 1] - subsequence_begins[gws_seq_indexes_correspondence[((uint)get_global_id(1))]])
 #define PA_BUFFERS , subsequence_begins, blocked_indexes_start, blocked_indexes_end, gws_seq_indexes_correspondence
 #define PA_BUFFERS_ARGS , const __global INPUT3_TYPE* subsequence_begins, const __global int* blocked_indexes_start, const __global int* blocked_indexes_end, const __global int* gws_seq_indexes_correspondence
#else
 #define PA_BUFFERS
 #define PA_BUFFERS_ARGS
#endif
#if HAS_ATTN_MASK_INPUT
 #define ATTN_MASK_BUFFER , attn_mask
 #define ATTN_MASK_BUFFER_ARG , const __global INPUT3_TYPE* attn_mask
#else
 #define ATTN_MASK_BUFFER
 #define ATTN_MASK_BUFFER_ARG
#endif
#if HAS_SCALE_INPUT
 #define ATTN_SCALE_BUFFER , scale
 #define ATTN_SCALE_BUFFER_ARG , const __global INPUT4_TYPE* scale
#else
 #define ATTN_SCALE_BUFFER
 #define ATTN_SCALE_BUFFER_ARG
#endif
#if IS_KV_COMPRESSED
#define APPLY_SCALES_TO_QUERY 1
#endif
#define MASK_VECTOR_TYPE MAKE_VECTOR_TYPE(INPUT0_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
inline MASK_VECTOR_TYPE FUNC(load_attn_mask)(OPTIONAL_SHAPE_INFO_ARG
 uint b0_idx,
 uint b1_idx,
 uint target_seq_idx,
 uint source_seq_idx
 ATTN_MASK_BUFFER_ARG
 ATTN_SCALE_BUFFER_ARG
 PA_BUFFERS_ARGS
 ) {
 MASK_VECTOR_TYPE mask_vec = INPUT0_VAL_ZERO;
#if !IS_CAUSAL && HAS_ATTN_MASK_INPUT
 const uint attn_mask_offset = INPUT3_GET_INDEX_SAFE(b0_idx, b1_idx, target_seq_idx, source_seq_idx);
 if (target_seq_idx >= (uint)TARGET_SEQ_LEN) {
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 mask_vec[i] = NAN;
 }
 } else {
 if (source_seq_idx + SUBGROUP_SIZE <= (uint)SOURCE_SEQ_LEN) {
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 const INPUT3_TYPE mask_val = attn_mask[attn_mask_offset + i];
 mask_vec[i] = mask_val;
 }
 } else {
 const uint max_mask_offset = min(source_seq_idx + SUBGROUP_SIZE, (uint)SOURCE_SEQ_LEN);
 for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 const INPUT3_TYPE mask_val = source_seq_idx + i < max_mask_offset ? attn_mask[attn_mask_offset + i] : NAN;
 mask_vec[i] = mask_val;
 }
 }
 }
#endif
#if !IS_CAUSAL && !HAS_ATTN_MASK_INPUT
 if (target_seq_idx >= (uint)TARGET_SEQ_LEN) {
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 mask_vec[i] = NAN;
 }
 } else {
 const uint max_mask_offset = min(source_seq_idx + SUBGROUP_SIZE, (uint)SOURCE_SEQ_LEN);
 for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 mask_vec[i] = source_seq_idx + i < max_mask_offset ? 0 : NAN;
 }
 }
#endif
#if IS_CAUSAL
 if (target_seq_idx >= (uint)TARGET_SEQ_LEN) {
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 mask_vec[i] = NAN;
 }
 } else {
 for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 if (source_seq_idx + i > target_seq_idx)
 mask_vec[i] = NAN;
 }
 }
#endif
#if HAS_SCALE_INPUT
 const OUTPUT_TYPE scale_val = OUTPUT_VAL_ONE / *scale;
#else
 const INPUT0_TYPE scale_val = TO_INPUT0_TYPE(STATIC_SCALE_VALUE_INV);
#endif
#if IS_CAUSAL || HAS_ATTN_MASK_INPUT
 mask_vec *= scale_val;
#endif
 return mask_vec;
}
#if IS_PAGED_ATTENTION && HAS_ALIBI
#if HAS_SCALE_INPUT
#define ALIBI_TYPE INPUT5_TYPE
#else
#define ALIBI_TYPE INPUT4_TYPE
#endif
#endif
REQD_SUB_GROUP_SIZE(SUBGROUP_SIZE)
KERNEL(sdpa_opt)(
 OPTIONAL_SHAPE_INFO_ARG
 const __global INPUT0_TYPE* query_input,
 const __global INPUT1_TYPE* key_input,
 const __global INPUT2_TYPE* value_input,
#if IS_PAGED_ATTENTION
 const __global INPUT3_TYPE* subsequence_begins,
#endif
#if HAS_ATTN_MASK_INPUT
 const __global INPUT3_TYPE* attn_mask,
#endif
#if HAS_SCALE_INPUT
 const __global INPUT4_TYPE* scale,
#endif
#if IS_PAGED_ATTENTION && HAS_ALIBI
 const __global ALIBI_TYPE* alibi_slopes,
#endif
 __global OUTPUT_TYPE* output,
#if IS_KV_COMPRESSED
 const __global KEY_COMPRESSION_SCALE_TYPE* key_scale,
 const __global VALUE_COMPRESSION_SCALE_TYPE* val_scale,
#endif
#ifdef BEAM_TABLE_TYPE
 const __global BEAM_TABLE_TYPE* beam_table,
#endif
#if IS_PAGED_ATTENTION
 const __global int* blocked_indexes_start,
 const __global int* blocked_indexes_end,
 const __global int* gws_seq_indexes_correspondence
#else
 __global SOFTMAX_ACCUMULATOR_TYPE* exp_sums,
 __global SOFTMAX_ACCUMULATOR_TYPE* max_logits,
 __global OUTPUT_TYPE* tmp_out
#endif
)
{
#if TARGET_SEQ_LEN_BLOCK_SIZE != 16
 #error sdpa_opt.cl: unsupported TARGET_SEQ_LEN_BLOCK_SIZE
#endif
 #define batch_idx ((uint)get_global_id(0))
 #define num_heads_dim ((uint)get_global_id(0))
 #define b0_idx (batch_idx / NUM_HEADS)
 #define b1_idx (batch_idx % NUM_HEADS)
 #define target_seq_dim ((uint)get_global_id(1))
 #define target_seq_idx ((uint)get_global_id(1) * TARGET_SEQ_LEN_BLOCK_SIZE)
 #define head_size_idx ((uint)get_local_id(2) % HEAD_SIZE)
 #define sglid (uint)get_sub_group_local_id()
 #define sgid (uint)get_sub_group_id()
 __local INPUT0_TYPE slm_query[HEAD_SIZE * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local OUTPUT_TYPE slm_qk_vals[SEQ_LEN_PARTITION_SIZE * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_qk_max_vals[SUBGROUPS_PER_WG * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_exp_sum_vals[SUBGROUPS_PER_WG * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_exp_sum_cur[TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_max_val_cur[TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_exp_sum_prev[TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_max_val_prev[TARGET_SEQ_LEN_BLOCK_SIZE];
 {
#if IS_PAGED_ATTENTION
 const uint block_start_pos = blocked_indexes_start[target_seq_dim];
 const uint block_end_pos = blocked_indexes_end[target_seq_dim];
 uint query_offset = INPUT0_OFFSET +
 block_start_pos * (HEAD_SIZE * NUM_HEADS + INPUT0_PAD_BEFORE_FEATURE_NUM + INPUT0_PAD_AFTER_FEATURE_NUM) +
 num_heads_dim * HEAD_SIZE + head_size_idx;
 const uint query_pitch = (HEAD_SIZE * NUM_HEADS + INPUT0_PAD_BEFORE_FEATURE_NUM + INPUT0_PAD_AFTER_FEATURE_NUM);
 const uint cur_target_seq_len_size = block_end_pos - block_start_pos;
#else
#ifdef INPUT0_DIMS_ORDER
 uint query_offset = FUNC_CALL(get_input0_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, target_seq_idx, (head_size_idx));
 uint query_offset_next_seq = FUNC_CALL(get_input0_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, target_seq_idx + 1, (head_size_idx));
 const uint query_pitch = query_offset_next_seq - query_offset;
#else
 uint query_offset = INPUT0_GET_INDEX(b0_idx, b1_idx, target_seq_idx, (head_size_idx));
 const uint query_pitch = HEAD_SIZE;
#endif
 const uint cur_target_seq_len_size = min(TARGET_SEQ_LEN - target_seq_idx, (uint)TARGET_SEQ_LEN_BLOCK_SIZE);
#endif
 uint query_local_offset = head_size_idx * TARGET_SEQ_LEN_BLOCK_SIZE;
#if APPLY_SCALES_TO_QUERY
#if HAS_SCALE_INPUT
 const INPUT0_TYPE scale_val = *scale;
#else
 const INPUT0_TYPE scale_val = TO_INPUT0_TYPE(STATIC_SCALE_VALUE);
#endif
#else
 const INPUT0_TYPE scale_val = INPUT0_VAL_ONE;
#endif
 if (cur_target_seq_len_size != TARGET_SEQ_LEN_BLOCK_SIZE) {
 if (sgid * SUBGROUP_SIZE < HEAD_SIZE) {
 for (uint seq_idx = 0; seq_idx < cur_target_seq_len_size; seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 }
 } else {
 #if SG_SCALE_FACTOR == 2
 if ((sgid < (SUBGROUPS_PER_WG / SG_SCALE_FACTOR))) {
 unroll_for (uint seq_idx = 0; seq_idx < (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR); seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 } else {
 query_local_offset += (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR);
 query_offset += query_pitch * (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR);
 unroll_for (uint seq_idx = 0; seq_idx < (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR); seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 }
 #elif SG_SCALE_FACTOR == 4
 query_local_offset += (sgid / (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) * (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR);
 query_offset += query_pitch * (sgid / (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) * (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR);
 unroll_for (uint seq_idx = 0; seq_idx < (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR); seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 #else
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 #endif
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 }
 {
 #if TARGET_SEQ_LEN_BLOCK_SIZE <= SUBGROUP_SIZE
 if (sgid == 0 && sglid < TARGET_SEQ_LEN_BLOCK_SIZE) {
 slm_max_val_prev[sglid] = SOFTMAX_ACCUMULATOR_VAL_MIN;
 slm_exp_sum_prev[sglid] = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 }
 #else
 #error sdpa_opt.cl: unsupported TARGET_SEQ_LEN_BLOCK_SIZE
 #endif
 }
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) output_acc = OUTPUT_VAL_ZERO;
 __attribute__((opencl_unroll_hint(1)))
 for (uint start_partition_idx = 0; start_partition_idx < SOURCE_SEQ_LEN; start_partition_idx += SEQ_LEN_PARTITION_SIZE) {
 SOFTMAX_ACCUMULATOR_TYPE qk_max = SOFTMAX_ACCUMULATOR_VAL_MIN;
 const uint seq_len = start_partition_idx + sgid * SUBGROUP_SIZE;
 const uint partition_seq_len = min((uint)SOURCE_SEQ_LEN - start_partition_idx, (uint)SEQ_LEN_PARTITION_SIZE);
#if IS_PAGED_ATTENTION
#ifdef BROADCAST_GROUP_SIZE
 const uint heads_dim = num_heads_dim / BROADCAST_GROUP_SIZE;
#else
 const uint heads_dim = num_heads_dim;
#endif
 #define KEY_SEQ_OFFSET subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]]
 const uint key_pitch = (HEAD_SIZE * NUM_KV_HEADS + INPUT1_PAD_BEFORE_FEATURE_NUM + INPUT1_PAD_AFTER_FEATURE_NUM);
 uint key_offset = INPUT1_OFFSET +
 KEY_SEQ_OFFSET * key_pitch +
 heads_dim * HEAD_SIZE +
 seq_len * key_pitch;
#else
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_key)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, seq_len + sglid, 0)];
 const uint key_offset = FUNC_CALL(get_input1_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, seq_len + sglid, 0);
#else
 const uint b_idx = b0_idx;
 #ifdef INPUT1_DIMS_ORDER
 uint key_offset = FUNC_CALL(get_input1_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, seq_len, 0);
 uint key_offset_next_seq = FUNC_CALL(get_input1_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, seq_len + 1, 0);
 const uint key_pitch = key_offset_next_seq - key_offset;
 #else
 uint key_offset = INPUT1_GET_INDEX(b0_idx, b1_idx, seq_len, 0);
 const uint key_pitch = HEAD_SIZE;
 #endif
#endif
#endif
 int seq_len_calc_size = min((int)(SOURCE_SEQ_LEN) - (int)seq_len, (int)SUBGROUP_SIZE);
 MAKE_VECTOR_TYPE(INPUT0_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) qk_acc;
 qk_acc = FUNC_CALL(load_attn_mask)(OPTIONAL_SHAPE_INFO_TENSOR
 b0_idx,
 b1_idx,
#if IS_PAGED_ATTENTION
 blocked_indexes_start[target_seq_dim] - subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]] + sglid,
#else
 target_seq_idx + sglid,
#endif
 seq_len
 ATTN_MASK_BUFFER
 ATTN_SCALE_BUFFER
 PA_BUFFERS);
 if (seq_len_calc_size >= SUBGROUP_SIZE) {
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(KEY_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, seq_len + sglid, 0);
 KEY_COMPRESSION_SCALE_TYPE comp_scale = key_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 KEY_COMPRESSION_SCALE_TYPE comp_zp = key_scale[comp_offset + 1];
#endif
#endif
 __attribute__((opencl_unroll_hint(1)))
 for (uint head_idx_index = 0; head_idx_index < HEAD_SIZE; head_idx_index += SUBGROUP_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, 1, ptr, offset);
 #define QUERY_VEC MAKE_VECTOR_TYPE(INPUT0_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
 QUERY_VEC queries_vec;
 uint query_local_offset = (head_idx_index * TARGET_SEQ_LEN_BLOCK_SIZE) + sglid;
 unroll_for (uint q_row_idx = 0; q_row_idx < TARGET_SEQ_LEN_BLOCK_SIZE; q_row_idx++) {
 queries_vec[q_row_idx] = slm_query[query_local_offset];
 query_local_offset += TARGET_SEQ_LEN_BLOCK_SIZE;
 }
 unroll_for (uint key_row_idx = 0; key_row_idx < TARGET_SEQ_LEN_BLOCK_SIZE; key_row_idx++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT1_TYPE key_packed = KEY_BLOCK_READ(key_input, sub_group_broadcast(key_offset, key_row_idx) + head_idx_index);
#else
 const INPUT1_TYPE key_packed = KEY_BLOCK_READ(key_input, key_offset + key_row_idx * key_pitch + head_idx_index);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_COMPRESSION_SCALE_TYPE key_vals = (TO_KEY_COMPRESSION_SCALE_TYPE(key_packed) - sub_group_broadcast(comp_zp, key_row_idx)) * sub_group_broadcast(comp_scale, key_row_idx);
#elif IS_KV_COMPRESSED
 KEY_COMPRESSION_SCALE_TYPE key_vals = (TO_KEY_COMPRESSION_SCALE_TYPE(key_packed) * sub_group_broadcast(comp_scale, key_row_idx));
#else
 INPUT1_TYPE key_vals = key_packed;
#endif
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 qk_acc[key_row_idx] = mad(sub_group_broadcast(key_vals, i), queries_vec[i], qk_acc[key_row_idx]);
 }
 }
 }
 } else if (seq_len_calc_size > 0) {
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(KEY_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, seq_len + min(sglid, (uint)seq_len_calc_size - 1), 0);
 KEY_COMPRESSION_SCALE_TYPE comp_scale = key_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 KEY_COMPRESSION_SCALE_TYPE comp_zp = key_scale[comp_offset + 1];
#endif
#endif
 __attribute__((opencl_unroll_hint(1)))
 for (uint head_idx_index = 0; head_idx_index < HEAD_SIZE; head_idx_index += SUBGROUP_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, 1, ptr, offset)
 #define QUERY_VEC_TYPE MAKE_VECTOR_TYPE(INPUT0_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
#if IS_KV_COMPRESSED
 #define KEY_UNPACKED_TYPE KEY_COMPRESSION_SCALE_TYPE
 #define KEY_UNPACKED_VEC_TYPE MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
 #define TO_KEY_UNPACKED_TYPE(val) TO_KEY_COMPRESSION_SCALE_TYPE(val)
#else
 #define KEY_UNPACKED_TYPE INPUT1_TYPE
 #define KEY_UNPACKED_VEC_TYPE MAKE_VECTOR_TYPE(INPUT1_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
 #define TO_KEY_UNPACKED_TYPE(val) TO_INPUT1_TYPE(val)
#endif
 QUERY_VEC_TYPE queries_vec;
 uint query_local_offset = (head_idx_index * TARGET_SEQ_LEN_BLOCK_SIZE) + sglid;
 unroll_for (uint q_row_idx = 0; q_row_idx < TARGET_SEQ_LEN_BLOCK_SIZE; q_row_idx++) {
 queries_vec[q_row_idx] = slm_query[query_local_offset];
 query_local_offset += TARGET_SEQ_LEN_BLOCK_SIZE;
 }
#ifndef LOAD_KEY_LEFTOVERS_IN_CALC_LOOP
 KEY_UNPACKED_VEC_TYPE key_vec = 0;
 unroll_for (uint key_row_idx = 0; key_row_idx < seq_len_calc_size; key_row_idx++) {
#ifdef BEAM_TABLE_TYPE
 key_vec[key_row_idx] = TO_KEY_UNPACKED_TYPE(KEY_BLOCK_READ(key_input, sub_group_broadcast(key_offset, key_row_idx) + head_idx_index));
#else
 key_vec[key_row_idx] = TO_KEY_UNPACKED_TYPE(KEY_BLOCK_READ(key_input, key_offset + key_row_idx * key_pitch + head_idx_index));
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 key_vec[key_row_idx] = (key_vec[key_row_idx] - sub_group_broadcast(comp_zp, key_row_idx)) * sub_group_broadcast(comp_scale, key_row_idx);
#elif IS_KV_COMPRESSED
 key_vec[key_row_idx] *= sub_group_broadcast(comp_scale, key_row_idx);
#endif
 }
#endif
 unroll_for (uint key_row_idx = 0; key_row_idx < TARGET_SEQ_LEN_BLOCK_SIZE; key_row_idx++) {
#ifdef LOAD_KEY_LEFTOVERS_IN_CALC_LOOP
 KEY_UNPACKED_TYPE key_vals = 0;
 if (key_row_idx < seq_len_calc_size) {
#ifdef BEAM_TABLE_TYPE
 key_vals = TO_KEY_UNPACKED_TYPE(KEY_BLOCK_READ(key_input, sub_group_broadcast(key_offset, key_row_idx) + head_idx_index));
#else
 key_vals = TO_KEY_UNPACKED_TYPE(KEY_BLOCK_READ(key_input, key_offset + key_row_idx * key_pitch + head_idx_index));
#endif
 }
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 key_vals = (key_vals - sub_group_broadcast(comp_zp, key_row_idx)) * sub_group_broadcast(comp_scale, key_row_idx);
#elif IS_KV_COMPRESSED
 key_vals *= sub_group_broadcast(comp_scale, key_row_idx);
#endif
#else
 #define key_vals key_vec[key_row_idx]
#endif
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 qk_acc[key_row_idx] = mad(sub_group_broadcast(key_vals, i), queries_vec[i], qk_acc[key_row_idx]);
 }
 }
 }
 }
 {
 unroll_for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
#if !APPLY_SCALES_TO_QUERY
#if HAS_SCALE_INPUT
 const OUTPUT_TYPE scale_val = *scale;
#else
 const OUTPUT_TYPE scale_val = TO_OUTPUT_TYPE(STATIC_SCALE_VALUE);
#endif
 qk_acc[i] *= scale_val;
#endif
#ifdef HAS_ALIBI
 const int alibi_val = (1 - SOURCE_SEQ_LEN) + seq_len + i;
 qk_acc[i] += alibi_slopes[num_heads_dim] * alibi_val;
#endif
 qk_acc[i] = INPUT0_MIN_FUNC(INPUT0_MAX_FUNC(qk_acc[i], INPUT0_VAL_MIN), INPUT0_VAL_MAX);
 qk_max = SOFTMAX_ACCUMULATOR_MAX_FUNC(qk_max, TO_SOFTMAX_ACCUMULATOR_TYPE(qk_acc[i]));
 }
 }
 {
 slm_qk_max_vals[sgid * SUBGROUP_SIZE + sglid] = qk_max;
 qk_max = SOFTMAX_ACCUMULATOR_VAL_MIN;
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 {
 SOFTMAX_ACCUMULATOR_TYPE qk_max_new = SOFTMAX_ACCUMULATOR_VAL_MIN;
 for (uint i = 0; i < SUBGROUPS_PER_WG; i++) {
 SOFTMAX_ACCUMULATOR_TYPE qk_max_val = slm_qk_max_vals[i * SUBGROUP_SIZE + sglid];
 qk_max_new = SOFTMAX_ACCUMULATOR_MAX_FUNC(qk_max_new, qk_max_val);
 }
 if (sgid == 0) {
 slm_max_val_cur[sglid] = qk_max_new;
 }
 SOFTMAX_ACCUMULATOR_TYPE exp_sum_new = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
 qk_acc[i] = native_exp(TO_SOFTMAX_ACCUMULATOR_TYPE(qk_acc[i]) - qk_max_new);
 exp_sum_new += qk_acc[i];
 }
 {
 slm_exp_sum_vals[sgid * SUBGROUP_SIZE + sglid] = exp_sum_new;
 }
 exp_sum_new = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 barrier(CLK_LOCAL_MEM_FENCE);
 for (uint i = 0; i < SUBGROUPS_PER_WG; i++) {
 SOFTMAX_ACCUMULATOR_TYPE exp_sum = slm_exp_sum_vals[i * SUBGROUP_SIZE + sglid];
 exp_sum_new += exp_sum;
 }
 for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
 qk_acc[i] = qk_acc[i] / exp_sum_new;
 }
 if (sgid == 0) {
 slm_exp_sum_cur[sglid] = exp_sum_new;
 }
 for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
 slm_qk_vals[sglid * SEQ_LEN_PARTITION_SIZE + sgid * TARGET_SEQ_LEN_BLOCK_SIZE + i] = qk_acc[i];
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 }
 {
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) acc_output_res = OUTPUT_VAL_ZERO;
#if IS_PAGED_ATTENTION
 const uint value_pitch = (HEAD_SIZE * NUM_KV_HEADS + INPUT2_PAD_BEFORE_FEATURE_NUM + INPUT2_PAD_AFTER_FEATURE_NUM);
#else
#ifdef INPUT2_DIMS_ORDER
 uint value_offset_base = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, 0, 0);
 uint value_offset_next_seq = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, 1, 0);
 const uint value_pitch = value_offset_next_seq - value_offset_base;
#else
 const uint value_pitch = HEAD_SIZE;
#endif
#endif
 if (partition_seq_len == SEQ_LEN_PARTITION_SIZE) {
 uint seq_len_start = (sgid / (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) * (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR);
 for (uint seq_len = seq_len_start; seq_len < seq_len_start + (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR); seq_len += SUBGROUP_SIZE) {
#if IS_PAGED_ATTENTION
#ifdef BROADCAST_GROUP_SIZE
 const uint heads_dim = num_heads_dim / BROADCAST_GROUP_SIZE;
#else
 const uint heads_dim = num_heads_dim;
#endif
 const uint value_seq_offset = subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]];
 uint value_offset = INPUT2_OFFSET +
 value_seq_offset * value_pitch +
 heads_dim * HEAD_SIZE +
 (start_partition_idx + (seq_len)) * value_pitch + head_size_idx;
#else
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len) + sglid, sgid * SUBGROUP_SIZE)];
 const uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + (seq_len) + sglid, sgid * SUBGROUP_SIZE);
#else
 const uint b_idx = b0_idx;
 #ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len), head_size_idx);
 #else
 uint value_offset = INPUT2_GET_INDEX(b0_idx, b1_idx, start_partition_idx + (seq_len), head_size_idx);
 #endif
#endif
#endif
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) qk_val;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = slm_qk_vals[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len + sglid];
 }
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + seq_len + sglid, 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, sub_group_broadcast(value_offset, i));
#else
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - sub_group_broadcast(comp_zp, i)) * sub_group_broadcast(comp_scale, i);
#elif IS_KV_COMPRESSED
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * sub_group_broadcast(comp_scale, i));
#else
 INPUT2_TYPE value_val = value_packed;
#endif
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc_output_res[seq_idx] = mad(sub_group_broadcast(qk_val[seq_idx], i), value_val, acc_output_res[seq_idx]);
 }
#ifndef BEAM_TABLE_TYPE
 value_offset += value_pitch;
#endif
 }
 }
 } else {
 const uint seq_len_start = (sgid / (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) * (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR);
 uint seq_len_end = 0;
 if (seq_len_start < partition_seq_len)
 seq_len_end = seq_len_start + min(partition_seq_len - seq_len_start, (uint)(SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR));;
 for (uint seq_len = seq_len_start / SUBGROUP_SIZE; seq_len < seq_len_end / SUBGROUP_SIZE; seq_len++) {
#if IS_PAGED_ATTENTION
#ifdef BROADCAST_GROUP_SIZE
 const uint heads_dim = num_heads_dim / BROADCAST_GROUP_SIZE;
#else
 const uint heads_dim = num_heads_dim;
#endif
 const uint value_seq_offset = subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]];
 uint value_offset = INPUT2_OFFSET +
 value_seq_offset * value_pitch +
 heads_dim * HEAD_SIZE +
 (start_partition_idx + (seq_len * SUBGROUP_SIZE)) * value_pitch + head_size_idx;
#else
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, sgid * SUBGROUP_SIZE)];
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, sgid * SUBGROUP_SIZE);
#else
 const uint b_idx = b0_idx;
 #ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE), head_size_idx);
 #else
 uint value_offset = INPUT2_GET_INDEX(b0_idx, b1_idx, start_partition_idx + (seq_len * SUBGROUP_SIZE), head_size_idx);
 #endif
#endif
#endif
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) qk_val;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = slm_qk_vals[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len * SUBGROUP_SIZE + sglid];
 }
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, sub_group_broadcast(value_offset, i));
#else
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - sub_group_broadcast(comp_zp, i)) * sub_group_broadcast(comp_scale, i);
#elif IS_KV_COMPRESSED
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * sub_group_broadcast(comp_scale, i));
#else
 INPUT2_TYPE value_val = value_packed;
#endif
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc_output_res[seq_idx] = mad(sub_group_broadcast(qk_val[seq_idx], i), value_val, acc_output_res[seq_idx]);
 }
#ifndef BEAM_TABLE_TYPE
 value_offset += value_pitch;
#endif
 }
 }
 const uint seq_len_leftovers_start = ((seq_len_end / SUBGROUP_SIZE) * SUBGROUP_SIZE);
 if (seq_len_leftovers_start != seq_len_end) {
 uint qk_offset = min(seq_len_leftovers_start + sglid, seq_len_end - 1);
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) qk_val;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = slm_qk_vals[qk_offset];
 qk_offset += SEQ_LEN_PARTITION_SIZE;
 }
#if IS_PAGED_ATTENTION
#ifdef BROADCAST_GROUP_SIZE
 const uint heads_dim = num_heads_dim / BROADCAST_GROUP_SIZE;
#else
 const uint heads_dim = num_heads_dim;
#endif
 const uint value_seq_offset = subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]];
 uint value_offset = INPUT2_OFFSET +
 value_seq_offset * value_pitch +
 heads_dim * HEAD_SIZE +
 (start_partition_idx + seq_len_leftovers_start) * value_pitch + head_size_idx;
#else
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + seq_len_leftovers_start + sglid, sgid * SUBGROUP_SIZE)];
 const uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + seq_len_leftovers_start + sglid, sgid * SUBGROUP_SIZE);
#else
 const uint b_idx = b0_idx;
 #ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + seq_len_leftovers_start, head_size_idx);
 #else
 uint value_offset = INPUT2_GET_INDEX(b0_idx, b1_idx, start_partition_idx + seq_len_leftovers_start, head_size_idx);
 #endif
#endif
#endif
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + min(seq_len_leftovers_start + sglid, seq_len_end - 1), 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 for (uint seq_len_idx = 0; seq_len_idx < partition_seq_len - seq_len_leftovers_start; seq_len_idx++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, sub_group_broadcast(value_offset, seq_len_idx));
#else
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - sub_group_broadcast(comp_zp, seq_len_idx)) * sub_group_broadcast(comp_scale, seq_len_idx);
#elif IS_KV_COMPRESSED
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * sub_group_broadcast(comp_scale, seq_len_idx));
#else
 INPUT2_TYPE value_val = value_packed;
#endif
 for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc_output_res[seq_idx] = mad(sub_group_broadcast(qk_val[seq_idx], seq_len_idx), value_val, acc_output_res[seq_idx]);
 }
#ifndef BEAM_TABLE_TYPE
 value_offset += value_pitch;
#endif
 }
 }
 }
 {
 SOFTMAX_ACCUMULATOR_TYPE exp_sum_prev = slm_exp_sum_prev[sglid];
 SOFTMAX_ACCUMULATOR_TYPE exp_sum_cur = slm_exp_sum_cur[sglid];
 SOFTMAX_ACCUMULATOR_TYPE max_val_prev = slm_max_val_prev[sglid];
 SOFTMAX_ACCUMULATOR_TYPE max_val_cur = slm_max_val_cur[sglid];
 barrier(CLK_LOCAL_MEM_FENCE);
#if IS_PAGED_ATTENTION
 const uint block_start_pos_new = blocked_indexes_start[target_seq_dim];
 const uint block_end_pos_new = blocked_indexes_end[target_seq_dim];
 const uint seq_idx_end = block_end_pos_new - block_start_pos_new;
#else
 const uint seq_idx_end = min(TARGET_SEQ_LEN - target_seq_idx, (uint)TARGET_SEQ_LEN_BLOCK_SIZE);
#endif
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE total_max = SOFTMAX_ACCUMULATOR_MAX_FUNC(sub_group_broadcast(max_val_prev, seq_idx), sub_group_broadcast(max_val_cur, seq_idx));
 SOFTMAX_ACCUMULATOR_TYPE updated_exp_sum_prev = sub_group_broadcast(exp_sum_prev, seq_idx) * native_exp(sub_group_broadcast(max_val_prev, seq_idx) - total_max);
 SOFTMAX_ACCUMULATOR_TYPE updated_exp_sum_cur = sub_group_broadcast(exp_sum_cur, seq_idx) * native_exp(sub_group_broadcast(max_val_cur, seq_idx) - total_max);
 SOFTMAX_ACCUMULATOR_TYPE updated_total_exp_sum = updated_exp_sum_prev + updated_exp_sum_cur;
 if (start_partition_idx > 0) {
 OUTPUT_TYPE updated_prev_res = TO_SOFTMAX_ACCUMULATOR_TYPE(output_acc[seq_idx]) * updated_exp_sum_prev / updated_total_exp_sum;;
 acc_output_res[seq_idx] *= updated_exp_sum_cur / updated_total_exp_sum;
 acc_output_res[seq_idx] += updated_prev_res;
 }
 output_acc[seq_idx] = acc_output_res[seq_idx];
 if (sgid == 0 && sglid == 0) {
 slm_exp_sum_prev[seq_idx] = updated_total_exp_sum;
 slm_max_val_prev[seq_idx] = total_max;
 }
 }
 }
 }
 }
 if (sgid >= (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) {
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 slm_qk_vals[seq_idx * SEQ_LEN_PARTITION_SIZE + (uint)get_local_id(2)] = output_acc[seq_idx];
 }
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 if (sgid < (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) {
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 unroll_for (uint i = 1; i < SG_SCALE_FACTOR; i++) {
 output_acc[seq_idx] += slm_qk_vals[seq_idx * SEQ_LEN_PARTITION_SIZE + (i * HEAD_SIZE) + head_size_idx];
 }
 }
#if IS_PAGED_ATTENTION
 const uint block_start_pos_new = blocked_indexes_start[target_seq_dim];
 const uint block_end_pos_new = blocked_indexes_end[target_seq_dim];
 uint output_offset = block_start_pos_new * HEAD_SIZE * NUM_HEADS + num_heads_dim * HEAD_SIZE + sgid * SUBGROUP_SIZE;
 const uint output_pitch = HEAD_SIZE * NUM_HEADS;
#else
 uint output_offset = OUTPUT_GET_INDEX(b0_idx, b1_idx, target_seq_idx, sgid * SUBGROUP_SIZE);
 const uint output_pitch = HEAD_SIZE;
#endif
#if IS_PAGED_ATTENTION
 if (block_start_pos_new + TARGET_SEQ_LEN_BLOCK_SIZE != block_end_pos_new) {
 const uint seq_idx_end = block_end_pos_new - block_start_pos_new;
#else
 if (get_global_id(1) == get_global_size(1) - 1) {
 const uint seq_idx_end = min((uint)TARGET_SEQ_LEN - target_seq_idx, (uint)TARGET_SEQ_LEN_BLOCK_SIZE);
#endif
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 OUTPUT_BLOCK_WRITE(output, output_offset, output_acc[seq_idx]);
 output_offset += output_pitch;
 }
 } else {
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 OUTPUT_BLOCK_WRITE(output, output_offset, output_acc[seq_idx]);
 output_offset += output_pitch;
 }
 }
 }
}
#endif
#endif
#ifdef SDPA_STAGE_1
#if SOFTMAX_ACCUMULATOR_TYPE_SIZE == 4
#define REG_VERSION_MAX_VALUES_PER_WI 24
#define REG_VERSION_MAX_VALUES_PER_WI_LOWER 8
#elif SOFTMAX_ACCUMULATOR_TYPE_SIZE == 2
#define REG_VERSION_MAX_VALUES_PER_WI 48
#define REG_VERSION_MAX_VALUES_PER_WI_LOWER 16
#else
#error Unexpected SOFTMAX_ACCUMULATOR data type size
#endif
REQD_SUB_GROUP_SIZE(SUBGROUP_SIZE)
KERNEL(sdpa_opt_finalization_stage)(
 OPTIONAL_SHAPE_INFO_ARG
 __global OUTPUT_TYPE* output,
 const __global SOFTMAX_ACCUMULATOR_TYPE* exp_sums,
 const __global SOFTMAX_ACCUMULATOR_TYPE* max_logits,
 const __global OUTPUT_TYPE* tmp_out,
 const uint num_of_partitions) {
 const uint batch_idx = get_global_id(0);
 const uint b0_idx = batch_idx / NUM_HEADS;
 const uint b1_idx = batch_idx % NUM_HEADS;
 const uint target_seq_idx = get_global_id(1);
 const uint sglid = get_sub_group_local_id();
 if (num_of_partitions <= SUBGROUP_SIZE * REG_VERSION_MAX_VALUES_PER_WI_LOWER) {
 SOFTMAX_ACCUMULATOR_TYPE exp_sum[REG_VERSION_MAX_VALUES_PER_WI_LOWER] = {SOFTMAX_ACCUMULATOR_VAL_ZERO};
 SOFTMAX_ACCUMULATOR_TYPE max_logit[REG_VERSION_MAX_VALUES_PER_WI_LOWER] = {SOFTMAX_ACCUMULATOR_VAL_MIN};
 SOFTMAX_ACCUMULATOR_TYPE local_exp_sum = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 SOFTMAX_ACCUMULATOR_TYPE local_max_logit = SOFTMAX_ACCUMULATOR_VAL_MIN;
 const uint iters_num = CEIL_DIV(num_of_partitions, SUBGROUP_SIZE);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 const uint max_logit_offset = exp_sums_offset;
 if (partition_idx < num_of_partitions) {
 exp_sum[i] = exp_sums[exp_sums_offset];
 max_logit[i] = max_logits[max_logit_offset];
 local_max_logit = SOFTMAX_ACCUMULATOR_MAX_FUNC(local_max_logit, max_logit[i]);
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_max = sub_group_reduce_max(local_max_logit);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 if (partition_idx < num_of_partitions) {
 exp_sum[i] = exp_sum[i] * native_exp(max_logit[i] - global_max);
 local_exp_sum += exp_sum[i];
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_sum = sub_group_reduce_add(local_exp_sum);
 for (uint head_size_idx = 0; head_size_idx < HEAD_SIZE / SUBGROUP_SIZE; head_size_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE acc = 0.0f;
 for (uint partition_idx = 0; partition_idx < num_of_partitions; partition_idx++) {
 const uint tmp_out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 target_seq_idx * (num_of_partitions * HEAD_SIZE) +
 partition_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 OUTPUT_TYPE out_val = tmp_out[tmp_out_offset];
 acc += TO_SOFTMAX_ACCUMULATOR_TYPE(out_val) *
 TO_SOFTMAX_ACCUMULATOR_TYPE(sub_group_broadcast(exp_sum[partition_idx / SUBGROUP_SIZE], partition_idx % SUBGROUP_SIZE)) /
 TO_SOFTMAX_ACCUMULATOR_TYPE(global_sum);
 }
 const uint out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * HEAD_SIZE) +
 target_seq_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 output[out_offset] = TO_OUTPUT_TYPE(acc);
 }
 } else if (num_of_partitions <= SUBGROUP_SIZE * REG_VERSION_MAX_VALUES_PER_WI) {
 SOFTMAX_ACCUMULATOR_TYPE exp_sum[REG_VERSION_MAX_VALUES_PER_WI] = {SOFTMAX_ACCUMULATOR_VAL_ZERO};
 SOFTMAX_ACCUMULATOR_TYPE max_logit[REG_VERSION_MAX_VALUES_PER_WI] = {SOFTMAX_ACCUMULATOR_VAL_MIN};
 SOFTMAX_ACCUMULATOR_TYPE local_exp_sum = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 SOFTMAX_ACCUMULATOR_TYPE local_max_logit = SOFTMAX_ACCUMULATOR_VAL_MIN;
 const uint iters_num = CEIL_DIV(num_of_partitions, SUBGROUP_SIZE);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 const uint max_logit_offset = exp_sums_offset;
 if (partition_idx < num_of_partitions) {
 exp_sum[i] = exp_sums[exp_sums_offset];
 max_logit[i] = max_logits[max_logit_offset];
 local_max_logit = SOFTMAX_ACCUMULATOR_MAX_FUNC(local_max_logit, max_logit[i]);
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_max = sub_group_reduce_max(local_max_logit);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 if (partition_idx < num_of_partitions) {
 exp_sum[i] = exp_sum[i] * native_exp(max_logit[i] - global_max);
 local_exp_sum += exp_sum[i];
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_sum = sub_group_reduce_add(local_exp_sum);
 for (uint head_size_idx = 0; head_size_idx < HEAD_SIZE / SUBGROUP_SIZE; head_size_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE acc = 0.0f;
 for (uint partition_idx = 0; partition_idx < num_of_partitions; partition_idx++) {
 const uint tmp_out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 target_seq_idx * (num_of_partitions * HEAD_SIZE) +
 partition_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 OUTPUT_TYPE out_val = tmp_out[tmp_out_offset];
 acc += TO_SOFTMAX_ACCUMULATOR_TYPE(out_val) *
 TO_SOFTMAX_ACCUMULATOR_TYPE(sub_group_broadcast(exp_sum[partition_idx / SUBGROUP_SIZE], partition_idx % SUBGROUP_SIZE)) /
 TO_SOFTMAX_ACCUMULATOR_TYPE(global_sum);
 }
 const uint out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * HEAD_SIZE) +
 target_seq_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 output[out_offset] = TO_OUTPUT_TYPE(acc);
 }
 } else {
 SOFTMAX_ACCUMULATOR_TYPE local_exp_sum = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 SOFTMAX_ACCUMULATOR_TYPE local_max_logit = SOFTMAX_ACCUMULATOR_VAL_MIN;
 const uint iters_num = CEIL_DIV(num_of_partitions, SUBGROUP_SIZE);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 const uint max_logit_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 if (partition_idx < num_of_partitions) {
 local_max_logit = SOFTMAX_ACCUMULATOR_MAX_FUNC(local_max_logit, max_logits[max_logit_offset]);
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_max = sub_group_reduce_max(local_max_logit);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 const uint max_logit_offset = exp_sums_offset;
 if (partition_idx < num_of_partitions) {
 local_exp_sum += exp_sums[exp_sums_offset] * native_exp(max_logits[max_logit_offset] - global_max);
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_sum = sub_group_reduce_add(local_exp_sum);
 for (uint head_size_idx = 0; head_size_idx < HEAD_SIZE / SUBGROUP_SIZE; head_size_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE acc = 0.0f;
 for (uint partition_idx = 0; partition_idx < num_of_partitions; partition_idx++) {
 const uint tmp_out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 target_seq_idx * (num_of_partitions * HEAD_SIZE) +
 partition_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 const uint max_logit_offset = exp_sums_offset;
 SOFTMAX_ACCUMULATOR_TYPE new_exp_sum = exp_sums[exp_sums_offset] * native_exp(max_logits[max_logit_offset] - global_max);
 OUTPUT_TYPE out_val = tmp_out[tmp_out_offset];
 acc += TO_SOFTMAX_ACCUMULATOR_TYPE(out_val) * new_exp_sum / TO_SOFTMAX_ACCUMULATOR_TYPE(global_sum);
 }
 const uint out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * HEAD_SIZE) +
 target_seq_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 output[out_offset] = TO_OUTPUT_TYPE(acc);
 }
 }
}
#endif
#ifdef OUTPUT_BLOCK_READ
#undef OUTPUT_BLOCK_READ
#endif
#ifdef OUTPUT_BLOCK_WRITE
#undef OUTPUT_BLOCK_WRITE
#endif
#ifdef VALUE_BLOCK_READ
#undef VALUE_BLOCK_READ
#endif
#ifdef SUBGROUPS_PER_WG
#undef SUBGROUPS_PER_WG
#endif
#ifdef GET_COMPRESSION_INDEX
#undef GET_COMPRESSION_INDEX
#endif
#ifdef GET_COMPRESSION_INDEX
#undef GET_COMPRESSION_INDEX
#endif
#ifdef QUERY_STEP_LOCAL
#undef QUERY_STEP_LOCAL
#endif
#ifdef QUERY_BLOCK_SIZE
#undef QUERY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_SIZE
#undef KEY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef KEY_BLOCK
#undef KEY_BLOCK
#endif
#ifdef KEY_BLOCK_UNCOMPRESSED
#undef KEY_BLOCK_UNCOMPRESSED
#endif
#ifdef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#undef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#endif
#ifdef QUERY_BLOCK
#undef QUERY_BLOCK
#endif
#ifdef KEY_BLOCK_SIZE
#undef KEY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef KEY_BLOCK
#undef KEY_BLOCK
#endif
#ifdef KEY_BLOCK_UNCOMPRESSED
#undef KEY_BLOCK_UNCOMPRESSED
#endif
#ifdef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#undef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#endif
#ifdef QUERY_BLOCK
#undef QUERY_BLOCK
#endif
#ifdef KEY_BLOCK_SIZE
#undef KEY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef KEY_BLOCK
#undef KEY_BLOCK
#endif
#ifdef KEY_BLOCK_UNCOMPRESSED
#undef KEY_BLOCK_UNCOMPRESSED
#endif
#ifdef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#undef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#endif
#ifdef QUERY_BLOCK
#undef QUERY_BLOCK
#endif
#ifdef KEY_BLOCK_SIZE
#undef KEY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef KEY_BLOCK
#undef KEY_BLOCK
#endif
#ifdef KEY_BLOCK_UNCOMPRESSED
#undef KEY_BLOCK_UNCOMPRESSED
#endif
#ifdef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#undef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#endif
#ifdef QUERY_BLOCK
#undef QUERY_BLOCK
#endif
#ifdef SOURCE_SEQ_LEN
#undef SOURCE_SEQ_LEN
#endif
#ifdef TARGET_SEQ_LEN
#undef TARGET_SEQ_LEN
#endif
#ifdef PA_BUFFERS
#undef PA_BUFFERS
#endif
#ifdef PA_BUFFERS_ARGS
#undef PA_BUFFERS_ARGS
#endif
#ifdef PA_BUFFERS
#undef PA_BUFFERS
#endif
#ifdef PA_BUFFERS_ARGS
#undef PA_BUFFERS_ARGS
#endif
#ifdef ATTN_MASK_BUFFER
#undef ATTN_MASK_BUFFER
#endif
#ifdef ATTN_MASK_BUFFER_ARG
#undef ATTN_MASK_BUFFER_ARG
#endif
#ifdef ATTN_MASK_BUFFER
#undef ATTN_MASK_BUFFER
#endif
#ifdef ATTN_MASK_BUFFER_ARG
#undef ATTN_MASK_BUFFER_ARG
#endif
#ifdef ATTN_SCALE_BUFFER
#undef ATTN_SCALE_BUFFER
#endif
#ifdef ATTN_SCALE_BUFFER_ARG
#undef ATTN_SCALE_BUFFER_ARG
#endif
#ifdef ATTN_SCALE_BUFFER
#undef ATTN_SCALE_BUFFER
#endif
#ifdef ATTN_SCALE_BUFFER_ARG
#undef ATTN_SCALE_BUFFER_ARG
#endif
#ifdef APPLY_SCALES_TO_QUERY
#undef APPLY_SCALES_TO_QUERY
#endif
#ifdef MASK_VECTOR_TYPE
#undef MASK_VECTOR_TYPE
#endif
#ifdef ALIBI_TYPE
#undef ALIBI_TYPE
#endif
#ifdef ALIBI_TYPE
#undef ALIBI_TYPE
#endif
#ifdef batch_idx
#undef batch_idx
#endif
#ifdef num_heads_dim
#undef num_heads_dim
#endif
#ifdef b0_idx
#undef b0_idx
#endif
#ifdef b1_idx
#undef b1_idx
#endif
#ifdef target_seq_dim
#undef target_seq_dim
#endif
#ifdef target_seq_idx
#undef target_seq_idx
#endif
#ifdef head_size_idx
#undef head_size_idx
#endif
#ifdef sglid
#undef sglid
#endif
#ifdef sgid
#undef sgid
#endif
#ifdef KEY_SEQ_OFFSET
#undef KEY_SEQ_OFFSET
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef QUERY_VEC
#undef QUERY_VEC
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef QUERY_VEC_TYPE
#undef QUERY_VEC_TYPE
#endif
#ifdef KEY_UNPACKED_TYPE
#undef KEY_UNPACKED_TYPE
#endif
#ifdef KEY_UNPACKED_VEC_TYPE
#undef KEY_UNPACKED_VEC_TYPE
#endif
#ifdef TO_KEY_UNPACKED_TYPE
#undef TO_KEY_UNPACKED_TYPE
#endif
#ifdef KEY_UNPACKED_TYPE
#undef KEY_UNPACKED_TYPE
#endif
#ifdef KEY_UNPACKED_VEC_TYPE
#undef KEY_UNPACKED_VEC_TYPE
#endif
#ifdef TO_KEY_UNPACKED_TYPE
#undef TO_KEY_UNPACKED_TYPE
#endif
#ifdef key_vals
#undef key_vals
#endif
#ifdef REG_VERSION_MAX_VALUES_PER_WI
#undef REG_VERSION_MAX_VALUES_PER_WI
#endif
#ifdef REG_VERSION_MAX_VALUES_PER_WI_LOWER
#undef REG_VERSION_MAX_VALUES_PER_WI_LOWER
#endif
#ifdef REG_VERSION_MAX_VALUES_PER_WI
#undef REG_VERSION_MAX_VALUES_PER_WI
#endif
#ifdef REG_VERSION_MAX_VALUES_PER_WI_LOWER
#undef REG_VERSION_MAX_VALUES_PER_WI_LOWER
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
#ifdef BROADCAST_GROUP_SIZE
#undef BROADCAST_GROUP_SIZE
#endif
#ifdef DO_BROADCAST_KEY_VALUE
#undef DO_BROADCAST_KEY_VALUE
#endif
#ifdef IS_CAUSAL
#undef IS_CAUSAL
#endif
#ifdef HAS_ATTN_MASK_INPUT
#undef HAS_ATTN_MASK_INPUT
#endif
#ifdef HAS_SCALE_INPUT
#undef HAS_SCALE_INPUT
#endif
#ifdef IS_KV_COMPRESSED
#undef IS_KV_COMPRESSED
#endif
#ifdef INPUT0_DIMS_ORDER
#undef INPUT0_DIMS_ORDER
#endif
#ifdef INPUT1_DIMS_ORDER
#undef INPUT1_DIMS_ORDER
#endif
#ifdef INPUT2_DIMS_ORDER
#undef INPUT2_DIMS_ORDER
#endif
#ifdef TARGET_SEQ_LEN
#undef TARGET_SEQ_LEN
#endif
#ifdef NUM_HEADS
#undef NUM_HEADS
#endif
#ifdef NUM_KV_HEADS
#undef NUM_KV_HEADS
#endif
#ifdef SOURCE_SEQ_LEN
#undef SOURCE_SEQ_LEN
#endif
#ifdef SOFTMAX_ACCUMULATOR_TYPE
#undef SOFTMAX_ACCUMULATOR_TYPE
#endif
#ifdef SOFTMAX_ACCUMULATOR_VAL_MAX
#undef SOFTMAX_ACCUMULATOR_VAL_MAX
#endif
#ifdef SOFTMAX_ACCUMULATOR_VAL_MIN
#undef SOFTMAX_ACCUMULATOR_VAL_MIN
#endif
#ifdef SOFTMAX_ACCUMULATOR_VAL_ONE
#undef SOFTMAX_ACCUMULATOR_VAL_ONE
#endif
#ifdef SOFTMAX_ACCUMULATOR_VAL_ZERO
#undef SOFTMAX_ACCUMULATOR_VAL_ZERO
#endif
#ifdef TO_SOFTMAX_ACCUMULATOR_TYPE
#undef TO_SOFTMAX_ACCUMULATOR_TYPE
#endif
#ifdef TO_SOFTMAX_ACCUMULATOR_TYPE_SAT
#undef TO_SOFTMAX_ACCUMULATOR_TYPE_SAT
#endif
#ifdef AS_SOFTMAX_ACCUMULATOR_TYPE
#undef AS_SOFTMAX_ACCUMULATOR_TYPE
#endif
#ifdef SOFTMAX_ACCUMULATOR_MAX_FUNC
#undef SOFTMAX_ACCUMULATOR_MAX_FUNC
#endif
#ifdef SOFTMAX_ACCUMULATOR_MIN_FUNC
#undef SOFTMAX_ACCUMULATOR_MIN_FUNC
#endif
#ifdef SOFTMAX_ACCUMULATOR_ABS_FUNC
#undef SOFTMAX_ACCUMULATOR_ABS_FUNC
#endif
#ifdef SOFTMAX_ACCUMULATOR_TYPE_SIZE
#undef SOFTMAX_ACCUMULATOR_TYPE_SIZE
#endif
#ifdef SOFTMAX_ACCUMULATOR_IS_FP
#undef SOFTMAX_ACCUMULATOR_IS_FP
#endif
#ifdef SUBGROUP_SIZE
#undef SUBGROUP_SIZE
#endif
#ifdef HEAD_SIZE
#undef HEAD_SIZE
#endif
#ifdef SEQ_LEN_PARTITION_SIZE
#undef SEQ_LEN_PARTITION_SIZE
#endif
#ifdef TARGET_SEQ_LEN_BLOCK_SIZE
#undef TARGET_SEQ_LEN_BLOCK_SIZE
#endif
#ifdef SDPA_STAGE_0
#undef SDPA_STAGE_0
#endif
#ifdef SG_SCALE_FACTOR
#undef SG_SCALE_FACTOR
#endif
#ifdef STATIC_SCALE_VALUE_INV
#undef STATIC_SCALE_VALUE_INV
#endif
#ifdef STATIC_SCALE_VALUE
#undef STATIC_SCALE_VALUE
#endif

//====================================================
// Kernel template: sdpa_opt_multi_tokens 
// Kernel name: sdpa_opt_multi_tokens_8660372428234100028_0_0__sa
#define KERNEL(name) __kernel void sdpa_opt_multi_tokens_8660372428234100028_0_0__sa
#define KERNEL_ID sdpa_opt_multi_tokens_8660372428234100028_0_0__sa
#define FUNC(name)  _##name##_sdpa_opt_multi_tokens_8660372428234100028_0_0__sa
#define FUNC_CALL(name)  _##name##_sdpa_opt_multi_tokens_8660372428234100028_0_0__sa
#define CONST_ARRAY_DECL(name) __constant size_t  _##name##_sdpa_opt_multi_tokens_8660372428234100028_0_0__sa []
#define CONST_ARRAY_REF(name)  _##name##_sdpa_opt_multi_tokens_8660372428234100028_0_0__sa
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
#define BROADCAST_GROUP_SIZE 7
#define DO_BROADCAST_KEY_VALUE f /= 7;
#define IS_CAUSAL 0
#define HAS_ATTN_MASK_INPUT 1
#define HAS_SCALE_INPUT 0
#define IS_KV_COMPRESSED 0
#define INPUT0_DIMS_ORDER b,f,w,z,y,x
#define INPUT1_DIMS_ORDER b,f,w,z,y,x
#define INPUT2_DIMS_ORDER b,f,w,z,y,x
#define TARGET_SEQ_LEN (shape_info[6] )
#define NUM_HEADS 28
#define NUM_KV_HEADS -1
#define SOURCE_SEQ_LEN (shape_info[14] )
#define SOFTMAX_ACCUMULATOR_TYPE float
#define SOFTMAX_ACCUMULATOR_VAL_MAX FLT_MAX
#define SOFTMAX_ACCUMULATOR_VAL_MIN -SOFTMAX_ACCUMULATOR_VAL_MAX
#define SOFTMAX_ACCUMULATOR_VAL_ONE 1.0f
#define SOFTMAX_ACCUMULATOR_VAL_ZERO 0.0f
#define TO_SOFTMAX_ACCUMULATOR_TYPE(v) convert_float(v)
#define TO_SOFTMAX_ACCUMULATOR_TYPE_SAT(v) convert_float(v)
#define AS_SOFTMAX_ACCUMULATOR_TYPE(v) as_float(v)
#define SOFTMAX_ACCUMULATOR_MAX_FUNC fmax
#define SOFTMAX_ACCUMULATOR_MIN_FUNC fmin
#define SOFTMAX_ACCUMULATOR_ABS_FUNC fabs
#define SOFTMAX_ACCUMULATOR_TYPE_SIZE 4
#define SOFTMAX_ACCUMULATOR_IS_FP 1
#define SUBGROUP_SIZE 16
#define HEAD_SIZE 128
#define SEQ_LEN_PARTITION_SIZE 256
#define TARGET_SEQ_LEN_BLOCK_SIZE 16
#define SDPA_STAGE_0 1
#define SG_SCALE_FACTOR 2
#define STATIC_SCALE_VALUE_INV as_float(0x413504f3)/*1.131371e+01*/
#define STATIC_SCALE_VALUE as_float(0x3db504f2)/*8.838834e-02*/


inline uint FUNC(get_input0_index_nt)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#if INPUT0_SIMPLE
 return GET_DATA_INDEX_6D(INPUT0, b, f, w, z, y, x);
#else
#if INPUT0_DIMS == 4
 return INPUT0_GET_INDEX(b, f, y, x);
#elif INPUT0_DIMS == 5
 return INPUT0_GET_INDEX(b, f, z, y, x);
#elif INPUT0_DIMS == 6
 return INPUT0_GET_INDEX(b, f, w, z, y, x);
#else
# error sdpa_opt.cl : Unsupported input 0 format
#endif
#endif
}
inline uint FUNC(get_input0_index)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef INPUT0_DIMS_ORDER
 return FUNC_CALL(get_input0_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT0_DIMS_ORDER);
#else
 return FUNC_CALL(get_input0_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR b, f, w, z, y, x);
#endif
}
inline uint FUNC(get_input1_index_nt)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef DO_BROADCAST_KEY_VALUE
 DO_BROADCAST_KEY_VALUE;
#endif
#if INPUT1_SIMPLE
 return GET_DATA_INDEX_6D(INPUT1, b, f, w, z, y, x);
#else
#if INPUT1_DIMS == 4
 return INPUT1_GET_INDEX(b, f, y, x);
#elif INPUT1_DIMS == 5
 return INPUT1_GET_INDEX(b, f, z, y, x);
#elif INPUT1_DIMS == 6
 return INPUT1_GET_INDEX(b, f, w, z, y, x);
#else
# error sdpa_opt.cl : Unsupported input 1 format
#endif
#endif
}
inline uint FUNC(get_input1_index)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef INPUT1_DIMS_ORDER
 return FUNC_CALL(get_input1_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT1_DIMS_ORDER);
#else
 return FUNC_CALL(get_input1_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR b, f, w, z, y, x);
#endif
}
inline uint FUNC(get_input2_index_nt)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef DO_BROADCAST_KEY_VALUE
 DO_BROADCAST_KEY_VALUE;
#endif
#if INPUT2_SIMPLE
 return GET_DATA_INDEX_6D_SAFE(INPUT2, b, f, w, z, y, x);
#else
#if INPUT2_DIMS == 4
 return INPUT2_GET_INDEX(b, f, y, x);
#elif INPUT2_DIMS == 5
 return INPUT2_GET_INDEX(b, f, z, y, x);
#elif INPUT2_DIMS == 6
 return INPUT2_GET_INDEX(b, f, w, z, y, x);
#else
# error sdpa_opt.cl : Unsupported input 1 format
#endif
#endif
}
inline uint FUNC(get_input2_index)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#ifdef INPUT2_DIMS_ORDER
 return FUNC_CALL(get_input2_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT2_DIMS_ORDER);
#else
 return FUNC_CALL(get_input2_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR b, f, w, z, y, x);
#endif
}
#ifdef BEAM_TABLE_TYPE
inline uint FUNC(get_bt_index_nt)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
#if BEAM_TABLE_SIMPLE
 return GET_DATA_INDEX_6D_SAFE(BEAM_TABLE, b, f, w, z, y, x);
#else
# error sdpa_opt.cl : Unsupported beam table format
#endif
}
inline uint FUNC(get_bt_index_key)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
 return FUNC_CALL(get_bt_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT1_DIMS_ORDER);
}
inline uint FUNC(get_bt_index_value)(OPTIONAL_SHAPE_INFO_ARG uint b, uint f, uint w, uint z, uint y, uint x) {
 return FUNC_CALL(get_bt_index_nt)(OPTIONAL_SHAPE_INFO_TENSOR INPUT2_DIMS_ORDER);
}
#endif
#define OUTPUT_BLOCK_READ(ptr, offset) BLOCK_READN(OUTPUT_TYPE, 1, ptr, offset)
#define OUTPUT_BLOCK_WRITE(ptr, offset, val) BLOCK_WRITEN(OUTPUT_TYPE, 1, ptr, offset, val)
#define VALUE_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT2_TYPE, 1, ptr, offset)
#define SUBGROUPS_PER_WG (HEAD_SIZE * SG_SCALE_FACTOR / SUBGROUP_SIZE)
#if IS_KV_COMPRESSED
#if COMPRESSED_PER_HEAD
 #define GET_COMPRESSION_INDEX(INPUT, b, f, y, x) GET_DATA_INDEX(INPUT, (b), (f), (y), (0));
#else
 #define GET_COMPRESSION_INDEX(INPUT, b, f, y, x) GET_DATA_INDEX(INPUT, (b), (0), (y), (0));
#endif
#endif
#ifdef SDPA_STAGE_0
#if TARGET_SEQ_LEN_BLOCK_SIZE == 1
REQD_SUB_GROUP_SIZE(SUBGROUP_SIZE)
__attribute__((reqd_work_group_size(1, 1, HEAD_SIZE * SG_SCALE_FACTOR)))
KERNEL(sdpa_opt)(
 OPTIONAL_SHAPE_INFO_ARG
 const __global INPUT0_TYPE* query_input,
 const __global INPUT1_TYPE* key_input,
 const __global INPUT2_TYPE* value_input,
#if HAS_ATTN_MASK_INPUT
 const __global INPUT3_TYPE* attn_mask,
#endif
#if HAS_SCALE_INPUT
 const __global INPUT4_TYPE* scale,
#endif
 __global OUTPUT_TYPE* output,
#if IS_KV_COMPRESSED
 const __global KEY_COMPRESSION_SCALE_TYPE* key_scale,
 const __global VALUE_COMPRESSION_SCALE_TYPE* val_scale,
#endif
#ifdef BEAM_TABLE_TYPE
 const __global BEAM_TABLE_TYPE* beam_table,
#endif
 __global SOFTMAX_ACCUMULATOR_TYPE* exp_sums,
 __global SOFTMAX_ACCUMULATOR_TYPE* max_logits,
 __global OUTPUT_TYPE* tmp_out
)
{
 const uint batch_idx = get_global_id(0);
 const uint b0_idx = batch_idx / NUM_HEADS;
 const uint b1_idx = batch_idx % NUM_HEADS;
 const uint target_seq_idx = get_global_id(1);
 const uint lid = get_local_id(2);
#if SG_SCALE_FACTOR == 2
 const uint head_size_idx = lid % HEAD_SIZE;
#elif SG_SCALE_FACTOR == 1
 const uint head_size_idx = lid;
#else
 #error "sdpa_opt.cl: Unsupported scale factor"
#endif
#if SUBGROUPS_PER_WG > SUBGROUP_SIZE
 #error "sdpa_opt.cl: Number of subgroups per work group should be less than subgroup_size
#endif
 const uint sgid = get_sub_group_id();
 const uint sglid = get_sub_group_local_id();
 const uint partition_idx = get_group_id(2);
 const uint num_of_partitions = get_num_groups(2);
 const uint wi_num_per_partition = get_local_size(2);
 const uint start_partition_idx = partition_idx * SEQ_LEN_PARTITION_SIZE;
 const uint partition_seq_len =
 ((partition_idx + 1) < num_of_partitions) ? (SEQ_LEN_PARTITION_SIZE)
 : (SOURCE_SEQ_LEN - partition_idx * SEQ_LEN_PARTITION_SIZE);
 __local INPUT0_TYPE query_local[HEAD_SIZE * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE qk_local[SEQ_LEN_PARTITION_SIZE * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE qk_max_vals[SUBGROUPS_PER_WG * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE qk_sum_vals[SUBGROUPS_PER_WG * TARGET_SEQ_LEN_BLOCK_SIZE];
 {
 SOFTMAX_ACCUMULATOR_TYPE qk_max[TARGET_SEQ_LEN_BLOCK_SIZE] = {SOFTMAX_ACCUMULATOR_VAL_MIN};
 for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
 qk_max[i] = SOFTMAX_ACCUMULATOR_VAL_MIN;
 }
 {
#if HAS_SCALE_INPUT
 const OUTPUT_TYPE scale_val = *scale;
#else
 const OUTPUT_TYPE scale_val = OUTPUT_VAL_ONE / sqrt(TO_OUTPUT_TYPE(HEAD_SIZE));
#endif
 {
 #define QUERY_STEP_LOCAL SUBGROUP_SIZE * SUBGROUPS_PER_WG
 uint query_local_offset = sgid * SUBGROUP_SIZE + sglid;
 const uint seq_idx_end = 1;
#ifdef INPUT0_DIMS_ORDER
 uint query_offset = FUNC_CALL(get_input0_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, target_seq_idx, (sgid * SUBGROUP_SIZE));
 uint query_offset_next_seq = FUNC_CALL(get_input0_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, target_seq_idx + 1, (sgid * SUBGROUP_SIZE));
 const uint query_pitch = query_offset_next_seq - query_offset;
#else
 uint query_offset = INPUT0_GET_INDEX(b0_idx, b1_idx, target_seq_idx, (sgid * SUBGROUP_SIZE));
 const uint query_pitch = QUERY_STEP_LOCAL;
#endif
#if SG_SCALE_FACTOR == 2
 if (sgid < HEAD_SIZE / SUBGROUP_SIZE) {
#else
 {
#endif
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 #define QUERY_BLOCK_SIZE 1
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, QUERY_BLOCK_SIZE, query_input, query_offset);
 query_local[query_local_offset] = val * scale_val;
 query_local_offset += QUERY_STEP_LOCAL;
 query_offset += query_pitch;
 }
 }
 #undef QUERY_BLOCK_SIZE
 #undef QUERY_STEP
 barrier(CLK_LOCAL_MEM_FENCE);
 }
 for (uint seq_len = sgid; seq_len < partition_seq_len; seq_len += (HEAD_SIZE / SUBGROUP_SIZE) * SG_SCALE_FACTOR) {
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_key)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + seq_len, 0)];
#else
 const uint b_idx = b0_idx;
#endif
#ifdef INPUT1_DIMS_ORDER
 uint key_offset = FUNC_CALL(get_input1_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + seq_len, 0);
#else
 uint key_offset = INPUT1_GET_INDEX(b_idx, b1_idx, start_partition_idx + seq_len, 0);
#endif
 SOFTMAX_ACCUMULATOR_TYPE acc[TARGET_SEQ_LEN_BLOCK_SIZE] = {SOFTMAX_ACCUMULATOR_VAL_ZERO};
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(KEY_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + seq_len, 0);
 KEY_COMPRESSION_SCALE_TYPE comp_scale = key_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 KEY_COMPRESSION_SCALE_TYPE comp_zp = key_scale[comp_offset + 1];
#endif
#endif
 uint head_idx_index = 0;
 #define KEY_BLOCK_SIZE 8
 for (; head_idx_index + (KEY_BLOCK_SIZE * SUBGROUP_SIZE) <= HEAD_SIZE; head_idx_index += SUBGROUP_SIZE * KEY_BLOCK_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, KEY_BLOCK_SIZE, ptr, offset);
 #define KEY_BLOCK MAKE_VECTOR_TYPE(INPUT1_TYPE, KEY_BLOCK_SIZE)
 #define KEY_BLOCK_UNCOMPRESSED MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, KEY_BLOCK_SIZE)
 #define TO_KEY_BLOCK_UNCOMPRESSED_TYPE(val) CAT(convert_, KEY_BLOCK_UNCOMPRESSED)(val)
 #define QUERY_BLOCK MAKE_VECTOR_TYPE(INPUT0_TYPE, KEY_BLOCK_SIZE)
 KEY_BLOCK key_vec_packed = KEY_BLOCK_READ(key_input, key_offset + head_idx_index);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed) - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed)) * comp_scale;
#else
 KEY_BLOCK key_vals = key_vec_packed;
#endif
 uint query_offset = head_idx_index + sglid;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 QUERY_BLOCK query_vals_reg;
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 query_vals_reg[i] = query_local[query_offset + i * SUBGROUP_SIZE];
 }
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 acc[seq_idx] = mad(TO_SOFTMAX_ACCUMULATOR_TYPE(query_vals_reg[i]), TO_SOFTMAX_ACCUMULATOR_TYPE(key_vals[i]), acc[seq_idx]);
 }
 query_offset += HEAD_SIZE;
 }
 }
 #define KEY_BLOCK_SIZE 4
 for (; head_idx_index + (KEY_BLOCK_SIZE * SUBGROUP_SIZE) <= HEAD_SIZE; head_idx_index += SUBGROUP_SIZE * KEY_BLOCK_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, KEY_BLOCK_SIZE, ptr, offset);
 #define KEY_BLOCK MAKE_VECTOR_TYPE(INPUT1_TYPE, KEY_BLOCK_SIZE)
 #define KEY_BLOCK_UNCOMPRESSED MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, KEY_BLOCK_SIZE)
 #define TO_KEY_BLOCK_UNCOMPRESSED_TYPE(val) CAT(convert_, KEY_BLOCK_UNCOMPRESSED)(val)
 #define QUERY_BLOCK MAKE_VECTOR_TYPE(INPUT0_TYPE, KEY_BLOCK_SIZE)
 KEY_BLOCK key_vec_packed = KEY_BLOCK_READ(key_input, key_offset + head_idx_index);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed) - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed)) * comp_scale;
#else
 KEY_BLOCK key_vals = key_vec_packed;
#endif
 uint query_offset = head_idx_index + sglid;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 QUERY_BLOCK query_vals_reg;
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 query_vals_reg[i] = query_local[query_offset + i * SUBGROUP_SIZE];
 }
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 acc[seq_idx] = mad(TO_SOFTMAX_ACCUMULATOR_TYPE(query_vals_reg[i]), TO_SOFTMAX_ACCUMULATOR_TYPE(key_vals[i]), acc[seq_idx]);
 }
 query_offset += HEAD_SIZE;
 }
 }
 #define KEY_BLOCK_SIZE 2
 for (; head_idx_index + (KEY_BLOCK_SIZE * SUBGROUP_SIZE) <= HEAD_SIZE; head_idx_index += SUBGROUP_SIZE * KEY_BLOCK_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, KEY_BLOCK_SIZE, ptr, offset);
 #define KEY_BLOCK MAKE_VECTOR_TYPE(INPUT1_TYPE, KEY_BLOCK_SIZE)
 #define KEY_BLOCK_UNCOMPRESSED MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, KEY_BLOCK_SIZE)
 #define TO_KEY_BLOCK_UNCOMPRESSED_TYPE(val) CAT(convert_, KEY_BLOCK_UNCOMPRESSED)(val)
 #define QUERY_BLOCK MAKE_VECTOR_TYPE(INPUT0_TYPE, KEY_BLOCK_SIZE)
 KEY_BLOCK key_vec_packed = KEY_BLOCK_READ(key_input, key_offset + head_idx_index);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed) - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed)) * comp_scale;
#else
 KEY_BLOCK key_vals = key_vec_packed;
#endif
 uint query_offset = head_idx_index + sglid;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 QUERY_BLOCK query_vals_reg;
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 query_vals_reg[i] = query_local[query_offset + i * SUBGROUP_SIZE];
 }
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 acc[seq_idx] = mad(TO_SOFTMAX_ACCUMULATOR_TYPE(query_vals_reg[i]), TO_SOFTMAX_ACCUMULATOR_TYPE(key_vals[i]), acc[seq_idx]);
 }
 query_offset += HEAD_SIZE;
 }
 }
 #define KEY_BLOCK_SIZE 1
 for (; head_idx_index + (KEY_BLOCK_SIZE * SUBGROUP_SIZE) <= HEAD_SIZE; head_idx_index += SUBGROUP_SIZE * KEY_BLOCK_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, KEY_BLOCK_SIZE, ptr, offset);
 #define KEY_BLOCK MAKE_VECTOR_TYPE(INPUT1_TYPE, KEY_BLOCK_SIZE)
 #define KEY_BLOCK_UNCOMPRESSED MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, KEY_BLOCK_SIZE)
 #define TO_KEY_BLOCK_UNCOMPRESSED_TYPE(val) CAT(convert_, KEY_BLOCK_UNCOMPRESSED)(val)
 #define QUERY_BLOCK MAKE_VECTOR_TYPE(INPUT0_TYPE, KEY_BLOCK_SIZE)
 KEY_BLOCK key_vec_packed = KEY_BLOCK_READ(key_input, key_offset + head_idx_index);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed) - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 KEY_BLOCK_UNCOMPRESSED key_vals = (TO_KEY_BLOCK_UNCOMPRESSED_TYPE(key_vec_packed)) * comp_scale;
#else
 KEY_BLOCK key_vals = key_vec_packed;
#endif
 uint query_offset = head_idx_index + sglid;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 QUERY_BLOCK query_vals_reg;
 unroll_for(uint i = 0; i < KEY_BLOCK_SIZE; i++) {
 query_vals_reg = query_local[query_offset + i * SUBGROUP_SIZE];
 }
 acc[seq_idx] = mad(TO_SOFTMAX_ACCUMULATOR_TYPE(query_vals_reg), TO_SOFTMAX_ACCUMULATOR_TYPE(key_vals), acc[seq_idx]);
 query_offset += HEAD_SIZE;
 }
 }
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc[seq_idx] = sub_group_reduce_add(acc[seq_idx]);
 qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len] = acc[seq_idx];
 }
 }
 {
 barrier(CLK_LOCAL_MEM_FENCE);
 SOFTMAX_ACCUMULATOR_TYPE qk_val[TARGET_SEQ_LEN_BLOCK_SIZE];
 const uint seq_idx_end = 1;
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 for (uint seq_len = sgid * SUBGROUP_SIZE + sglid; seq_len < partition_seq_len; seq_len += (HEAD_SIZE * SG_SCALE_FACTOR)) {
 qk_val[seq_idx] = qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len];
#if IS_CAUSAL
 if (start_partition_idx + seq_len > target_seq_idx + seq_idx)
 qk_val[seq_idx] += INPUT0_VAL_MIN;
#elif !IS_CAUSAL && HAS_ATTN_MASK_INPUT
 const uint attn_mask_offset = INPUT3_GET_INDEX_SAFE(b0_idx, b1_idx, target_seq_idx + seq_idx, start_partition_idx + seq_len);
 qk_val[seq_idx] += attn_mask[attn_mask_offset];
#endif
 qk_max[seq_idx] = SOFTMAX_ACCUMULATOR_MAX_FUNC(qk_max[seq_idx], TO_SOFTMAX_ACCUMULATOR_TYPE(qk_val[seq_idx]));
 qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len] = qk_val[seq_idx];
 }
 }
 }
 }
 {
 const uint seq_idx_end = 1;
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 qk_max[seq_idx] = sub_group_reduce_max(qk_max[seq_idx]);
 }
 if (sglid == 0) {
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 qk_max_vals[seq_idx * SUBGROUPS_PER_WG + sgid] = qk_max[seq_idx];
 }
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_max[seq_idx] = SOFTMAX_ACCUMULATOR_VAL_MIN;
 if (sglid < SUBGROUPS_PER_WG)
 qk_max[seq_idx] = qk_max_vals[seq_idx * SUBGROUPS_PER_WG + sglid];
 qk_max[seq_idx] = sub_group_reduce_max(qk_max[seq_idx]);
 }
 SOFTMAX_ACCUMULATOR_TYPE exp_sum[TARGET_SEQ_LEN_BLOCK_SIZE] = {SOFTMAX_ACCUMULATOR_VAL_ZERO};
 const uint qk_num_per_wi = CEIL_DIV(partition_seq_len, SUBGROUPS_PER_WG * SUBGROUP_SIZE);
 for (uint qk_idx = 0; qk_idx < qk_num_per_wi; qk_idx++) {
 const uint local_data_idx = qk_idx * (SUBGROUPS_PER_WG * SUBGROUP_SIZE) + lid;
 if (local_data_idx < partition_seq_len) {
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE qk_new = native_exp(TO_SOFTMAX_ACCUMULATOR_TYPE(qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + local_data_idx]) - qk_max[seq_idx]);
 qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + local_data_idx] = TO_OUTPUT_TYPE(qk_new);
 exp_sum[seq_idx] += qk_new;
 }
 }
 }
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 exp_sum[seq_idx] = sub_group_reduce_add(exp_sum[seq_idx]);
 if (sglid == 0)
 qk_sum_vals[seq_idx * SUBGROUPS_PER_WG + sgid] = exp_sum[seq_idx];
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 exp_sum[seq_idx] = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 if (sglid < SUBGROUPS_PER_WG)
 exp_sum[seq_idx] = qk_sum_vals[seq_idx * SUBGROUPS_PER_WG + sglid];
 exp_sum[seq_idx] = sub_group_reduce_add(exp_sum[seq_idx]);
 }
 for (uint qk_idx = 0; qk_idx < qk_num_per_wi; qk_idx++) {
 const uint local_data_idx = qk_idx * (SUBGROUPS_PER_WG * SUBGROUP_SIZE) + lid;
 if (local_data_idx < partition_seq_len) {
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE qk_new = TO_SOFTMAX_ACCUMULATOR_TYPE(qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + local_data_idx]) / exp_sum[seq_idx];
 qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + local_data_idx] = TO_OUTPUT_TYPE(qk_new);
 }
 }
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 {
 if (num_of_partitions > 1 && lid == 0) {
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 (seq_idx + target_seq_idx) * (num_of_partitions) +
 partition_idx;
 exp_sums[exp_sums_offset] = exp_sum[seq_idx];
 const uint max_logits_offset = exp_sums_offset;
 max_logits[max_logits_offset] = qk_max[seq_idx];
 }
 }
 }
 }
 }
 {
 OUTPUT_TYPE acc[TARGET_SEQ_LEN_BLOCK_SIZE] = {OUTPUT_VAL_ZERO};
#ifndef BEAM_TABLE_TYPE
#ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, 0, 0);
 uint value_offset_next_seq = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, 1, 0);
 const uint value_pitch = value_offset_next_seq - value_offset;
#else
 const uint value_pitch = HEAD_SIZE;
#endif
#endif
#if SG_SCALE_FACTOR > 1
 const uint seq_len_start = (sgid / (HEAD_SIZE / SUBGROUP_SIZE)) * (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR / SUBGROUP_SIZE);
 const uint seq_len_end = min(seq_len_start + (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR / SUBGROUP_SIZE), partition_seq_len / SUBGROUP_SIZE);
#else
 const uint seq_len_start = 0;
 const uint seq_len_end = partition_seq_len / SUBGROUP_SIZE;
#endif
 for (uint seq_len = seq_len_start; seq_len < seq_len_end; seq_len++) {
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, sgid * SUBGROUP_SIZE)];
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, sgid * SUBGROUP_SIZE);
#else
 const uint b_idx = b0_idx;
#ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE), head_size_idx);
#else
 uint value_offset = INPUT2_GET_INDEX(b_idx, b1_idx, start_partition_idx + (seq_len * SUBGROUP_SIZE), head_size_idx);
#endif
#endif
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 OUTPUT_TYPE qk_val[TARGET_SEQ_LEN_BLOCK_SIZE];
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len * SUBGROUP_SIZE + sglid];
 }
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, sub_group_broadcast(value_offset, i));
#else
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - sub_group_broadcast(comp_zp, i)) * sub_group_broadcast(comp_scale, i);
#elif IS_KV_COMPRESSED
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * sub_group_broadcast(comp_scale, i));
#else
 INPUT2_TYPE value_val = value_packed;
#endif
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc[seq_idx] = mad(sub_group_broadcast(qk_val[seq_idx], i), value_val, acc[seq_idx]);
 }
#ifndef BEAM_TABLE_TYPE
 value_offset += value_pitch;
#endif
 }
 }
#if SG_SCALE_FACTOR > 1
 if (sgid >= HEAD_SIZE / SUBGROUP_SIZE) {
#endif
 for (uint seq_len = (partition_seq_len / SUBGROUP_SIZE) * SUBGROUP_SIZE; seq_len < partition_seq_len; seq_len++) {
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + seq_len, head_size_idx)];
#else
 const uint b_idx = b0_idx;
#endif
#ifdef INPUT2_DIMS_ORDER
 const uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + seq_len, head_size_idx);
#else
 const uint value_offset = INPUT2_GET_INDEX(b_idx, b1_idx, start_partition_idx + seq_len, head_size_idx);
#endif
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + seq_len, 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 OUTPUT_TYPE qk_val[TARGET_SEQ_LEN_BLOCK_SIZE];
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = qk_local[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len];
 }
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 const VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - comp_zp) * comp_scale;
#elif IS_KV_COMPRESSED
 const VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * comp_scale);
#else
 const INPUT2_TYPE value_val = value_packed;
#endif
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc[seq_idx] = mad(qk_val[seq_idx], value_val, acc[seq_idx]);
 }
 }
#if SG_SCALE_FACTOR > 1
 }
#endif
#if SG_SCALE_FACTOR > 1
 if ((partition_seq_len > (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR)) || (partition_seq_len % SUBGROUP_SIZE != 0)) {
 if (sgid >= HEAD_SIZE / SUBGROUP_SIZE) {
 query_local[head_size_idx] = acc[0];
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 if (sgid < HEAD_SIZE / SUBGROUP_SIZE) {
 acc[0] += query_local[head_size_idx];
 }
 }
#endif
#if SG_SCALE_FACTOR > 1
 if (sgid < HEAD_SIZE / SUBGROUP_SIZE) {
#endif
 if (num_of_partitions > 1) {
 const uint seq_idx_end = 1;
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 const uint tmp_out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 (target_seq_idx + seq_idx) * (num_of_partitions * HEAD_SIZE) +
 partition_idx * (HEAD_SIZE) +
 head_size_idx;
 tmp_out[tmp_out_offset] = acc[seq_idx];
 }
 } else {
 const uint seq_idx_end = 1;
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 const uint output_offset = OUTPUT_GET_INDEX(b0_idx, b1_idx, target_seq_idx + seq_idx, head_size_idx);
 output[output_offset] = acc[seq_idx];
 }
 }
#if SG_SCALE_FACTOR > 1
 }
#endif
 }
}
#else
#if IS_PAGED_ATTENTION
 #define SOURCE_SEQ_LEN (subsequence_begins[gws_seq_indexes_correspondence[((uint)get_global_id(1))] + 1] - subsequence_begins[gws_seq_indexes_correspondence[((uint)get_global_id(1))]])
 #define TARGET_SEQ_LEN (subsequence_begins[gws_seq_indexes_correspondence[((uint)get_global_id(1))] + 1] - subsequence_begins[gws_seq_indexes_correspondence[((uint)get_global_id(1))]])
 #define PA_BUFFERS , subsequence_begins, blocked_indexes_start, blocked_indexes_end, gws_seq_indexes_correspondence
 #define PA_BUFFERS_ARGS , const __global INPUT3_TYPE* subsequence_begins, const __global int* blocked_indexes_start, const __global int* blocked_indexes_end, const __global int* gws_seq_indexes_correspondence
#else
 #define PA_BUFFERS
 #define PA_BUFFERS_ARGS
#endif
#if HAS_ATTN_MASK_INPUT
 #define ATTN_MASK_BUFFER , attn_mask
 #define ATTN_MASK_BUFFER_ARG , const __global INPUT3_TYPE* attn_mask
#else
 #define ATTN_MASK_BUFFER
 #define ATTN_MASK_BUFFER_ARG
#endif
#if HAS_SCALE_INPUT
 #define ATTN_SCALE_BUFFER , scale
 #define ATTN_SCALE_BUFFER_ARG , const __global INPUT4_TYPE* scale
#else
 #define ATTN_SCALE_BUFFER
 #define ATTN_SCALE_BUFFER_ARG
#endif
#if IS_KV_COMPRESSED
#define APPLY_SCALES_TO_QUERY 1
#endif
#define MASK_VECTOR_TYPE MAKE_VECTOR_TYPE(INPUT0_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
inline MASK_VECTOR_TYPE FUNC(load_attn_mask)(OPTIONAL_SHAPE_INFO_ARG
 uint b0_idx,
 uint b1_idx,
 uint target_seq_idx,
 uint source_seq_idx
 ATTN_MASK_BUFFER_ARG
 ATTN_SCALE_BUFFER_ARG
 PA_BUFFERS_ARGS
 ) {
 MASK_VECTOR_TYPE mask_vec = INPUT0_VAL_ZERO;
#if !IS_CAUSAL && HAS_ATTN_MASK_INPUT
 const uint attn_mask_offset = INPUT3_GET_INDEX_SAFE(b0_idx, b1_idx, target_seq_idx, source_seq_idx);
 if (target_seq_idx >= (uint)TARGET_SEQ_LEN) {
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 mask_vec[i] = NAN;
 }
 } else {
 if (source_seq_idx + SUBGROUP_SIZE <= (uint)SOURCE_SEQ_LEN) {
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 const INPUT3_TYPE mask_val = attn_mask[attn_mask_offset + i];
 mask_vec[i] = mask_val;
 }
 } else {
 const uint max_mask_offset = min(source_seq_idx + SUBGROUP_SIZE, (uint)SOURCE_SEQ_LEN);
 for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 const INPUT3_TYPE mask_val = source_seq_idx + i < max_mask_offset ? attn_mask[attn_mask_offset + i] : NAN;
 mask_vec[i] = mask_val;
 }
 }
 }
#endif
#if !IS_CAUSAL && !HAS_ATTN_MASK_INPUT
 if (target_seq_idx >= (uint)TARGET_SEQ_LEN) {
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 mask_vec[i] = NAN;
 }
 } else {
 const uint max_mask_offset = min(source_seq_idx + SUBGROUP_SIZE, (uint)SOURCE_SEQ_LEN);
 for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 mask_vec[i] = source_seq_idx + i < max_mask_offset ? 0 : NAN;
 }
 }
#endif
#if IS_CAUSAL
 if (target_seq_idx >= (uint)TARGET_SEQ_LEN) {
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 mask_vec[i] = NAN;
 }
 } else {
 for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 if (source_seq_idx + i > target_seq_idx)
 mask_vec[i] = NAN;
 }
 }
#endif
#if HAS_SCALE_INPUT
 const OUTPUT_TYPE scale_val = OUTPUT_VAL_ONE / *scale;
#else
 const INPUT0_TYPE scale_val = TO_INPUT0_TYPE(STATIC_SCALE_VALUE_INV);
#endif
#if IS_CAUSAL || HAS_ATTN_MASK_INPUT
 mask_vec *= scale_val;
#endif
 return mask_vec;
}
#if IS_PAGED_ATTENTION && HAS_ALIBI
#if HAS_SCALE_INPUT
#define ALIBI_TYPE INPUT5_TYPE
#else
#define ALIBI_TYPE INPUT4_TYPE
#endif
#endif
REQD_SUB_GROUP_SIZE(SUBGROUP_SIZE)
KERNEL(sdpa_opt)(
 OPTIONAL_SHAPE_INFO_ARG
 const __global INPUT0_TYPE* query_input,
 const __global INPUT1_TYPE* key_input,
 const __global INPUT2_TYPE* value_input,
#if IS_PAGED_ATTENTION
 const __global INPUT3_TYPE* subsequence_begins,
#endif
#if HAS_ATTN_MASK_INPUT
 const __global INPUT3_TYPE* attn_mask,
#endif
#if HAS_SCALE_INPUT
 const __global INPUT4_TYPE* scale,
#endif
#if IS_PAGED_ATTENTION && HAS_ALIBI
 const __global ALIBI_TYPE* alibi_slopes,
#endif
 __global OUTPUT_TYPE* output,
#if IS_KV_COMPRESSED
 const __global KEY_COMPRESSION_SCALE_TYPE* key_scale,
 const __global VALUE_COMPRESSION_SCALE_TYPE* val_scale,
#endif
#ifdef BEAM_TABLE_TYPE
 const __global BEAM_TABLE_TYPE* beam_table,
#endif
#if IS_PAGED_ATTENTION
 const __global int* blocked_indexes_start,
 const __global int* blocked_indexes_end,
 const __global int* gws_seq_indexes_correspondence
#else
 __global SOFTMAX_ACCUMULATOR_TYPE* exp_sums,
 __global SOFTMAX_ACCUMULATOR_TYPE* max_logits,
 __global OUTPUT_TYPE* tmp_out
#endif
)
{
#if TARGET_SEQ_LEN_BLOCK_SIZE != 16
 #error sdpa_opt.cl: unsupported TARGET_SEQ_LEN_BLOCK_SIZE
#endif
 #define batch_idx ((uint)get_global_id(0))
 #define num_heads_dim ((uint)get_global_id(0))
 #define b0_idx (batch_idx / NUM_HEADS)
 #define b1_idx (batch_idx % NUM_HEADS)
 #define target_seq_dim ((uint)get_global_id(1))
 #define target_seq_idx ((uint)get_global_id(1) * TARGET_SEQ_LEN_BLOCK_SIZE)
 #define head_size_idx ((uint)get_local_id(2) % HEAD_SIZE)
 #define sglid (uint)get_sub_group_local_id()
 #define sgid (uint)get_sub_group_id()
 __local INPUT0_TYPE slm_query[HEAD_SIZE * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local OUTPUT_TYPE slm_qk_vals[SEQ_LEN_PARTITION_SIZE * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_qk_max_vals[SUBGROUPS_PER_WG * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_exp_sum_vals[SUBGROUPS_PER_WG * TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_exp_sum_cur[TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_max_val_cur[TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_exp_sum_prev[TARGET_SEQ_LEN_BLOCK_SIZE];
 __local SOFTMAX_ACCUMULATOR_TYPE slm_max_val_prev[TARGET_SEQ_LEN_BLOCK_SIZE];
 {
#if IS_PAGED_ATTENTION
 const uint block_start_pos = blocked_indexes_start[target_seq_dim];
 const uint block_end_pos = blocked_indexes_end[target_seq_dim];
 uint query_offset = INPUT0_OFFSET +
 block_start_pos * (HEAD_SIZE * NUM_HEADS + INPUT0_PAD_BEFORE_FEATURE_NUM + INPUT0_PAD_AFTER_FEATURE_NUM) +
 num_heads_dim * HEAD_SIZE + head_size_idx;
 const uint query_pitch = (HEAD_SIZE * NUM_HEADS + INPUT0_PAD_BEFORE_FEATURE_NUM + INPUT0_PAD_AFTER_FEATURE_NUM);
 const uint cur_target_seq_len_size = block_end_pos - block_start_pos;
#else
#ifdef INPUT0_DIMS_ORDER
 uint query_offset = FUNC_CALL(get_input0_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, target_seq_idx, (head_size_idx));
 uint query_offset_next_seq = FUNC_CALL(get_input0_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, target_seq_idx + 1, (head_size_idx));
 const uint query_pitch = query_offset_next_seq - query_offset;
#else
 uint query_offset = INPUT0_GET_INDEX(b0_idx, b1_idx, target_seq_idx, (head_size_idx));
 const uint query_pitch = HEAD_SIZE;
#endif
 const uint cur_target_seq_len_size = min(TARGET_SEQ_LEN - target_seq_idx, (uint)TARGET_SEQ_LEN_BLOCK_SIZE);
#endif
 uint query_local_offset = head_size_idx * TARGET_SEQ_LEN_BLOCK_SIZE;
#if APPLY_SCALES_TO_QUERY
#if HAS_SCALE_INPUT
 const INPUT0_TYPE scale_val = *scale;
#else
 const INPUT0_TYPE scale_val = TO_INPUT0_TYPE(STATIC_SCALE_VALUE);
#endif
#else
 const INPUT0_TYPE scale_val = INPUT0_VAL_ONE;
#endif
 if (cur_target_seq_len_size != TARGET_SEQ_LEN_BLOCK_SIZE) {
 if (sgid * SUBGROUP_SIZE < HEAD_SIZE) {
 for (uint seq_idx = 0; seq_idx < cur_target_seq_len_size; seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 }
 } else {
 #if SG_SCALE_FACTOR == 2
 if ((sgid < (SUBGROUPS_PER_WG / SG_SCALE_FACTOR))) {
 unroll_for (uint seq_idx = 0; seq_idx < (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR); seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 } else {
 query_local_offset += (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR);
 query_offset += query_pitch * (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR);
 unroll_for (uint seq_idx = 0; seq_idx < (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR); seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 }
 #elif SG_SCALE_FACTOR == 4
 query_local_offset += (sgid / (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) * (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR);
 query_offset += query_pitch * (sgid / (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) * (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR);
 unroll_for (uint seq_idx = 0; seq_idx < (TARGET_SEQ_LEN_BLOCK_SIZE / SG_SCALE_FACTOR); seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 #else
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 INPUT0_TYPE val = BLOCK_READN(INPUT0_TYPE, 1, query_input, query_offset);
 slm_query[query_local_offset] = val * scale_val;
 query_offset += query_pitch;
 query_local_offset++;
 }
 #endif
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 }
 {
 #if TARGET_SEQ_LEN_BLOCK_SIZE <= SUBGROUP_SIZE
 if (sgid == 0 && sglid < TARGET_SEQ_LEN_BLOCK_SIZE) {
 slm_max_val_prev[sglid] = SOFTMAX_ACCUMULATOR_VAL_MIN;
 slm_exp_sum_prev[sglid] = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 }
 #else
 #error sdpa_opt.cl: unsupported TARGET_SEQ_LEN_BLOCK_SIZE
 #endif
 }
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) output_acc = OUTPUT_VAL_ZERO;
 __attribute__((opencl_unroll_hint(1)))
 for (uint start_partition_idx = 0; start_partition_idx < SOURCE_SEQ_LEN; start_partition_idx += SEQ_LEN_PARTITION_SIZE) {
 SOFTMAX_ACCUMULATOR_TYPE qk_max = SOFTMAX_ACCUMULATOR_VAL_MIN;
 const uint seq_len = start_partition_idx + sgid * SUBGROUP_SIZE;
 const uint partition_seq_len = min((uint)SOURCE_SEQ_LEN - start_partition_idx, (uint)SEQ_LEN_PARTITION_SIZE);
#if IS_PAGED_ATTENTION
#ifdef BROADCAST_GROUP_SIZE
 const uint heads_dim = num_heads_dim / BROADCAST_GROUP_SIZE;
#else
 const uint heads_dim = num_heads_dim;
#endif
 #define KEY_SEQ_OFFSET subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]]
 const uint key_pitch = (HEAD_SIZE * NUM_KV_HEADS + INPUT1_PAD_BEFORE_FEATURE_NUM + INPUT1_PAD_AFTER_FEATURE_NUM);
 uint key_offset = INPUT1_OFFSET +
 KEY_SEQ_OFFSET * key_pitch +
 heads_dim * HEAD_SIZE +
 seq_len * key_pitch;
#else
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_key)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, seq_len + sglid, 0)];
 const uint key_offset = FUNC_CALL(get_input1_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, seq_len + sglid, 0);
#else
 const uint b_idx = b0_idx;
 #ifdef INPUT1_DIMS_ORDER
 uint key_offset = FUNC_CALL(get_input1_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, seq_len, 0);
 uint key_offset_next_seq = FUNC_CALL(get_input1_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, seq_len + 1, 0);
 const uint key_pitch = key_offset_next_seq - key_offset;
 #else
 uint key_offset = INPUT1_GET_INDEX(b0_idx, b1_idx, seq_len, 0);
 const uint key_pitch = HEAD_SIZE;
 #endif
#endif
#endif
 int seq_len_calc_size = min((int)(SOURCE_SEQ_LEN) - (int)seq_len, (int)SUBGROUP_SIZE);
 MAKE_VECTOR_TYPE(INPUT0_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) qk_acc;
 qk_acc = FUNC_CALL(load_attn_mask)(OPTIONAL_SHAPE_INFO_TENSOR
 b0_idx,
 b1_idx,
#if IS_PAGED_ATTENTION
 blocked_indexes_start[target_seq_dim] - subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]] + sglid,
#else
 target_seq_idx + sglid,
#endif
 seq_len
 ATTN_MASK_BUFFER
 ATTN_SCALE_BUFFER
 PA_BUFFERS);
 if (seq_len_calc_size >= SUBGROUP_SIZE) {
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(KEY_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, seq_len + sglid, 0);
 KEY_COMPRESSION_SCALE_TYPE comp_scale = key_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 KEY_COMPRESSION_SCALE_TYPE comp_zp = key_scale[comp_offset + 1];
#endif
#endif
 __attribute__((opencl_unroll_hint(1)))
 for (uint head_idx_index = 0; head_idx_index < HEAD_SIZE; head_idx_index += SUBGROUP_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, 1, ptr, offset);
 #define QUERY_VEC MAKE_VECTOR_TYPE(INPUT0_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
 QUERY_VEC queries_vec;
 uint query_local_offset = (head_idx_index * TARGET_SEQ_LEN_BLOCK_SIZE) + sglid;
 unroll_for (uint q_row_idx = 0; q_row_idx < TARGET_SEQ_LEN_BLOCK_SIZE; q_row_idx++) {
 queries_vec[q_row_idx] = slm_query[query_local_offset];
 query_local_offset += TARGET_SEQ_LEN_BLOCK_SIZE;
 }
 unroll_for (uint key_row_idx = 0; key_row_idx < TARGET_SEQ_LEN_BLOCK_SIZE; key_row_idx++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT1_TYPE key_packed = KEY_BLOCK_READ(key_input, sub_group_broadcast(key_offset, key_row_idx) + head_idx_index);
#else
 const INPUT1_TYPE key_packed = KEY_BLOCK_READ(key_input, key_offset + key_row_idx * key_pitch + head_idx_index);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 KEY_COMPRESSION_SCALE_TYPE key_vals = (TO_KEY_COMPRESSION_SCALE_TYPE(key_packed) - sub_group_broadcast(comp_zp, key_row_idx)) * sub_group_broadcast(comp_scale, key_row_idx);
#elif IS_KV_COMPRESSED
 KEY_COMPRESSION_SCALE_TYPE key_vals = (TO_KEY_COMPRESSION_SCALE_TYPE(key_packed) * sub_group_broadcast(comp_scale, key_row_idx));
#else
 INPUT1_TYPE key_vals = key_packed;
#endif
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 qk_acc[key_row_idx] = mad(sub_group_broadcast(key_vals, i), queries_vec[i], qk_acc[key_row_idx]);
 }
 }
 }
 } else if (seq_len_calc_size > 0) {
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(KEY_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, seq_len + min(sglid, (uint)seq_len_calc_size - 1), 0);
 KEY_COMPRESSION_SCALE_TYPE comp_scale = key_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 KEY_COMPRESSION_SCALE_TYPE comp_zp = key_scale[comp_offset + 1];
#endif
#endif
 __attribute__((opencl_unroll_hint(1)))
 for (uint head_idx_index = 0; head_idx_index < HEAD_SIZE; head_idx_index += SUBGROUP_SIZE) {
 #define KEY_BLOCK_READ(ptr, offset) BLOCK_READN(INPUT1_TYPE, 1, ptr, offset)
 #define QUERY_VEC_TYPE MAKE_VECTOR_TYPE(INPUT0_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
#if IS_KV_COMPRESSED
 #define KEY_UNPACKED_TYPE KEY_COMPRESSION_SCALE_TYPE
 #define KEY_UNPACKED_VEC_TYPE MAKE_VECTOR_TYPE(KEY_COMPRESSION_SCALE_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
 #define TO_KEY_UNPACKED_TYPE(val) TO_KEY_COMPRESSION_SCALE_TYPE(val)
#else
 #define KEY_UNPACKED_TYPE INPUT1_TYPE
 #define KEY_UNPACKED_VEC_TYPE MAKE_VECTOR_TYPE(INPUT1_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE)
 #define TO_KEY_UNPACKED_TYPE(val) TO_INPUT1_TYPE(val)
#endif
 QUERY_VEC_TYPE queries_vec;
 uint query_local_offset = (head_idx_index * TARGET_SEQ_LEN_BLOCK_SIZE) + sglid;
 unroll_for (uint q_row_idx = 0; q_row_idx < TARGET_SEQ_LEN_BLOCK_SIZE; q_row_idx++) {
 queries_vec[q_row_idx] = slm_query[query_local_offset];
 query_local_offset += TARGET_SEQ_LEN_BLOCK_SIZE;
 }
#ifndef LOAD_KEY_LEFTOVERS_IN_CALC_LOOP
 KEY_UNPACKED_VEC_TYPE key_vec = 0;
 unroll_for (uint key_row_idx = 0; key_row_idx < seq_len_calc_size; key_row_idx++) {
#ifdef BEAM_TABLE_TYPE
 key_vec[key_row_idx] = TO_KEY_UNPACKED_TYPE(KEY_BLOCK_READ(key_input, sub_group_broadcast(key_offset, key_row_idx) + head_idx_index));
#else
 key_vec[key_row_idx] = TO_KEY_UNPACKED_TYPE(KEY_BLOCK_READ(key_input, key_offset + key_row_idx * key_pitch + head_idx_index));
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 key_vec[key_row_idx] = (key_vec[key_row_idx] - sub_group_broadcast(comp_zp, key_row_idx)) * sub_group_broadcast(comp_scale, key_row_idx);
#elif IS_KV_COMPRESSED
 key_vec[key_row_idx] *= sub_group_broadcast(comp_scale, key_row_idx);
#endif
 }
#endif
 unroll_for (uint key_row_idx = 0; key_row_idx < TARGET_SEQ_LEN_BLOCK_SIZE; key_row_idx++) {
#ifdef LOAD_KEY_LEFTOVERS_IN_CALC_LOOP
 KEY_UNPACKED_TYPE key_vals = 0;
 if (key_row_idx < seq_len_calc_size) {
#ifdef BEAM_TABLE_TYPE
 key_vals = TO_KEY_UNPACKED_TYPE(KEY_BLOCK_READ(key_input, sub_group_broadcast(key_offset, key_row_idx) + head_idx_index));
#else
 key_vals = TO_KEY_UNPACKED_TYPE(KEY_BLOCK_READ(key_input, key_offset + key_row_idx * key_pitch + head_idx_index));
#endif
 }
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 key_vals = (key_vals - sub_group_broadcast(comp_zp, key_row_idx)) * sub_group_broadcast(comp_scale, key_row_idx);
#elif IS_KV_COMPRESSED
 key_vals *= sub_group_broadcast(comp_scale, key_row_idx);
#endif
#else
 #define key_vals key_vec[key_row_idx]
#endif
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
 qk_acc[key_row_idx] = mad(sub_group_broadcast(key_vals, i), queries_vec[i], qk_acc[key_row_idx]);
 }
 }
 }
 }
 {
 unroll_for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
#if !APPLY_SCALES_TO_QUERY
#if HAS_SCALE_INPUT
 const OUTPUT_TYPE scale_val = *scale;
#else
 const OUTPUT_TYPE scale_val = TO_OUTPUT_TYPE(STATIC_SCALE_VALUE);
#endif
 qk_acc[i] *= scale_val;
#endif
#ifdef HAS_ALIBI
 const int alibi_val = (1 - SOURCE_SEQ_LEN) + seq_len + i;
 qk_acc[i] += alibi_slopes[num_heads_dim] * alibi_val;
#endif
 qk_acc[i] = INPUT0_MIN_FUNC(INPUT0_MAX_FUNC(qk_acc[i], INPUT0_VAL_MIN), INPUT0_VAL_MAX);
 qk_max = SOFTMAX_ACCUMULATOR_MAX_FUNC(qk_max, TO_SOFTMAX_ACCUMULATOR_TYPE(qk_acc[i]));
 }
 }
 {
 slm_qk_max_vals[sgid * SUBGROUP_SIZE + sglid] = qk_max;
 qk_max = SOFTMAX_ACCUMULATOR_VAL_MIN;
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 {
 SOFTMAX_ACCUMULATOR_TYPE qk_max_new = SOFTMAX_ACCUMULATOR_VAL_MIN;
 for (uint i = 0; i < SUBGROUPS_PER_WG; i++) {
 SOFTMAX_ACCUMULATOR_TYPE qk_max_val = slm_qk_max_vals[i * SUBGROUP_SIZE + sglid];
 qk_max_new = SOFTMAX_ACCUMULATOR_MAX_FUNC(qk_max_new, qk_max_val);
 }
 if (sgid == 0) {
 slm_max_val_cur[sglid] = qk_max_new;
 }
 SOFTMAX_ACCUMULATOR_TYPE exp_sum_new = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
 qk_acc[i] = native_exp(TO_SOFTMAX_ACCUMULATOR_TYPE(qk_acc[i]) - qk_max_new);
 exp_sum_new += qk_acc[i];
 }
 {
 slm_exp_sum_vals[sgid * SUBGROUP_SIZE + sglid] = exp_sum_new;
 }
 exp_sum_new = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 barrier(CLK_LOCAL_MEM_FENCE);
 for (uint i = 0; i < SUBGROUPS_PER_WG; i++) {
 SOFTMAX_ACCUMULATOR_TYPE exp_sum = slm_exp_sum_vals[i * SUBGROUP_SIZE + sglid];
 exp_sum_new += exp_sum;
 }
 for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
 qk_acc[i] = qk_acc[i] / exp_sum_new;
 }
 if (sgid == 0) {
 slm_exp_sum_cur[sglid] = exp_sum_new;
 }
 for (uint i = 0; i < TARGET_SEQ_LEN_BLOCK_SIZE; i++) {
 slm_qk_vals[sglid * SEQ_LEN_PARTITION_SIZE + sgid * TARGET_SEQ_LEN_BLOCK_SIZE + i] = qk_acc[i];
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 }
 {
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) acc_output_res = OUTPUT_VAL_ZERO;
#if IS_PAGED_ATTENTION
 const uint value_pitch = (HEAD_SIZE * NUM_KV_HEADS + INPUT2_PAD_BEFORE_FEATURE_NUM + INPUT2_PAD_AFTER_FEATURE_NUM);
#else
#ifdef INPUT2_DIMS_ORDER
 uint value_offset_base = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, 0, 0);
 uint value_offset_next_seq = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, 1, 0);
 const uint value_pitch = value_offset_next_seq - value_offset_base;
#else
 const uint value_pitch = HEAD_SIZE;
#endif
#endif
 if (partition_seq_len == SEQ_LEN_PARTITION_SIZE) {
 uint seq_len_start = (sgid / (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) * (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR);
 for (uint seq_len = seq_len_start; seq_len < seq_len_start + (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR); seq_len += SUBGROUP_SIZE) {
#if IS_PAGED_ATTENTION
#ifdef BROADCAST_GROUP_SIZE
 const uint heads_dim = num_heads_dim / BROADCAST_GROUP_SIZE;
#else
 const uint heads_dim = num_heads_dim;
#endif
 const uint value_seq_offset = subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]];
 uint value_offset = INPUT2_OFFSET +
 value_seq_offset * value_pitch +
 heads_dim * HEAD_SIZE +
 (start_partition_idx + (seq_len)) * value_pitch + head_size_idx;
#else
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len) + sglid, sgid * SUBGROUP_SIZE)];
 const uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + (seq_len) + sglid, sgid * SUBGROUP_SIZE);
#else
 const uint b_idx = b0_idx;
 #ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len), head_size_idx);
 #else
 uint value_offset = INPUT2_GET_INDEX(b0_idx, b1_idx, start_partition_idx + (seq_len), head_size_idx);
 #endif
#endif
#endif
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) qk_val;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = slm_qk_vals[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len + sglid];
 }
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + seq_len + sglid, 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, sub_group_broadcast(value_offset, i));
#else
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - sub_group_broadcast(comp_zp, i)) * sub_group_broadcast(comp_scale, i);
#elif IS_KV_COMPRESSED
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * sub_group_broadcast(comp_scale, i));
#else
 INPUT2_TYPE value_val = value_packed;
#endif
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc_output_res[seq_idx] = mad(sub_group_broadcast(qk_val[seq_idx], i), value_val, acc_output_res[seq_idx]);
 }
#ifndef BEAM_TABLE_TYPE
 value_offset += value_pitch;
#endif
 }
 }
 } else {
 const uint seq_len_start = (sgid / (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) * (SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR);
 uint seq_len_end = 0;
 if (seq_len_start < partition_seq_len)
 seq_len_end = seq_len_start + min(partition_seq_len - seq_len_start, (uint)(SEQ_LEN_PARTITION_SIZE / SG_SCALE_FACTOR));;
 for (uint seq_len = seq_len_start / SUBGROUP_SIZE; seq_len < seq_len_end / SUBGROUP_SIZE; seq_len++) {
#if IS_PAGED_ATTENTION
#ifdef BROADCAST_GROUP_SIZE
 const uint heads_dim = num_heads_dim / BROADCAST_GROUP_SIZE;
#else
 const uint heads_dim = num_heads_dim;
#endif
 const uint value_seq_offset = subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]];
 uint value_offset = INPUT2_OFFSET +
 value_seq_offset * value_pitch +
 heads_dim * HEAD_SIZE +
 (start_partition_idx + (seq_len * SUBGROUP_SIZE)) * value_pitch + head_size_idx;
#else
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, sgid * SUBGROUP_SIZE)];
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, sgid * SUBGROUP_SIZE);
#else
 const uint b_idx = b0_idx;
 #ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + (seq_len * SUBGROUP_SIZE), head_size_idx);
 #else
 uint value_offset = INPUT2_GET_INDEX(b0_idx, b1_idx, start_partition_idx + (seq_len * SUBGROUP_SIZE), head_size_idx);
 #endif
#endif
#endif
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + (seq_len * SUBGROUP_SIZE) + sglid, 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) qk_val;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = slm_qk_vals[seq_idx * SEQ_LEN_PARTITION_SIZE + seq_len * SUBGROUP_SIZE + sglid];
 }
 unroll_for (uint i = 0; i < SUBGROUP_SIZE; i++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, sub_group_broadcast(value_offset, i));
#else
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - sub_group_broadcast(comp_zp, i)) * sub_group_broadcast(comp_scale, i);
#elif IS_KV_COMPRESSED
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * sub_group_broadcast(comp_scale, i));
#else
 INPUT2_TYPE value_val = value_packed;
#endif
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc_output_res[seq_idx] = mad(sub_group_broadcast(qk_val[seq_idx], i), value_val, acc_output_res[seq_idx]);
 }
#ifndef BEAM_TABLE_TYPE
 value_offset += value_pitch;
#endif
 }
 }
 const uint seq_len_leftovers_start = ((seq_len_end / SUBGROUP_SIZE) * SUBGROUP_SIZE);
 if (seq_len_leftovers_start != seq_len_end) {
 uint qk_offset = min(seq_len_leftovers_start + sglid, seq_len_end - 1);
 MAKE_VECTOR_TYPE(OUTPUT_TYPE, TARGET_SEQ_LEN_BLOCK_SIZE) qk_val;
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 qk_val[seq_idx] = slm_qk_vals[qk_offset];
 qk_offset += SEQ_LEN_PARTITION_SIZE;
 }
#if IS_PAGED_ATTENTION
#ifdef BROADCAST_GROUP_SIZE
 const uint heads_dim = num_heads_dim / BROADCAST_GROUP_SIZE;
#else
 const uint heads_dim = num_heads_dim;
#endif
 const uint value_seq_offset = subsequence_begins[gws_seq_indexes_correspondence[target_seq_dim]];
 uint value_offset = INPUT2_OFFSET +
 value_seq_offset * value_pitch +
 heads_dim * HEAD_SIZE +
 (start_partition_idx + seq_len_leftovers_start) * value_pitch + head_size_idx;
#else
#ifdef BEAM_TABLE_TYPE
 const uint b_idx = beam_table[FUNC_CALL(get_bt_index_value)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + seq_len_leftovers_start + sglid, sgid * SUBGROUP_SIZE)];
 const uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b_idx, b1_idx, 0, 0, start_partition_idx + seq_len_leftovers_start + sglid, sgid * SUBGROUP_SIZE);
#else
 const uint b_idx = b0_idx;
 #ifdef INPUT2_DIMS_ORDER
 uint value_offset = FUNC_CALL(get_input2_index)(OPTIONAL_SHAPE_INFO_TENSOR b0_idx, b1_idx, 0, 0, start_partition_idx + seq_len_leftovers_start, head_size_idx);
 #else
 uint value_offset = INPUT2_GET_INDEX(b0_idx, b1_idx, start_partition_idx + seq_len_leftovers_start, head_size_idx);
 #endif
#endif
#endif
#if IS_KV_COMPRESSED
 const uint comp_offset = GET_COMPRESSION_INDEX(VALUE_COMPRESSION_SCALE, b_idx, b1_idx / BROADCAST_GROUP_SIZE, start_partition_idx + min(seq_len_leftovers_start + sglid, seq_len_end - 1), 0);
 VALUE_COMPRESSION_SCALE_TYPE comp_scale = val_scale[comp_offset];
#if USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE comp_zp = val_scale[comp_offset + 1];
#endif
#endif
 for (uint seq_len_idx = 0; seq_len_idx < partition_seq_len - seq_len_leftovers_start; seq_len_idx++) {
#ifdef BEAM_TABLE_TYPE
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, sub_group_broadcast(value_offset, seq_len_idx));
#else
 const INPUT2_TYPE value_packed = VALUE_BLOCK_READ(value_input, value_offset);
#endif
#if IS_KV_COMPRESSED && USE_ASYMMETRIC_QUANTIZATION
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed - sub_group_broadcast(comp_zp, seq_len_idx)) * sub_group_broadcast(comp_scale, seq_len_idx);
#elif IS_KV_COMPRESSED
 VALUE_COMPRESSION_SCALE_TYPE value_val = (value_packed * sub_group_broadcast(comp_scale, seq_len_idx));
#else
 INPUT2_TYPE value_val = value_packed;
#endif
 for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 acc_output_res[seq_idx] = mad(sub_group_broadcast(qk_val[seq_idx], seq_len_idx), value_val, acc_output_res[seq_idx]);
 }
#ifndef BEAM_TABLE_TYPE
 value_offset += value_pitch;
#endif
 }
 }
 }
 {
 SOFTMAX_ACCUMULATOR_TYPE exp_sum_prev = slm_exp_sum_prev[sglid];
 SOFTMAX_ACCUMULATOR_TYPE exp_sum_cur = slm_exp_sum_cur[sglid];
 SOFTMAX_ACCUMULATOR_TYPE max_val_prev = slm_max_val_prev[sglid];
 SOFTMAX_ACCUMULATOR_TYPE max_val_cur = slm_max_val_cur[sglid];
 barrier(CLK_LOCAL_MEM_FENCE);
#if IS_PAGED_ATTENTION
 const uint block_start_pos_new = blocked_indexes_start[target_seq_dim];
 const uint block_end_pos_new = blocked_indexes_end[target_seq_dim];
 const uint seq_idx_end = block_end_pos_new - block_start_pos_new;
#else
 const uint seq_idx_end = min(TARGET_SEQ_LEN - target_seq_idx, (uint)TARGET_SEQ_LEN_BLOCK_SIZE);
#endif
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE total_max = SOFTMAX_ACCUMULATOR_MAX_FUNC(sub_group_broadcast(max_val_prev, seq_idx), sub_group_broadcast(max_val_cur, seq_idx));
 SOFTMAX_ACCUMULATOR_TYPE updated_exp_sum_prev = sub_group_broadcast(exp_sum_prev, seq_idx) * native_exp(sub_group_broadcast(max_val_prev, seq_idx) - total_max);
 SOFTMAX_ACCUMULATOR_TYPE updated_exp_sum_cur = sub_group_broadcast(exp_sum_cur, seq_idx) * native_exp(sub_group_broadcast(max_val_cur, seq_idx) - total_max);
 SOFTMAX_ACCUMULATOR_TYPE updated_total_exp_sum = updated_exp_sum_prev + updated_exp_sum_cur;
 if (start_partition_idx > 0) {
 OUTPUT_TYPE updated_prev_res = TO_SOFTMAX_ACCUMULATOR_TYPE(output_acc[seq_idx]) * updated_exp_sum_prev / updated_total_exp_sum;;
 acc_output_res[seq_idx] *= updated_exp_sum_cur / updated_total_exp_sum;
 acc_output_res[seq_idx] += updated_prev_res;
 }
 output_acc[seq_idx] = acc_output_res[seq_idx];
 if (sgid == 0 && sglid == 0) {
 slm_exp_sum_prev[seq_idx] = updated_total_exp_sum;
 slm_max_val_prev[seq_idx] = total_max;
 }
 }
 }
 }
 }
 if (sgid >= (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) {
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 slm_qk_vals[seq_idx * SEQ_LEN_PARTITION_SIZE + (uint)get_local_id(2)] = output_acc[seq_idx];
 }
 }
 barrier(CLK_LOCAL_MEM_FENCE);
 if (sgid < (SUBGROUPS_PER_WG / SG_SCALE_FACTOR)) {
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 unroll_for (uint i = 1; i < SG_SCALE_FACTOR; i++) {
 output_acc[seq_idx] += slm_qk_vals[seq_idx * SEQ_LEN_PARTITION_SIZE + (i * HEAD_SIZE) + head_size_idx];
 }
 }
#if IS_PAGED_ATTENTION
 const uint block_start_pos_new = blocked_indexes_start[target_seq_dim];
 const uint block_end_pos_new = blocked_indexes_end[target_seq_dim];
 uint output_offset = block_start_pos_new * HEAD_SIZE * NUM_HEADS + num_heads_dim * HEAD_SIZE + sgid * SUBGROUP_SIZE;
 const uint output_pitch = HEAD_SIZE * NUM_HEADS;
#else
 uint output_offset = OUTPUT_GET_INDEX(b0_idx, b1_idx, target_seq_idx, sgid * SUBGROUP_SIZE);
 const uint output_pitch = HEAD_SIZE;
#endif
#if IS_PAGED_ATTENTION
 if (block_start_pos_new + TARGET_SEQ_LEN_BLOCK_SIZE != block_end_pos_new) {
 const uint seq_idx_end = block_end_pos_new - block_start_pos_new;
#else
 if (get_global_id(1) == get_global_size(1) - 1) {
 const uint seq_idx_end = min((uint)TARGET_SEQ_LEN - target_seq_idx, (uint)TARGET_SEQ_LEN_BLOCK_SIZE);
#endif
 for (uint seq_idx = 0; seq_idx < seq_idx_end; seq_idx++) {
 OUTPUT_BLOCK_WRITE(output, output_offset, output_acc[seq_idx]);
 output_offset += output_pitch;
 }
 } else {
 unroll_for (uint seq_idx = 0; seq_idx < TARGET_SEQ_LEN_BLOCK_SIZE; seq_idx++) {
 OUTPUT_BLOCK_WRITE(output, output_offset, output_acc[seq_idx]);
 output_offset += output_pitch;
 }
 }
 }
}
#endif
#endif
#ifdef SDPA_STAGE_1
#if SOFTMAX_ACCUMULATOR_TYPE_SIZE == 4
#define REG_VERSION_MAX_VALUES_PER_WI 24
#define REG_VERSION_MAX_VALUES_PER_WI_LOWER 8
#elif SOFTMAX_ACCUMULATOR_TYPE_SIZE == 2
#define REG_VERSION_MAX_VALUES_PER_WI 48
#define REG_VERSION_MAX_VALUES_PER_WI_LOWER 16
#else
#error Unexpected SOFTMAX_ACCUMULATOR data type size
#endif
REQD_SUB_GROUP_SIZE(SUBGROUP_SIZE)
KERNEL(sdpa_opt_finalization_stage)(
 OPTIONAL_SHAPE_INFO_ARG
 __global OUTPUT_TYPE* output,
 const __global SOFTMAX_ACCUMULATOR_TYPE* exp_sums,
 const __global SOFTMAX_ACCUMULATOR_TYPE* max_logits,
 const __global OUTPUT_TYPE* tmp_out,
 const uint num_of_partitions) {
 const uint batch_idx = get_global_id(0);
 const uint b0_idx = batch_idx / NUM_HEADS;
 const uint b1_idx = batch_idx % NUM_HEADS;
 const uint target_seq_idx = get_global_id(1);
 const uint sglid = get_sub_group_local_id();
 if (num_of_partitions <= SUBGROUP_SIZE * REG_VERSION_MAX_VALUES_PER_WI_LOWER) {
 SOFTMAX_ACCUMULATOR_TYPE exp_sum[REG_VERSION_MAX_VALUES_PER_WI_LOWER] = {SOFTMAX_ACCUMULATOR_VAL_ZERO};
 SOFTMAX_ACCUMULATOR_TYPE max_logit[REG_VERSION_MAX_VALUES_PER_WI_LOWER] = {SOFTMAX_ACCUMULATOR_VAL_MIN};
 SOFTMAX_ACCUMULATOR_TYPE local_exp_sum = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 SOFTMAX_ACCUMULATOR_TYPE local_max_logit = SOFTMAX_ACCUMULATOR_VAL_MIN;
 const uint iters_num = CEIL_DIV(num_of_partitions, SUBGROUP_SIZE);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 const uint max_logit_offset = exp_sums_offset;
 if (partition_idx < num_of_partitions) {
 exp_sum[i] = exp_sums[exp_sums_offset];
 max_logit[i] = max_logits[max_logit_offset];
 local_max_logit = SOFTMAX_ACCUMULATOR_MAX_FUNC(local_max_logit, max_logit[i]);
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_max = sub_group_reduce_max(local_max_logit);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 if (partition_idx < num_of_partitions) {
 exp_sum[i] = exp_sum[i] * native_exp(max_logit[i] - global_max);
 local_exp_sum += exp_sum[i];
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_sum = sub_group_reduce_add(local_exp_sum);
 for (uint head_size_idx = 0; head_size_idx < HEAD_SIZE / SUBGROUP_SIZE; head_size_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE acc = 0.0f;
 for (uint partition_idx = 0; partition_idx < num_of_partitions; partition_idx++) {
 const uint tmp_out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 target_seq_idx * (num_of_partitions * HEAD_SIZE) +
 partition_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 OUTPUT_TYPE out_val = tmp_out[tmp_out_offset];
 acc += TO_SOFTMAX_ACCUMULATOR_TYPE(out_val) *
 TO_SOFTMAX_ACCUMULATOR_TYPE(sub_group_broadcast(exp_sum[partition_idx / SUBGROUP_SIZE], partition_idx % SUBGROUP_SIZE)) /
 TO_SOFTMAX_ACCUMULATOR_TYPE(global_sum);
 }
 const uint out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * HEAD_SIZE) +
 target_seq_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 output[out_offset] = TO_OUTPUT_TYPE(acc);
 }
 } else if (num_of_partitions <= SUBGROUP_SIZE * REG_VERSION_MAX_VALUES_PER_WI) {
 SOFTMAX_ACCUMULATOR_TYPE exp_sum[REG_VERSION_MAX_VALUES_PER_WI] = {SOFTMAX_ACCUMULATOR_VAL_ZERO};
 SOFTMAX_ACCUMULATOR_TYPE max_logit[REG_VERSION_MAX_VALUES_PER_WI] = {SOFTMAX_ACCUMULATOR_VAL_MIN};
 SOFTMAX_ACCUMULATOR_TYPE local_exp_sum = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 SOFTMAX_ACCUMULATOR_TYPE local_max_logit = SOFTMAX_ACCUMULATOR_VAL_MIN;
 const uint iters_num = CEIL_DIV(num_of_partitions, SUBGROUP_SIZE);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 const uint max_logit_offset = exp_sums_offset;
 if (partition_idx < num_of_partitions) {
 exp_sum[i] = exp_sums[exp_sums_offset];
 max_logit[i] = max_logits[max_logit_offset];
 local_max_logit = SOFTMAX_ACCUMULATOR_MAX_FUNC(local_max_logit, max_logit[i]);
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_max = sub_group_reduce_max(local_max_logit);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 if (partition_idx < num_of_partitions) {
 exp_sum[i] = exp_sum[i] * native_exp(max_logit[i] - global_max);
 local_exp_sum += exp_sum[i];
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_sum = sub_group_reduce_add(local_exp_sum);
 for (uint head_size_idx = 0; head_size_idx < HEAD_SIZE / SUBGROUP_SIZE; head_size_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE acc = 0.0f;
 for (uint partition_idx = 0; partition_idx < num_of_partitions; partition_idx++) {
 const uint tmp_out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 target_seq_idx * (num_of_partitions * HEAD_SIZE) +
 partition_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 OUTPUT_TYPE out_val = tmp_out[tmp_out_offset];
 acc += TO_SOFTMAX_ACCUMULATOR_TYPE(out_val) *
 TO_SOFTMAX_ACCUMULATOR_TYPE(sub_group_broadcast(exp_sum[partition_idx / SUBGROUP_SIZE], partition_idx % SUBGROUP_SIZE)) /
 TO_SOFTMAX_ACCUMULATOR_TYPE(global_sum);
 }
 const uint out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * HEAD_SIZE) +
 target_seq_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 output[out_offset] = TO_OUTPUT_TYPE(acc);
 }
 } else {
 SOFTMAX_ACCUMULATOR_TYPE local_exp_sum = SOFTMAX_ACCUMULATOR_VAL_ZERO;
 SOFTMAX_ACCUMULATOR_TYPE local_max_logit = SOFTMAX_ACCUMULATOR_VAL_MIN;
 const uint iters_num = CEIL_DIV(num_of_partitions, SUBGROUP_SIZE);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 const uint max_logit_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 if (partition_idx < num_of_partitions) {
 local_max_logit = SOFTMAX_ACCUMULATOR_MAX_FUNC(local_max_logit, max_logits[max_logit_offset]);
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_max = sub_group_reduce_max(local_max_logit);
 for (uint i = 0; i < iters_num; i++) {
 const uint partition_idx = i * SUBGROUP_SIZE + sglid;
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 const uint max_logit_offset = exp_sums_offset;
 if (partition_idx < num_of_partitions) {
 local_exp_sum += exp_sums[exp_sums_offset] * native_exp(max_logits[max_logit_offset] - global_max);
 }
 }
 SOFTMAX_ACCUMULATOR_TYPE global_sum = sub_group_reduce_add(local_exp_sum);
 for (uint head_size_idx = 0; head_size_idx < HEAD_SIZE / SUBGROUP_SIZE; head_size_idx++) {
 SOFTMAX_ACCUMULATOR_TYPE acc = 0.0f;
 for (uint partition_idx = 0; partition_idx < num_of_partitions; partition_idx++) {
 const uint tmp_out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions * HEAD_SIZE) +
 target_seq_idx * (num_of_partitions * HEAD_SIZE) +
 partition_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 const uint exp_sums_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * num_of_partitions) +
 b1_idx * (TARGET_SEQ_LEN * num_of_partitions) +
 target_seq_idx * (num_of_partitions) +
 partition_idx;
 const uint max_logit_offset = exp_sums_offset;
 SOFTMAX_ACCUMULATOR_TYPE new_exp_sum = exp_sums[exp_sums_offset] * native_exp(max_logits[max_logit_offset] - global_max);
 OUTPUT_TYPE out_val = tmp_out[tmp_out_offset];
 acc += TO_SOFTMAX_ACCUMULATOR_TYPE(out_val) * new_exp_sum / TO_SOFTMAX_ACCUMULATOR_TYPE(global_sum);
 }
 const uint out_offset = b0_idx * (NUM_HEADS * TARGET_SEQ_LEN * HEAD_SIZE) +
 b1_idx * (TARGET_SEQ_LEN * HEAD_SIZE) +
 target_seq_idx * (HEAD_SIZE) +
 (head_size_idx * SUBGROUP_SIZE + sglid);
 output[out_offset] = TO_OUTPUT_TYPE(acc);
 }
 }
}
#endif
#ifdef OUTPUT_BLOCK_READ
#undef OUTPUT_BLOCK_READ
#endif
#ifdef OUTPUT_BLOCK_WRITE
#undef OUTPUT_BLOCK_WRITE
#endif
#ifdef VALUE_BLOCK_READ
#undef VALUE_BLOCK_READ
#endif
#ifdef SUBGROUPS_PER_WG
#undef SUBGROUPS_PER_WG
#endif
#ifdef GET_COMPRESSION_INDEX
#undef GET_COMPRESSION_INDEX
#endif
#ifdef GET_COMPRESSION_INDEX
#undef GET_COMPRESSION_INDEX
#endif
#ifdef QUERY_STEP_LOCAL
#undef QUERY_STEP_LOCAL
#endif
#ifdef QUERY_BLOCK_SIZE
#undef QUERY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_SIZE
#undef KEY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef KEY_BLOCK
#undef KEY_BLOCK
#endif
#ifdef KEY_BLOCK_UNCOMPRESSED
#undef KEY_BLOCK_UNCOMPRESSED
#endif
#ifdef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#undef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#endif
#ifdef QUERY_BLOCK
#undef QUERY_BLOCK
#endif
#ifdef KEY_BLOCK_SIZE
#undef KEY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef KEY_BLOCK
#undef KEY_BLOCK
#endif
#ifdef KEY_BLOCK_UNCOMPRESSED
#undef KEY_BLOCK_UNCOMPRESSED
#endif
#ifdef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#undef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#endif
#ifdef QUERY_BLOCK
#undef QUERY_BLOCK
#endif
#ifdef KEY_BLOCK_SIZE
#undef KEY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef KEY_BLOCK
#undef KEY_BLOCK
#endif
#ifdef KEY_BLOCK_UNCOMPRESSED
#undef KEY_BLOCK_UNCOMPRESSED
#endif
#ifdef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#undef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#endif
#ifdef QUERY_BLOCK
#undef QUERY_BLOCK
#endif
#ifdef KEY_BLOCK_SIZE
#undef KEY_BLOCK_SIZE
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef KEY_BLOCK
#undef KEY_BLOCK
#endif
#ifdef KEY_BLOCK_UNCOMPRESSED
#undef KEY_BLOCK_UNCOMPRESSED
#endif
#ifdef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#undef TO_KEY_BLOCK_UNCOMPRESSED_TYPE
#endif
#ifdef QUERY_BLOCK
#undef QUERY_BLOCK
#endif
#ifdef SOURCE_SEQ_LEN
#undef SOURCE_SEQ_LEN
#endif
#ifdef TARGET_SEQ_LEN
#undef TARGET_SEQ_LEN
#endif
#ifdef PA_BUFFERS
#undef PA_BUFFERS
#endif
#ifdef PA_BUFFERS_ARGS
#undef PA_BUFFERS_ARGS
#endif
#ifdef PA_BUFFERS
#undef PA_BUFFERS
#endif
#ifdef PA_BUFFERS_ARGS
#undef PA_BUFFERS_ARGS
#endif
#ifdef ATTN_MASK_BUFFER
#undef ATTN_MASK_BUFFER
#endif
#ifdef ATTN_MASK_BUFFER_ARG
#undef ATTN_MASK_BUFFER_ARG
#endif
#ifdef ATTN_MASK_BUFFER
#undef ATTN_MASK_BUFFER
#endif
#ifdef ATTN_MASK_BUFFER_ARG
#undef ATTN_MASK_BUFFER_ARG
#endif
#ifdef ATTN_SCALE_BUFFER
#undef ATTN_SCALE_BUFFER
#endif
#ifdef ATTN_SCALE_BUFFER_ARG
#undef ATTN_SCALE_BUFFER_ARG
#endif
#ifdef ATTN_SCALE_BUFFER
#undef ATTN_SCALE_BUFFER
#endif
#ifdef ATTN_SCALE_BUFFER_ARG
#undef ATTN_SCALE_BUFFER_ARG
#endif
#ifdef APPLY_SCALES_TO_QUERY
#undef APPLY_SCALES_TO_QUERY
#endif
#ifdef MASK_VECTOR_TYPE
#undef MASK_VECTOR_TYPE
#endif
#ifdef ALIBI_TYPE
#undef ALIBI_TYPE
#endif
#ifdef ALIBI_TYPE
#undef ALIBI_TYPE
#endif
#ifdef batch_idx
#undef batch_idx
#endif
#ifdef num_heads_dim
#undef num_heads_dim
#endif
#ifdef b0_idx
#undef b0_idx
#endif
#ifdef b1_idx
#undef b1_idx
#endif
#ifdef target_seq_dim
#undef target_seq_dim
#endif
#ifdef target_seq_idx
#undef target_seq_idx
#endif
#ifdef head_size_idx
#undef head_size_idx
#endif
#ifdef sglid
#undef sglid
#endif
#ifdef sgid
#undef sgid
#endif
#ifdef KEY_SEQ_OFFSET
#undef KEY_SEQ_OFFSET
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef QUERY_VEC
#undef QUERY_VEC
#endif
#ifdef KEY_BLOCK_READ
#undef KEY_BLOCK_READ
#endif
#ifdef QUERY_VEC_TYPE
#undef QUERY_VEC_TYPE
#endif
#ifdef KEY_UNPACKED_TYPE
#undef KEY_UNPACKED_TYPE
#endif
#ifdef KEY_UNPACKED_VEC_TYPE
#undef KEY_UNPACKED_VEC_TYPE
#endif
#ifdef TO_KEY_UNPACKED_TYPE
#undef TO_KEY_UNPACKED_TYPE
#endif
#ifdef KEY_UNPACKED_TYPE
#undef KEY_UNPACKED_TYPE
#endif
#ifdef KEY_UNPACKED_VEC_TYPE
#undef KEY_UNPACKED_VEC_TYPE
#endif
#ifdef TO_KEY_UNPACKED_TYPE
#undef TO_KEY_UNPACKED_TYPE
#endif
#ifdef key_vals
#undef key_vals
#endif
#ifdef REG_VERSION_MAX_VALUES_PER_WI
#undef REG_VERSION_MAX_VALUES_PER_WI
#endif
#ifdef REG_VERSION_MAX_VALUES_PER_WI_LOWER
#undef REG_VERSION_MAX_VALUES_PER_WI_LOWER
#endif
#ifdef REG_VERSION_MAX_VALUES_PER_WI
#undef REG_VERSION_MAX_VALUES_PER_WI
#endif
#ifdef REG_VERSION_MAX_VALUES_PER_WI_LOWER
#undef REG_VERSION_MAX_VALUES_PER_WI_LOWER
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
#ifdef BROADCAST_GROUP_SIZE
#undef BROADCAST_GROUP_SIZE
#endif
#ifdef DO_BROADCAST_KEY_VALUE
#undef DO_BROADCAST_KEY_VALUE
#endif
#ifdef IS_CAUSAL
#undef IS_CAUSAL
#endif
#ifdef HAS_ATTN_MASK_INPUT
#undef HAS_ATTN_MASK_INPUT
#endif
#ifdef HAS_SCALE_INPUT
#undef HAS_SCALE_INPUT
#endif
#ifdef IS_KV_COMPRESSED
#undef IS_KV_COMPRESSED
#endif
#ifdef INPUT0_DIMS_ORDER
#undef INPUT0_DIMS_ORDER
#endif
#ifdef INPUT1_DIMS_ORDER
#undef INPUT1_DIMS_ORDER
#endif
#ifdef INPUT2_DIMS_ORDER
#undef INPUT2_DIMS_ORDER
#endif
#ifdef TARGET_SEQ_LEN
#undef TARGET_SEQ_LEN
#endif
#ifdef NUM_HEADS
#undef NUM_HEADS
#endif
#ifdef NUM_KV_HEADS
#undef NUM_KV_HEADS
#endif
#ifdef SOURCE_SEQ_LEN
#undef SOURCE_SEQ_LEN
#endif
#ifdef SOFTMAX_ACCUMULATOR_TYPE
#undef SOFTMAX_ACCUMULATOR_TYPE
#endif
#ifdef SOFTMAX_ACCUMULATOR_VAL_MAX
#undef SOFTMAX_ACCUMULATOR_VAL_MAX
#endif
#ifdef SOFTMAX_ACCUMULATOR_VAL_MIN
#undef SOFTMAX_ACCUMULATOR_VAL_MIN
#endif
#ifdef SOFTMAX_ACCUMULATOR_VAL_ONE
#undef SOFTMAX_ACCUMULATOR_VAL_ONE
#endif
#ifdef SOFTMAX_ACCUMULATOR_VAL_ZERO
#undef SOFTMAX_ACCUMULATOR_VAL_ZERO
#endif
#ifdef TO_SOFTMAX_ACCUMULATOR_TYPE
#undef TO_SOFTMAX_ACCUMULATOR_TYPE
#endif
#ifdef TO_SOFTMAX_ACCUMULATOR_TYPE_SAT
#undef TO_SOFTMAX_ACCUMULATOR_TYPE_SAT
#endif
#ifdef AS_SOFTMAX_ACCUMULATOR_TYPE
#undef AS_SOFTMAX_ACCUMULATOR_TYPE
#endif
#ifdef SOFTMAX_ACCUMULATOR_MAX_FUNC
#undef SOFTMAX_ACCUMULATOR_MAX_FUNC
#endif
#ifdef SOFTMAX_ACCUMULATOR_MIN_FUNC
#undef SOFTMAX_ACCUMULATOR_MIN_FUNC
#endif
#ifdef SOFTMAX_ACCUMULATOR_ABS_FUNC
#undef SOFTMAX_ACCUMULATOR_ABS_FUNC
#endif
#ifdef SOFTMAX_ACCUMULATOR_TYPE_SIZE
#undef SOFTMAX_ACCUMULATOR_TYPE_SIZE
#endif
#ifdef SOFTMAX_ACCUMULATOR_IS_FP
#undef SOFTMAX_ACCUMULATOR_IS_FP
#endif
#ifdef SUBGROUP_SIZE
#undef SUBGROUP_SIZE
#endif
#ifdef HEAD_SIZE
#undef HEAD_SIZE
#endif
#ifdef SEQ_LEN_PARTITION_SIZE
#undef SEQ_LEN_PARTITION_SIZE
#endif
#ifdef TARGET_SEQ_LEN_BLOCK_SIZE
#undef TARGET_SEQ_LEN_BLOCK_SIZE
#endif
#ifdef SDPA_STAGE_0
#undef SDPA_STAGE_0
#endif
#ifdef SG_SCALE_FACTOR
#undef SG_SCALE_FACTOR
#endif
#ifdef STATIC_SCALE_VALUE_INV
#undef STATIC_SCALE_VALUE_INV
#endif
#ifdef STATIC_SCALE_VALUE
#undef STATIC_SCALE_VALUE
#endif

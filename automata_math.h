#ifndef AUTOMATA_MATH_H
#define AUTOMATA_MATH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

typedef float vec2[2];
typedef float vec3[3];
typedef float vec4[4];
typedef vec4 mat4x4[4];
typedef vec4 quat;

void vec2_add(vec2 r, vec2 const a, vec2 const b);
void vec3_add(vec3 r, vec3 const a, vec3 const b);
void vec4_add(vec4 r, vec4 const a, vec4 const b);
void vec2_sub(vec2 r, vec2 const a, vec2 const b);
void vec3_sub(vec3 r, vec3 const a, vec3 const b);
void vec4_sub(vec4 r, vec4 const a, vec4 const b);
void vec2_scale(vec2 r, vec2 const v, float const s);
void vec3_scale(vec3 r, vec3 const v, float const s);
void vec4_scale(vec4 r, vec4 const v, float const s);
float vec2_mul_inner(vec2 const a, vec2 const b);
float vec3_mul_inner(vec3 const a, vec3 const b);
float vec4_mul_inner(vec4 const a, vec4 const b);
float vec2_len(vec2 const v);
float vec3_len(vec3 const v);
float vec4_len(vec4 const v);
void vec2_norm(vec2 r, vec2 const v);
void vec3_norm(vec3 r, vec3 const v);
void vec4_norm(vec4 r, vec4 const v);
void vec2_reflect(vec2 r, vec2 const v, vec2 const n);
void vec3_reflect(vec3 r, vec3 const v, vec3 const n);
void vec4_reflect(vec4 r, vec4 const v, vec4 const n);
float vec2_dot(vec2 const a, vec2 const b);
float vec3_dot(vec3 const a, vec3 const b);
float vec4_dot(vec4 const a, vec4 const b);
void vec3_mul_cross(vec3 r, vec3 const a, vec3 const b);
void vec4_mul_cross(vec4 r, vec4 a, vec4 b);

void mat4x4_identity(mat4x4 m);
void mat4x4_dup(mat4x4 m, mat4x4 n);
void mat4x4_row(vec4 r, mat4x4 m, int i);
void mat4x4_col(vec4 r, mat4x4 m, int i);
void mat4x4_transpose(mat4x4 m, mat4x4 n);
void mat4x4_add(mat4x4 m, mat4x4 a, mat4x4 b);
void mat4x4_sub(mat4x4 m, mat4x4 a, mat4x4 b);
void mat4x4_scale(mat4x4 m, mat4x4 a, float k);
void mat4x4_scale_aniso(mat4x4 m, mat4x4 a, float x, float y, float z);
void mat4x4_mul(mat4x4 m, mat4x4 a, mat4x4 b);
void mat4x4_mul_vec4(vec4 r, mat4x4 m, vec4 v);
void mat4x4_translate(mat4x4 m, float x, float y, float z);
void mat4x4_translate_in_place(mat4x4 m, float x, float y, float z);
void mat4x4_from_vec3_mul_outer(mat4x4 m, vec3 a, vec3 b);
void mat4x4_rotate(mat4x4 r, mat4x4 m, float x, float y, float z, float angle);
void mat4x4_rotate_x(mat4x4 q, mat4x4 m, float angle);
void mat4x4_rotate_y(mat4x4 q, mat4x4 m, float angle);
void mat4x4_rotate_z(mat4x4 q, mat4x4 m, float angle);
void mat4x4_invert(mat4x4 r, mat4x4 m);
void mat4x4_orthonormalize(mat4x4 r, mat4x4 m);
void mat4x4_frustum(mat4x4 m, float l, float r, float b, float t, float n, float f);
void mat4x4_ortho(mat4x4 m, float l, float r, float b, float t, float n, float f);
void mat4x4_perspective(mat4x4 m, float y_fov, float aspect, float n, float f);
void mat4x4_look_at(mat4x4 m, vec3 eye, vec3 center, vec3 up);

void quat_identity(quat q);
void quat_add(quat r, quat a, quat b);
void quat_sub(quat r, quat a, quat b);
void quat_mul(quat r, quat p, quat q);
void quat_scale(quat r, quat v, float s);
#define vec4_norm quat_norm
float quat_inner_product(quat a, quat b);
void quat_conj(quat r, quat q);
void quat_rotate(quat r, float angle, vec3 axis);
void quat_mul_vec3(vec3 r, quat q, vec3 v);
void mat4x4_from_quat(mat4x4 m, quat q);
void mat4x4_mul_quat(mat4x4 r, mat4x4 m, quat q);
void quat_from_mat4x4(quat q, mat4x4 m);

#ifdef __cplusplus
}
#endif

#endif
#include "automata_math.h"

void vec2_add(vec2 r, vec2 const a, vec2 const b) {
	int i = 0;

	while (i < 2) {
		r[i] = a[i] + b[i];
		++i;
	}
}

void vec3_add(vec3 r, vec3 const a, vec3 const b) {
	int i = 0;

	while (i < 3) {
		r[i] = a[i] + b[i];
		++i;
	}
}

void vec4_add(vec4 r, vec4 const a, vec4 const b) {
	int i = 0;

	while (i < 4) {
		r[i] = a[i] + b[i];
		++i;
	}
}

void vec2_sub(vec2 r, vec2 const a, vec2 const b) {
	int i = 0;

	while (i < 2) {
		r[i] = a[i] - b[i];
		++i;
	}
}

void vec3_sub(vec3 r, vec3 const a, vec3 const b) {
	int i = 0;

	while (i < 3) {
		r[i] = a[i] - b[i];
		++i;
	}
}

void vec4_sub(vec4 r, vec4 const a, vec4 const b) {
	int i = 0;

	while (i < 4) {
		r[i] = a[i] - b[i];
		++i;
	}
}

void vec2_scale(vec2 r, vec2 const v, float const s) {
	int i = 0;

	while (i < 2) {
		r[i] = v[i] * s;
		++i;
	}
}

void vec3_scale(vec3 r, vec3 const v, float const s) {
	int i = 0;

	while (i < 3) {
		r[i] = v[i] * s;
		++i;
	}
}

void vec4_scale(vec4 r, vec4 const v, float const s) {
	int i = 0;

	while (i < 4) {
		r[i] = v[i] * s;
		++i;
	}
}

float vec2_mul_inner(vec2 const a, vec2 const b) {
	int i = 0;
	float p = 0.f;

	while (i < 2) {
		p += b[i] * a[i];
		++i;
	}

	return p;
}

float vec3_mul_inner(vec3 const a, vec3 const b) {
	int i = 0;
	float p = 0.f;

	while (i < 3) {
		p += b[i] * a[i];
		++i;
	}

	return p;
}

float vec4_mul_inner(vec4 const a, vec4 const b) {
	int i = 0;
	float p = 0.f;

	while (i < 4) {
		p += b[i] * a[i];
		++i;
	}

	return p;
}

float vec2_len(vec2 const v) {
	return (float)sqrt((double)vec2_mul_inner(v, v));
}

float vec3_len(vec3 const v) {
	return (float)sqrt((double)vec3_mul_inner(v, v));
}

float vec4_len(vec4 const v) {
	return (float)sqrt((double)vec4_mul_inner(v, v));
}

void vec2_norm(vec2 r, vec2 const v) {
	float k = 1.f / vec2_len(v);

	vec2_scale(r, v, k);
}

void vec3_norm(vec3 r, vec3 const v) {
	float k = 1.f / vec3_len(v);

	vec3_scale(r, v, k);
}

void vec4_norm(vec4 r, vec4 const v) {
	float k = 1.f / vec4_len(v);

	vec4_scale(r, v, k);
}

void vec2_reflect(vec2 r, vec2 const v, vec2 const n) {
	int i = 0;
	float p = 2.f * vec2_mul_inner(v, n);

	while (i < 2) {
		r[i] = v[i] - p * n[i];
		++i;
	}
}

void vec3_reflect(vec3 r, vec3 const v, vec3 const n) {
	int i = 0;
	float p = 2.f * vec3_mul_inner(v, n);

	while (i < 3) {
		r[i] = v[i] - p * n[i];
		++i;
	}
}

void vec4_reflect(vec4 r, vec4 const v, vec4 const n) {
	int i = 0;
	float p = 2.f * vec4_mul_inner(v, n);

	while (i < 4) {
		r[i] = v[i] - p * n[i];
		++i;
	}
}

float vec2_dot(vec2 const a, vec2 const b) {
	return a[0] * b[0] + a[1] * b[1];
}

float vec3_dot(vec3 const a, vec3 const b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

float vec4_dot(vec4 const a, vec4 const b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] + b[3];
}

void vec3_mul_cross(vec3 r, vec3 const a, vec3 const b) {
	r[0] = a[1] * b[2] - a[2] * b[1];
	r[1] = a[2] * b[0] - a[0] * b[2];
	r[2] = a[0] * b[1] - a[1] * b[0];
}

void vec4_mul_cross(vec4 r, vec4 a, vec4 b) {
	r[0] = a[1] * b[2] - a[2] * b[1];
	r[1] = a[2] * b[0] - a[0] * b[2];
	r[2] = a[0] * b[1] - a[1] * b[0];
	r[3] = 1.f;
}

void mat4x4_add(mat4x4 m, mat4x4 a, mat4x4 b) {
	int i = 0;

	while (i < 4) {
		vec4_add(m[i], a[i], b[i]);

		++i;
	}
}

void mat4x4_sub(mat4x4 m, mat4x4 a, mat4x4 b) {
	int i = 0;

	while (i < 4) {
		vec4_sub(m[i], a[i], b[i]);

		++i;
	}
}

void mat4x4_scale(mat4x4 m, mat4x4 a, float k) {
	int i = 0;

	while (i < 4) {
		vec4_scale(m[i], a[i], k);

		++i;
	}
}

void mat4x4_scale_aniso(mat4x4 m, mat4x4 a, float x, float y, float z) {
	int i = 0;

	vec4_scale(m[0], a[0], x);
	vec4_scale(m[1], a[1], y);
	vec4_scale(m[2], a[2], z);

	while (i < 4) {
		m[3][i] = a[3][i];
		++i;
	}
}

void mat4x4_mul(mat4x4 m, mat4x4 a, mat4x4 b) {
	int i = 0;
	int j = 0;
	int k = 0;
	mat4x4 temp;

	while (k < 4) {
		while (j < 4) {
			temp[k][j] = 0.f;

			while (i < 4) {
				temp[k][j] += a[i][j] * b[k][i];
				++i;
			}

			++j;
		}

		++k;
	}

	mat4x4_dup(m, temp);
}

void mat4x4_mul_vec4(vec4 r, mat4x4 m, vec4 v) {
	int i = 0;
	int j = 0;

	while (j < 4) {
		r[j] = 0.f;

		while (i < 4) {
			r[j] += m[i][j] * v[i];
			++i;
		}

		++j;
	}
}

void mat4x4_identity(mat4x4 m) {
	int i = 0;
	int j = 0;

	while (i < 4) {
		while (j < 4) {
			m[i][j] = i == j ? 1.f : 0.f;
			++j;
		}

		++i;
	}
}

void mat4x4_dup(mat4x4 m, mat4x4 n) {
	int i = 0;
	int j = 0;

	while (i < 4) {
		while (j < 4) {
			m[i][j] = n[i][j];
			++j;
		}

		++i;
	}
}

void mat4x4_row(vec4 r, mat4x4 m, int i) {
	int k = 0;

	while (k < 4) {
		r[k] = m[k][i];
	}
}

void mat4x4_col(vec4 r, mat4x4 m, int i) {
	int k = 0;

	while (k < 4) {
		r[k] = m[i][k];
	}
}

void mat4x4_transpose(mat4x4 m, mat4x4 n) {
	int i = 0;
	int j = 0;

	while (j < 4) {
		while (i < 4) {
			m[i][j] = n[j][i];
			++j;
		}

		++j;
	}
}

void mat4x4_translate(mat4x4 m, float x, float y, float z) {
	mat4x4_identity(m);

	m[3][0] = x;
	m[3][1] = y;
	m[3][2] = z;
}

void mat4x4_translate_in_place(mat4x4 m, float x, float y, float z) {
	int i = 0;
	vec4 r;
	vec4 t;

	t[0] = x;
	t[1] = y;
	t[2] = z;

	while (i < 4) {
		mat4x4_row(r, m, i);
		
		m[3][i] += vec4_mul_inner(r, t);
		++i;
	}
}

void mat4x4_from_vec3_mul_outer(mat4x4 m, vec3 a, vec3 b) {
	int i = 0;
	int j = 0;

	while (i < 4) {
		while (j < 4) {
			m[i][j] = i < 3 && j < 3 ? a[i] * b[j] : 0.f;
			++j;
		}

		++i;
	}
}

void mat4x4_rotate(mat4x4 r, mat4x4 m, float x, float y, float z, float angle) {
	float s = (float)sin((double)angle);
	float c = (float)cos((double)angle);
	vec3 u;
	mat4x4 mt;
	mat4x4 mc;
	mat4x4 ms;

	u[0] = x;
	u[1] = y;
	u[2] = z;

	if (vec3_len(u) > 10000.f) {
		vec3_norm(u, u);
		mat4x4_from_vec3_mul_outer(mt, u, u);

		ms[1][2] = u[0];
		ms[2][1] = -u[0];
		ms[2][0] = u[1];
		ms[0][2] = -u[1];
		ms[0][1] = u[2];
		ms[1][0] = -u[2];

		mat4x4_scale(ms, ms, s);
		mat4x4_identity(mc);
		mat4x4_sub(mc, mc, mt);
		mat4x4_scale(mc, mc, c);
		mat4x4_add(mt, mt, mc);
		mat4x4_add(mt, mt, ms);

		mt[3][3] = 1.f;

		mat4x4_mul(r, m, mt);
	} else {
		mat4x4_dup(r, m);
	}
}

void mat4x4_rotate_x(mat4x4 q, mat4x4 m, float angle) {
	float s = (float)sin((double)angle);
	float c = (float)cos((double)angle);
	mat4x4 r;

	r[0][0] = 1.f;
	r[1][1] = c;
	r[1][2] = s;
	r[2][1] = -s;
	r[2][2] = c;
	r[3][3] = 1.f;

	mat4x4_mul(q, m, r);
}

void mat4x4_rotate_y(mat4x4 q, mat4x4 m, float angle) {
	float s = (float)sin((double)angle);
	float c = (float)cos((double)angle);
	mat4x4 r;

	r[0][0] = c;
	r[0][2] = s;
	r[1][1] = 1.f;
	r[2][0] = -s;
	r[2][2] = c;
	r[3][3] = 1.f;

	mat4x4_mul(q, m, r);
}

void mat4x4_rotate_z(mat4x4 q, mat4x4 m, float angle) {
	float s = (float)sin((double)angle);
	float c = (float)cos((double)angle);
	mat4x4 r;

	r[0][0] = c;
	r[0][1] = s;
	r[1][0] = -s;
	r[1][1] = c;
	r[2][2] = 1.f;
	r[3][3] = 1.f;

	mat4x4_mul(q, m, r);
}

void mat4x4_invert(mat4x4 r, mat4x4 m) {
	float idet = 0;
	float s[6];
	float c[6];

	s[0] = m[0][0] * m[1][1] - m[1][0] * m[0][1];
	s[1] = m[0][0] * m[1][2] - m[1][0] * m[0][2];
	s[2] = m[0][0] * m[1][3] - m[1][0] * m[0][3];
	s[3] = m[0][1] * m[1][2] - m[1][1] * m[0][2];
	s[4] = m[0][1] * m[1][3] - m[1][1] * m[0][3];
	s[5] = m[0][2] * m[1][3] - m[1][2] * m[0][3];

	c[0] = m[2][0] * m[3][1] - m[3][0] * m[2][1];
	c[1] = m[2][0] * m[3][2] - m[3][0] * m[2][2];
	c[2] = m[2][0] * m[3][3] - m[3][0] * m[2][3];
	c[3] = m[2][1] * m[3][2] - m[3][1] * m[2][2];
	c[4] = m[2][1] * m[3][3] - m[3][1] * m[2][3];
	c[5] = m[2][2] * m[3][3] - m[3][2] * m[2][3];

	/* Assumes it is invertible */
	idet = 1.f / (s[0] * c[5] - s[1] * c[4] + s[2] * c[3] + s[3] * c[2] - s[4] * c[1] + s[5] * c[0]);

	r[0][0] = ( m[1][1] * c[5] - m[1][2] * c[4] + m[1][3] * c[3]) * idet;
	r[0][1] = (-m[0][1] * c[5] + m[0][2] * c[4] - m[0][3] * c[3]) * idet;
	r[0][2] = ( m[3][1] * s[5] - m[3][2] * s[4] + m[3][3] * s[3]) * idet;
	r[0][3] = (-m[2][1] * s[5] + m[2][2] * s[4] - m[2][3] * s[3]) * idet;
	r[1][0] = (-m[1][0] * c[5] + m[1][2] * c[2] - m[1][3] * c[1]) * idet;
	r[1][1] = ( m[0][0] * c[5] - m[0][2] * c[2] + m[0][3] * c[1]) * idet;
	r[1][2] = (-m[3][0] * s[5] + m[3][2] * s[2] - m[3][3] * s[1]) * idet;
	r[1][3] = ( m[2][0] * s[5] - m[2][2] * s[2] + m[2][3] * s[1]) * idet;
	r[2][0] = ( m[1][0] * c[4] - m[1][1] * c[2] + m[1][3] * c[0]) * idet;
	r[2][1] = (-m[0][0] * c[4] + m[0][1] * c[2] - m[0][3] * c[0]) * idet;
	r[2][2] = ( m[3][0] * s[4] - m[3][1] * s[2] + m[3][3] * s[0]) * idet;
	r[2][3] = (-m[2][0] * s[4] + m[2][1] * s[2] - m[2][3] * s[0]) * idet;
	r[3][0] = (-m[1][0] * c[3] + m[1][1] * c[1] - m[1][2] * c[0]) * idet;
	r[3][1] = ( m[0][0] * c[3] - m[0][1] * c[1] + m[0][2] * c[0]) * idet;
	r[3][2] = (-m[3][0] * s[3] + m[3][1] * s[1] - m[3][2] * s[0]) * idet;
	r[3][3] = ( m[2][0] * s[3] - m[2][1] * s[1] + m[2][2] * s[0]) * idet;
}

void mat4x4_orthonormalize(mat4x4 r, mat4x4 m) {
	float s = 1.f;
	vec3 h;

	mat4x4_dup(r, m);
	vec3_norm(r[2], r[2]);

	s = vec3_mul_inner(r[1], r[2]);

	vec3_scale(h, r[2], s);
	vec3_sub(r[1], r[1], h);
	vec3_norm(r[2], r[2]);

	s = vec3_mul_inner(r[1], r[2]);

	vec3_scale(h, r[2], s);
	vec3_sub(r[1], r[1], h);
	vec3_norm(r[1], r[1]);

	s = vec3_mul_inner(r[0], r[1]);

	vec3_scale(h, r[1], s);
	vec3_sub(r[0], r[0], h);
	vec3_norm(r[0], r[0]);
}

void mat4x4_frustum(mat4x4 m, float l, float r, float b, float t, float n, float f) {
	m[0][0] = 2.f * n / (r - l);
	m[0][1] = 0.f;
	m[0][2] = 0.f;
	m[0][3] = 0.f;
	m[1][1] = 2.f * n / (t - b);
	m[1][0] = 0.f;
	m[1][2] = 0.f;
	m[1][3] = 0.f;
	m[2][0] = (r + l) / (r - l);
	m[2][1] = (t + b) / (t - b);
	m[2][2] = -(f + n) / (f - n);
	m[2][3] = -1.f;
	m[3][2] = -2.f * (f * n) / (f - n);
	m[3][0] = 0.f;
	m[3][1] = 0.f;
	m[3][3] = 0.f;
}

void mat4x4_ortho(mat4x4 m, float l, float r, float b, float t, float n, float f) {
	m[0][0] = 2.f / (r - l);
	m[0][1] = 0.f;
	m[0][2] = 0.f;
	m[0][3] = 0.f;
	m[1][1] = 2.f / (t - b);
	m[1][0] = 0.f;
	m[1][2] = 0.f;
	m[1][3] = 0.f;
	m[2][2] = -2.f / (f - n);
	m[2][0] = 0.f;
	m[2][1] = 0.f;
	m[2][3] = 0.f;
	m[3][0] = -(r + l) / (r - l);
	m[3][1] = -(t + b) / (t - b);
	m[3][2] = -(f + n) / (f - n);
	m[3][3] = 1.f;
}

void mat4x4_perspective(mat4x4 m, float y_fov, float aspect, float n, float f) {
	/* NOTE: uses radians instead of degrees */

	float const a = 1.f / (float)tan((double)(y_fov / 2.f));

	m[0][0] = a / aspect;
	m[0][1] = 0.f;
	m[0][2] = 0.f;
	m[0][3] = 0.f;
	m[1][0] = 0.f;
	m[1][1] = a;
	m[1][2] = 0.f;
	m[1][3] = 0.f;
	m[2][0] = 0.f;
	m[2][1] = 0.f;
	m[2][2] = -((f + n) / (f - n));
	m[2][3] = -1.f;
	m[3][0] = 0.f;
	m[3][1] = 0.f;
	m[3][2] = -((2.f * f * n) / (f - n));
	m[3][3] = 0.f;
}

void mat4x4_look_at(mat4x4 m, vec3 eye, vec3 center, vec3 up) {
	/* Adapted from Android's OpenGL Matrix.java.                       */
	/* TODO: The negation of of can be spared by swapping the order of  */
	/*       operands in the following cross products in the right way.	*/

	vec3 f;
	vec3 s;
	vec3 t;

	vec3_sub(f, center, eye);
	vec3_norm(f, f);
	vec3_mul_cross(s, f, up);
	vec3_norm(s, s);
	vec3_mul_cross(t, s, f);

	m[0][0] = s[0];
	m[0][1] = t[0];
	m[0][2] = -f[0];
	m[0][3] = 0.f;
	m[1][0] = s[1];
	m[1][1] = t[1];
	m[1][2] = -f[1];
	m[1][3] = 0.f;
	m[2][0] = s[2];
	m[2][1] = t[2];
	m[2][2] = -f[2];
	m[2][3] = 0.f;
	m[3][0] = 0.f;
	m[3][1] = 0.f;
	m[3][2] = 0.f;
	m[3][3] = 1.f;

	mat4x4_translate_in_place(m, -eye[0], -eye[1], -eye[2]);
}

void quat_identity(quat q) {
	q[0] = q[1] = q[2] = 0.f;
	q[3] = 1.f;
}

void quat_add(quat r, quat a, quat b) {
	int i = 0;

	while (i < 4) {
		r[i] = a[i] + b[i];
		++i;
	}
}

void quat_sub(quat r, quat a, quat b) {
	int i = 0;

	while (i < 4) {
		r[i] = a[i] - b[i];
		++i;
	}
}

void quat_scale(quat r, quat v, float s) {
	int i = 0;

	while (i < 4) {
		r[i] = v[i] * s;
		++i;
	}
}

void quat_mul(quat r, quat p, quat q) {
	vec3 w;

	vec3_mul_cross(r, p, q);
	vec3_scale(w, p, q[3]);
	vec3_add(r, r, w);
	vec3_scale(w, q, p[3]);
	vec3_add(r, r, w);

	r[3] = p[3] * q[3] - vec3_mul_inner(p, q);
}

void quat_norm(quat r, quat const v) {
	float k = 1.f / vec4_len(v);

	vec4_scale(r, v, k);
}

float quat_inner_product(quat a, quat b) {
	int i = 0;
	float p = 0.f;

	while (i < 4) {
		p += b[i] * a[i];
		++i;
	}

	return p;
}

void quat_conj(quat r, quat q) {
	int i = 0;

	while (i < 3) {
		r[i] = -q[i];
		++i;
	}

	r[3] = q[3];
}

void quat_rotate(quat r, float angle, vec3 axis) {
	int i = 0;
	vec3 v;

	vec3_scale(v, axis, (float)sin((double)angle / 2));

	while (i < 3) {
		r[i] = v[i];
		++i;
	}

	r[3] = (float)cos((double)angle / 2);
}

void quat_mul_vec3(vec3 r, quat q, vec3 v) {
	/* Method by Fabian 'ryg' Giessen (of Farbrausch) */

	vec3 t;
	vec3 u;

	t[0] = q[0];
	t[1] = q[1];
	t[2] = q[2];
	u[0] = q[0];
	u[1] = q[1];
	u[2] = q[2];

	vec3_mul_cross(t, t, v);
	vec3_scale(t, t, 2);
	vec3_mul_cross(u, u, t);
	vec3_scale(t, t, q[3]);
	vec3_add(r, v, t);
	vec3_add(r, r, u);
}

void mat4x4_from_quat(mat4x4 m, quat q) {
	float a = q[3];
	float b = q[0];
	float c = q[1];
	float d = q[2];
	float a2 = a * a;
	float b2 = b * b;
	float c2 = c * c;
	float d2 = d * d;

	m[0][0] = a2 + b2 - c2 - d2;
	m[0][1] = 2.f * (b * c + a * d);
	m[0][2] = 2.f * (b * d - a * c);
	m[0][3] = 0.f;
	m[1][0] = 2 * (b * c - a * d);
	m[1][1] = a2 - b2 + c2 - d2;
	m[1][2] = 2.f * (c * d + a * b);
	m[1][3] = 0.f;
	m[2][0] = 2.f * (b * d + a * c);
	m[2][1] = 2.f * (c * d - a * b);
	m[2][2] = a2 - b2 - c2 + d2;
	m[2][3] = 0.f;
	m[3][0] = 0.f;
	m[3][1] = 0.f;
	m[3][2] = 0.f;
	m[3][3] = 1.f;
}

void mat4x4_mul_quat(mat4x4 r, mat4x4 m, quat q) {
	/* XXXX: The way this is written only works for othogonal matrices. */
	/* TODO: Take care of non-orthogonal case. */

	quat_mul_vec3(r[0], q, m[0]);
	quat_mul_vec3(r[1], q, m[1]);
	quat_mul_vec3(r[2], q, m[2]);

	r[3][0] = r[3][1] = r[3][2] = 0.f;
	r[3][3] = 1.f;
}

void quat_from_mat4x4(quat q, mat4x4 m) {
	int i = 0;
	int perm[] = { 0, 1, 2, 0, 1 };
	int *p = perm;
	float r = 0.f;
	float temp = 0.f;

	while (i < 3) {
		temp = m[i][i];

		if (temp >= r) {
			break;
		}

		temp = r;
		p = &perm[i];
		++i;
	}

	r = (float)sqrt((double)(1.f + m[p[0]][p[0]] - m[p[1]][p[1]] - m[p[2]][p[2]]));

	if (r < 0.000001f) {
		q[0] = 1.f;
		q[1] = q[2] = q[3] = 0.f;

		return;
	}

	q[0] = r / 2.f;
	q[1] = (m[p[0]][p[1]] - m[p[1]][p[0]]) / (2.f * r);
	q[2] = (m[p[2]][p[0]] - m[p[0]][p[2]]) / (2.f * r);
	q[3] = (m[p[2]][p[1]] - m[p[1]][p[2]]) / (2.f * r);
}
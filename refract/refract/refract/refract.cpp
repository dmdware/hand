// refract.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

struct v3f
{
	float x, y, z;
};

typedef struct v3f v3f;

struct v4f
{
	float x, y, z, w;
};

typedef struct v4f v4f;

void v3fmul(v3f *v, v3f f, float s)
{
	v->x = f.x * s;
	v->y = f.y * s;
	v->z = f.z * s;
}

void v3fsub(v3f *v, v3f f, v3f s)
{
	v->x = f.x - s.x;
	v->y = f.y - s.y;
	v->z = f.z - s.z;
}

void v3fadd(v3f *v, v3f f, v3f s)
{
	v->x = f.x + s.x;
	v->y = f.y + s.y;
	v->z = f.z + s.z;
}

float mag3f(v3f v)
{
	return (float)sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

float mag4f(v4f v)
{
	return (float)sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z) + (v.w * v.w));
}

float dot3f(v3f v1, v3f v2)
{
	return ((v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z));
}

v3f cross3f(v3f v1, v3f v2)
{
	v3f n;

	n.x = ((v1.y * v2.z) - (v1.z * v2.y));
	n.y = ((v1.z * v2.x) - (v1.x * v2.z));
	n.z = ((v1.x * v2.y) - (v1.y * v2.x));

	return n;
}

float vol3f2(v3f a, v3f b, v3f c)
{
	float f;
	b = cross3f(b, c);
	f = dot3f(a, b);
	f = fabs(f);
	f *= 1.0 / (2.0 * 3.0);
	return f;
}

float vol3f(v3f a, v3f b, v3f c, v3f d)
{
	v3fsub(&a, a, d);
	v3fsub(&b, b, d);
	v3fsub(&c, c, d);
	return vol3f2(a, b, c);
}

float sa3f(v3f a, v3f b, v3f c)
{
	float s, a1, a2, a3;
	v3f ab, ac, bc;

	v3fsub(&ab, a, b);
	v3fsub(&ac, a, c);
	v3fsub(&bc, b, c);

	a1 = mag3f(ab);
	a2 = mag3f(ac);
	a3 = mag3f(bc);

	s = (a1 + a2 + a3) / 2.0f;

	return sqrtf(s * (s - a1) * (s - a2) * (s - a3));
}

float vol4f2(float s1, float s2, float s3, float s4, float a1, float a2, float a3, float a4, float a5, float a6)
{
	return sqrtf(s1*(s1 - a1)*(s1 - a2)*(s1 - a4)) +
		sqrtf(s2*(s2 - a2)*(s2 - a3)*(s2 - a5)) +
		sqrtf(s3*(s3 - a3)*(s3 - a6)*(s3 - a1)) +
		sqrtf(s4*(s4 - a4)*(s4 - a5)*(s4 - a6));
}

void v4fsub(v4f *v, v4f v1, v4f v2)
{
	v->x = v1.x - v2.x;
	v->y = v1.y - v2.y;
	v->z = v1.z - v2.z;
	v->w = v1.w - v2.w;
}

float vol4f(v4f a, v4f b, v4f c, v4f d)
{
	float s1, s2, s3, s4,
		a1, a2, a3, a4, a5, a6;
	v4f ab, ac, ad, bc, bd, cd;

	v4fsub(&ab, a, b);
	v4fsub(&ac, a, c);
	v4fsub(&ad, a, d);
	v4fsub(&bc, b, c);
	v4fsub(&bd, b, d);
	v4fsub(&cd, c, d);

	a1 = mag4f(ab);
	a2 = mag4f(ac);
	a3 = mag4f(ad);
	a4 = mag4f(bc);
	a5 = mag4f(bd);
	a6 = mag4f(cd);

	s1 = (a1 + a2 + a4) / 2.0f;
	s1 = (a2 + a3 + a5) / 2.0f;
	s1 = (a3 + a6 + a1) / 2.0f;
	s1 = (a4 + a5 + a6) / 2.0f;

	return vol4f2(s1, s2, s3, s4, a1, a2, a3, a4, a5, a6);
}

float sa4f(v4f a, v4f b, v4f c)
{
	float s, a1, a2, a3;
	v4f ab, ac, bc;

	v4fsub(&ab, a, b);
	v4fsub(&ac, a, c);
	v4fsub(&bc, b, c);

	a1 = mag4f(ab);
	a2 = mag4f(ac);
	a3 = mag4f(bc);

	s = (a1 + a2 + a3) / 2.0f;

	return sqrtf(s * (s - a1) * (s - a2) * (s - a3));
}

v3f toxy2(v3f vi, float wx, float wy, v3f p[8])
{
	float v[12];
	int i;
	v[0] = vol3f(vi, p[0], p[1], p[2]) / sa3f(p[1], p[2], p[0]);
	v[1] = vol3f(vi, p[0], p[2], p[3]) / sa3f(p[0], p[3], p[2]);
	v[2] = vol3f(vi, p[1], p[2], p[6]) / sa3f(p[1], p[2], p[6]);
	v[3] = vol3f(vi, p[1], p[6], p[5]) / sa3f(p[1], p[6], p[5]);
	v[4] = vol3f(vi, p[2], p[7], p[6]) / sa3f(p[2], p[6], p[7]);
	v[5] = vol3f(vi, p[3], p[7], p[2]) / sa3f(p[7], p[3], p[2]);
	v[6] = vol3f(vi, p[3], p[0], p[7]) / sa3f(p[7], p[3], p[0]);
	v[7] = vol3f(vi, p[0], p[7], p[4]) / sa3f(p[7], p[4], p[0]);
	v[8] = vol3f(vi, p[1], p[5], p[4]) / sa3f(p[1], p[5], p[4]);
	v[9] = vol3f(vi, p[0], p[1], p[4]) / sa3f(p[0], p[1], p[4]);
	v[10] = vol3f(vi, p[5], p[6], p[7]) / sa3f(p[5], p[6], p[7]);
	v[11] = vol3f(vi, p[4], p[5], p[7]) / sa3f(p[4], p[5], p[7]);
	for (i = 0; i < 12; i++)
		printf("v[%d]=%f\r\n", i, v[i]);
	vi.x = wx * (v[4] + v[5]) / (v[8] + v[9] + v[4] + v[5]);
	vi.y = wy * (v[2] + v[3]) / (v[6] + v[7] + v[2] + v[3]);
	vi.z = ((v[10] + v[11]) / (v[0] + v[1] + v[10] + v[11]));
	return vi;
}

v3f norm3f(v3f n)
{
	float m;

	m = mag3f(n);

	n.x /= m;
	n.y /= m;
	n.z /= m;

	return n;
}



void v4freset(v4f *v)
{
	v->x = v->y = v->z = v->w = 0;
}

void v4fset(v4f *v, float x, float y, float z, float w)
{
	v->x = x;
	v->y = y;
	v->z = z;
	v->w = w;
}

void v4fset2(v4f *v, const float *values)
{
	v->x = values[0];
	v->y = values[1];
	v->z = values[2];
	v->w = values[3];
}

void v4fadd(v4f *v, v4f v1, v4f v2)
{
	v->x = v1.x + v2.x;
	v->y = v1.y + v2.y;
	v->z = v1.z + v2.z;
	v->w = v1.w + v2.w;
}

void v4fmul(v4f *v, v4f v1, const float num)
{
	v->x = v1.x * num;
	v->y = v1.y * num;
	v->z = v1.z * num;
	v->w = v1.w * num;
}

void v4fdiv(v4f *v, const v4f v1, const float num)
{
	v->x = v1.x / num;
	v->y = v1.y / num;
	v->z = v1.z / num;
	v->w = v1.w / num;
}

void trilinear(v3f f, v4f p[2][2][2], v4f* out)
{
	/*
	Vxyz =	V000 (1 - x) (1 - y) (1 - z) +
	V100 x (1 - y) (1 - z) +
	V010 (1 - x) y (1 - z) +
	V001 (1 - x) (1 - y) z +
	V101 x (1 - y) z +
	V011 (1 - x) y z +
	V110 x y (1 - z) +
	V111 x y z
	*/
	v4f temp;

	v4freset(out);

	v4fmul(&temp, p[0][0][0], (1 - f.x) * (1 - f.y) * (1 - f.z));
	v4fadd(out, *out, temp);
	v4fmul(&temp, p[1][0][0], (f.x) * (1 - f.y) * (1 - f.z));
	v4fadd(out, *out, temp);
	v4fmul(&temp, p[0][1][0], (1 - f.x) * (f.y) * (1 - f.z));
	v4fadd(out, *out, temp);
	v4fmul(&temp, p[0][0][1], (1 - f.x) * (1 - f.y) * (f.z));
	v4fadd(out, *out, temp);
	v4fmul(&temp, p[1][0][1], (f.x) * (1 - f.y) * (f.z));
	v4fadd(out, *out, temp);
	v4fmul(&temp, p[0][1][1], (1 - f.x) * (f.y) * (f.z));
	v4fadd(out, *out, temp);
	v4fmul(&temp, p[1][1][0], (f.x) * (f.y) * (1 - f.z));
	v4fadd(out, *out, temp);
	v4fmul(&temp, p[1][1][1], (f.x) * (f.y) * (f.z));
	v4fadd(out, *out, temp);
}

v3f rot3f(v3f v, float rad, float x, float y, float z)
{
	v3f newv;
	float costheta, sintheta;
	costheta = (float)cos(rad);
	sintheta = (float)sin(rad);

	newv.x = (costheta + (1 - costheta) * x * x)		 *v.x;
	newv.x += ((1 - costheta) * x * y - z * sintheta)	 *v.y;
	newv.x += ((1 - costheta) * x * z + y * sintheta)	 *v.z;

	newv.y = ((1 - costheta) * x * y + z * sintheta)	 *v.x;
	newv.y += (costheta + (1 - costheta) * y * y)		 *v.y;
	newv.y += ((1 - costheta) * y * z - x * sintheta)	 *v.z;

	newv.z = ((1 - costheta) * x * z - y * sintheta)	 *v.x;
	newv.z += ((1 - costheta) * y * z + x * sintheta)	 *v.y;
	newv.z += (costheta + (1 - costheta) * z * z)		 *v.z;

	return newv;
}

//brane node of mesh
struct brode
{
	v4f p;
};

typedef struct brode brode;

static int	brosz[3];
#define BROSZ	brosz
brode ***g_brode;
//MS3DModel g_md;

void to4(v3f in, v4f *out)
{
	int nearp[3][2];	//nearest grid point coordinate component / index [x/y/z][lower/higher] 
	v3f nearf;	//mid-way between nearest grid point fraction
	v4f nearc[2][2][2];	//nearest grid point 4-space coordinate	[x][y][z]
	int x, y, z;

	while (in.x < 0)
		in.x += BROSZ[0];
	while (in.y < 0)
		in.y += BROSZ[1];
	while (in.z < 0)
		in.z += BROSZ[2];

	nearf.x = in.x - (float)(int)in.x;
	nearf.y = in.y - (float)(int)in.y;
	nearf.z = in.z - (float)(int)in.z;

	while (in.x >= BROSZ[0])
		in.x -= BROSZ[0];
	while (in.y >= BROSZ[1])
		in.y -= BROSZ[1];
	while (in.z >= BROSZ[2])
		in.z -= BROSZ[2];

	nearp[0][0] = (int)in.x;
	nearp[0][1] = ((int)in.x + 1) % BROSZ[0];
	nearp[1][0] = (int)in.y;
	nearp[1][1] = ((int)in.y + 1) % BROSZ[1];
	nearp[2][0] = (int)in.z;
	nearp[2][1] = ((int)in.z + 1) % BROSZ[2];

	for (x = 0; x < 2; ++x)
	{
		for (y = 0; y < 2; ++y)
		{
			for (z = 0; z < 2; ++z)
			{
				//printf("nearp[0,1,2][%d,%d,%d] = %d,%d,%d\r\n",
				///	x, y, z, nearp[0][x], nearp[1][y], nearp[2][z]);
				nearc[x][y][z] = g_brode[nearp[0][x]][nearp[1][y]][nearp[2][z]].p;
				//printf("nearc[%d][%d][%d]=%f,%f,%f,%f\r\n",
				//	x, y, z,
				//	nearc[x][y][z].x,
				//	nearc[x][y][z].y,
				//	nearc[x][y][z].z,
				//	nearc[x][y][z].w);
			}
		}
	}

	trilinear(nearf, nearc, out);
}


void refract(float x, float y, float z, float dx, float dy, float dz, float d, v3f* d3);

void refract2(float x, float y, float z, float dx, float dy, float dz, float d, v3f* d3)
{
	//refract(x, y, z, dx, dy, dz, d, d3);
	//return;
	//d3->x = dx/d;
	//d3->y = dy/d;
	//d3->z = dz/d;
	//return;

	float dens[3][2];	//[x/y/z][+/-]
	v4f v4[4];
	v3f v3[5];
	float dx2 = 0, dy2 = 0, dz2;
	float d2;
	float s[10];
	int i;

	//"z"
	v3[0].x = dx;
	v3[0].y = dy;
	v3[0].z = dz;
	//v3[0] = norm3f(v3[0]);

	v3[1].x = 0;
	v3[1].y = 0;
	v3[1].z = 0;

	v3[2].x = x;
	v3[2].y = y;
	v3[2].z = z;

	v3[3].x = 0;
	v3[3].y = 0;
	v3[3].z = 0;

	v3fadd(&v3[4], v3[0], v3[2]);
	to4(v3[4], &v4[3]);

	if (fabs(dx) < 1)
		v3[1].x = 1;
	else if (fabs(dy) < 1)
		v3[1].y = 1;
	else
		v3[1].z = 1;

	v3[1] = cross3f(v3[1], v3[0]);	//side dir 1 "x"
	v3[3] = cross3f(v3[1], v3[0]);	//side dir 2 "y"

	v3[1] = norm3f(v3[1]);
	v3[3] = norm3f(v3[3]);

	//printf("v3[0]=%f,%f,%f,   dp %f   %f,%f,%f\r\n", v3[0].x, v3[0].y, v3[0].z, d, dx/d, dy/d, dz/d);
	//printf("v3[1]=%f,%f,%f,   x1\r\n", v3[1].x, v3[1].y, v3[1].z);
	//printf("v3[3]=%f,%f,%f,   y1\r\n", v3[3].x, v3[3].y, v3[3].z);
	//printf("v3[2]=%f,%f,%f,   p\r\n", v3[2].x, v3[2].y, v3[2].z);

	v3fmul(&v3[1], v3[1], d);
	v3fmul(&v3[3], v3[3], d);

	//printf("v3[1]=%f,%f,%f,   x1\r\n", v3[1].x, v3[1].y, v3[1].z);
	//printf("v3[3]=%f,%f,%f,   y1\r\n", v3[3].x, v3[3].y, v3[3].z);


	to4(v3[2], &v4[0]);
	//printf("v3[2]=%f,%f,%f xyz=>v4[0]=%f,%f,%f,%f\r\n", v3[2].x, v3[2].y, v3[2].z, v4[0].x, v4[0].y, v4[0].z, v4[0].w);

	//"x"+
	v3fadd(&v3[4], v3[2], v3[1]);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[0], v4[1]);
	dens[0][0] = mag4f(v4[2]);
	//printf("v3[4]=%f,%f,%f x+=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,		dens[0][0]);
	//"x+"+z
	v3fadd(&v3[4], v3[0], v3[4]);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[3], v4[1]);
	dens[0][0] = mag4f(v4[2]) / dens[0][0];
	//printf("v3[4]=%f,%f,%f x++z=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[0][0]);

	//"x"-
	v3fsub(&v3[4], v3[2], v3[1]);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[0], v4[1]);
	dens[0][1] = mag4f(v4[2]);
	//printf("v3[4]=%f,%f,%f x-=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[0][1]);
	//"x-"+z
	v3fadd(&v3[4], v3[0], v3[4]);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[3], v4[1]);
	dens[0][1] = mag4f(v4[2]) / dens[0][1];
	//printf("v3[4]=%f,%f,%f x-+z=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[0][1]);

	//"y"+
	v3fadd(&v3[4], v3[2], v3[3]);
	//printf("v3[3]=%f,%f,%f,\r\n", v3[3].x, v3[3].y, v3[3].z);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[0], v4[1]);
	dens[1][0] = mag4f(v4[2]);
	//printf("v3[4]=%f,%f,%f y+=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[1][0]);
	//"y+"+z
	v3fadd(&v3[4], v3[0], v3[4]);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[3], v4[1]);
	dens[1][0] = mag4f(v4[2]) / dens[1][0];
	//printf("v3[4]=%f,%f,%f y++z=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[1][0]);

	//"y"-
	v3fsub(&v3[4], v3[2], v3[3]);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[0], v4[1]);
	dens[1][1] = mag4f(v4[2]);
	//printf("v3[4]=%f,%f,%f y-=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[1][1]);
	//"y-"+z
	v3fadd(&v3[4], v3[0], v3[4]);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[3], v4[1]);
	dens[1][1] = mag4f(v4[2]) / dens[1][1];
	//printf("v3[4]=%f,%f,%f y-+z=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[1][1]);

	//"z"+
	//v3fadd(&v3[4], v3[2], v3[0]);
	//to4(v3[4], &v4[1]);
	//v4fsub(&v4[2], v4[0], v4[1]);
	//dens[2][0] = mag4f(v4[2]);

	//"z"-
	//v3fsub(&v3[4], v3[2], v3[0]);
	//to4(v3[4], &v4[1]);
	//v4fsub(&v4[2], v4[0], v4[1]);
	//dens[2][1] = mag4f(v4[2]);

	//dx2 = (dens[0][1] - dens[0][0]) / 6.0f;	// /1.0f because [-0.5,0.5]=1 not [-1,1]=3     eg [-1,1]with 1,1+1 / 3 = 1     [-0.5,0.5]with 1,1+1.......
	//d(x)=d(min)+(d(max)-d(min))/(max-min) * x          x' = min+(x-min)/(max-min) 
	for (i = 0; i < 10; ++i)
	{
		s[i] = dens[0][1] + (dens[0][0] - dens[0][1]) / 9.0f;
		dx2 += s[i];
	}
	dx2 = (s[4] + s[5]) / dx2 / 2.0f;

	//x'=(d(x-1)-d(x+1))/3+x, y' = y + 1 / d(y + 1)
	//dy2 = (dens[1][1] - dens[1][0]) / 6.0f; 
	for (i = 0; i < 10; ++i)
	{
		s[i] = dens[1][1] + (dens[1][0] - dens[1][1]) / 9.0f;
		dy2 += s[i];
	}
	dy2 = (s[4] + s[5]) / dy2 / 2.0f;

	//dz2 = (dens[2][0] / dens[2][1]);

	//d2 = sqrtf(dx2*dx2 + dy2*dy2);

	//dx2 /= d2;
	//dy2 /= d2;

	//1            2      3         4        5            6          7           8        9                    sample points 1->9  starting coordinate space squares
	//2          1.875   1.75     1.625    1.5          1.375       1.25       1.125      1          13.5      linear interp, sum at end, of proper space per coordinate space squares in new place
	//.148148    .13888  .12962   .12037   .11111      .101851     .092592    .08333     0.074074              ratio of total
	//1.333332   1.24992  1.16658  1.08333  .9999      0.916659    0.833328   0.749997   0.666666    8.999712  ratio times max (9) and sum at end, so that the coordinate squares push out from left
	//1.333332   2.583252 3.749832 4.833162 5.833062   6.749721    7.583049   8.333046   8.999712              adding up the above coordinate squares from left to right, to give their position
	//1.333374   2.583334         dividing above coordinate squares by max (8.999712 or actually 9) 
	//5=>5.833062
	//do with different interval of samples to check


	//1           2          3           4       5        6     7
	//2        1.833       1.66        1.5     1.33    1.166    1     10.4999995
	//.1904    .1746      .1587      .1428    .1269    .1111 .0952
	//         0.365      .5237      .6665    .7934    .9045 .9997
	//1.3328  2.555       3.6659    4.6655   5.5538   6.3378

	//d2 = d * dz2;
	//d2 = sqrtf(dx*dx + dy*dy + dz*dz);
	//d2 = sqrtf(dx2*dx2 + dy2*dy2 + dz2*dz2);
	//dx2 /= d2;
	//dy2 /= d2;
	//dz2 /= d2;
	//dx2 *= d;

	//v3fmul(&v3[0], v3[0], dz2);
	v3fmul(&v3[1], v3[1], dx2);
	v3fmul(&v3[3], v3[3], dy2);

	v3fadd(&v3[0], v3[0], v3[1]);
	v3fadd(&v3[0], v3[0], v3[3]);

	v3[0] = norm3f(v3[0]);
	v3fmul(&v3[0], v3[0], d);

	//"z"+
	v3fadd(&v3[4], v3[2], v3[0]);
	to4(v3[4], &v4[1]);
	v4f aaa = v4[1];
	v3f aaaa = v3[4];
	v4fsub(&v4[2], v4[0], v4[1]);
	dens[2][0] = mag4f(v4[2]);
	//printf("v3[4]=%f,%f,%f +z=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[2][0]);

	//"z"-
	v3fsub(&v3[4], v3[2], v3[0]);
	to4(v3[4], &v4[1]);
	v4fsub(&v4[2], v4[0], v4[1]);
	dens[2][1] = mag4f(v4[2]);
	//printf("v3[4]=%f,%f,%f -z=>v4[1]=%f,%f,%f,%f  v4[2]=%f,%f,%f,%f mag%f\r\n", v3[4].x, v3[4].y, v3[4].z, v4[1].x, v4[1].y, v4[1].z, v4[1].w, v4[2].x, v4[2].y, v4[2].z, v4[2].w,	dens[2][1]);

	//dz2 = dens[2][0] / dens[2][1];
	dz2 = dens[2][1] / dens[2][0];

	v3fmul(d3, v3[0], dz2);

	//printf("d3%f,%f,%fl%f,%f,%f\r\n", dx / d, dy / d, dz / d, d3->x, d3->y, d3->z);

	//d3->x = dx / d * 0.95f + 0.05f * d3->x;
	//d3->y = dy / d * 0.95f + 0.05f * d3->y;
	//d3->z = dz / d * 0.95f + 0.05f * d3->z;

	//printf("[%f,%f][%f,%f,%f,%f][%f,%f,%f,%f][%f,%f,%f]ddd%fasdasd%flsjlasdkj%fasldjasd(%f,%f,%f)(%f,%f,%f)", dens[2][0], dens[2][1], v4[0].x, v4[0].y, v4[0].z, v4[0].w,	aaa.x, aaa.y, aaa.z, aaa.w, aaaa.x, aaaa.y, aaaa.z, d*dz2, d, dz2, d3->x, d3->y, d3->z, v3[0].x, v3[0].y, v3[0].z);
}

double angle3f(v3f v1, v3f v2)
{
	float dot, vmag;
	double angle;

	dot = dot3f(v1, v2);
	vmag = mag3f(v1) * mag3f(v2);
	angle = acos(dot / vmag);

	//if (_isnan(angle))
#define ISNAN(a)	(a!=a)
	if(ISNAN(angle))
		return 0;

	return(angle);
}

void refract(float x, float y, float z, float dx, float dy, float dz, float d, v3f* d3)
{
	float tdr;
	int x1, x2, y1, y2, z1, z2;
	float mx, my, mz;

	float kx1, ky1, kz1, kw1;
	float kx2, ky2, kz2, kw2;

	v4f densc[4];
	float dens[2];
	float dx2, dy2, dz2;
	float d2;
	v3f in;
	v4f off;
	float trav[2];
	float ddens[3];
	v4f ddensc[3][2];
	v3f n;
	float ran[2];
	//v3f d3;
	v3f bv;
	v3f ax;
	float nd;
	float ri[2];
	float rir;
	float jump = 1;

	//calculate light traversal time through density, compared to previous position's density

	//printf("dx=%f,y=%f,z=%f,=%f\r\n", dx, dy, dz, d);
	//if (d < 0.00001f)
	//	system("pause");

	dx2 = dx / d;
	dy2 = dy / d;
	dz2 = dz / d;
	d2 = 1;

	do
	{

		//printf("dx2=%f,y2=%f,z2=%f,2=%f\r\n", dx2, dy2, dz2, d2);

		in.x = x - dx2 / jump;
		in.y = y - dy2 / jump;
		in.z = z - dz2 / jump;
		to4(in, &densc[0]);

		//printf("in1 %f,%f,%f \r\n", in.x, in.y, in.z);

		in.x = x + dx2 / jump;
		in.y = y + dy2 / jump;
		in.z = z + dz2 / jump;
		to4(in, &densc[1]);

		//printf("in2 %f,%f,%f \r\n", in.x, in.y, in.z);

		v4fsub(&off, densc[1], densc[0]);
		dens[0] = sqrtf(off.x*off.x + off.y*off.y + off.z*off.z + off.w*off.w);

		//printf("densc[0] %f,%f,%f     densc[1] %f,%f,%f \r\n",
		////	densc[0].x, densc[0].y, densc[0].z,
		//	densc[1].y, densc[1].y, densc[1].z);
		//printf("dens%f\r\n", dens[0]);

		//if (dens[0] < 0.00001f)
		//	system("pause");

		trav[0] = dens[0] * jump;
	} while (trav[0] < 0.001f && (jump /= 10.0f));

	//slow down light by density of space
	dx = dx2 / trav[0];
	dy = dy2 / trav[0];
	dz = dz2 / trav[0];
	d = jump / trav[0];

	//second point's density, from two more 4-coordinates
	in.x = x - dx2 / jump + dx;
	in.y = y - dy2 / jump + dy;
	in.z = z - dz2 / jump + dz;
	to4(in, &densc[2]);
	in.x = x + dx2 / jump + dx;
	in.y = y + dy2 / jump + dy;
	in.z = z + dz2 / jump + dz;
	to4(in, &densc[3]);
	v4fsub(&off, densc[3], densc[2]);
	dens[1] = sqrtf(off.x*off.x + off.y*off.y + off.z*off.z + off.w*off.w);
	trav[1] = dens[1] * jump;

	//find surface normal of increasing density

	//x
	in.x = x;
	in.y = y;
	in.z = z;
	to4(in, &ddensc[0][0]);
	in.x = x + fabs(dx2) / jump;
	to4(in, &ddensc[0][1]);
	v4fsub(&off, ddensc[0][1], ddensc[0][0]);
	ddens[0] = sqrtf(off.x*off.x + off.y*off.y + off.z*off.z + off.w*off.w);
	//x2
	in.x = x - fabs(dx2) / jump;
	to4(in, &ddensc[0][0]);
	in.x = x;
	to4(in, &ddensc[0][1]);
	v4fsub(&off, ddensc[0][1], ddensc[0][0]);
	ddens[0] -= sqrtf(off.x*off.x + off.y*off.y + off.z*off.z + off.w*off.w);

	//y
	in.x = x;
	in.y = y;
	in.z = z;
	to4(in, &ddensc[1][0]);
	in.y = y + fabs(dy2) / jump;
	to4(in, &ddensc[1][1]);
	v4fsub(&off, ddensc[1][1], ddensc[1][0]);
	ddens[1] = sqrtf(off.x*off.x + off.y*off.y + off.z*off.z + off.w*off.w);
	//y2
	in.y = y - fabs(dy2) / jump;
	to4(in, &ddensc[1][0]);
	in.y = y;
	to4(in, &ddensc[1][1]);
	v4fsub(&off, ddensc[1][1], ddensc[1][0]);
	ddens[1] -= sqrtf(off.x*off.x + off.y*off.y + off.z*off.z + off.w*off.w);

	//z
	in.x = x;
	in.y = y;
	in.z = z;
	to4(in, &ddensc[2][0]);
	in.z = z + fabs(dz2) / jump;
	to4(in, &ddensc[2][1]);
	v4fsub(&off, ddensc[2][1], ddensc[2][0]);
	ddens[2] = sqrtf(off.x*off.x + off.y*off.y + off.z*off.z + off.w*off.w);
	//z2
	in.z = z - fabs(dz2) / jump;
	to4(in, &ddensc[2][0]);
	in.z = z;
	to4(in, &ddensc[2][1]);
	v4fsub(&off, ddensc[2][1], ddensc[2][0]);
	ddens[2] -= sqrtf(off.x*off.x + off.y*off.y + off.z*off.z + off.w*off.w);

	//normal of surface
	n.x = -ddens[0];
	n.y = -ddens[1];
	n.z = -ddens[2];
	nd = sqrtf(n.x*n.x + n.y*n.y + n.z*n.z);
	n.x /= nd;
	n.y /= nd;
	n.z /= nd;

	//get the refractive index
	ri[0] = trav[0];
	ri[1] = trav[1];
	rir = ri[0] / ri[1];
	//get angle of incidence upon surface
	d3->x = dx;
	d3->y = dy;
	d3->z = dz;
	ran[0] = angle3f(*d3, n);
	//refractive angle
	//n1/n2 = sino2/sino1
	//sino2 = sino1 n1 / n2
	//o2 = asin( sino1 n1 / n2 )
	//printf("asin sin(%f) %f = %f => %f\r\n", ri[0], rir, sin(ran[0]) * rir, asin(sin(ran[0]) * rir));

	if (ri[1] >= ri[0])
		ran[1] = asin(sin(ran[0]) * rir);
	else
	{
		//n1/n2 = sino2/sino1
		//sino1 = sino2 n2 / n1
		//o1 = asin( sino2 n2 / n1 )
		ran[1] = asin(sin(ran[0]) / rir);
	}

	//refract by using an axis-angle rotation of the direction vector
	bv = cross3f(*d3, n);
	//ax = cross3f(bv, d3);
	ax = bv;
	//printf("rot %f,%f,%f,%f\r\n%f,%f,%f\r\n", ran[1], ax.x, ax.y, ax.z, d3->x, d3->y, d3->z);
	*d3 = rot3f(*d3, ran[1], ax.x, ax.y, ax.z);

	//check for intersections between [x,y,z] and [x+d3.x,y+d3.y,z+d3.z]
	//...
}

#if 0
int trace(float x, float y, float z, float dx, float dy, float dz, float d, unsigned char* r, unsigned char* g, unsigned char *b, v3f* d3)
{
	v3f line[2];
	v3f tri[3];
	int i;
	int j;
	v3f n;
	float o;
	v3f v;

	refract(x, y, z, dx, dy, dz, d, d3);

	//check for intersections between [x,y,z] and [x+d3.x,y+d3.y,z+d3.z]
	//...

	line[0].x = x;
	line[0].y = y;
	line[0].z = z;
	line[1].x = x + d3->x;
	line[1].y = y + d3->y;
	line[1].z = z + d3->z;

	for (i = 0; i < g_md.m_numTriangles; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			tri[j] = *(v3f*)g_md.m_pVertices[g_md.m_pTriangles[i].m_vertexIndices[j]].m_location;
		}

		if (intpg(tri, line, &n, &o, &v, 3))
		{
			*r = 0;
			*g = 0;
			*b = 0;
			return 1;
		}
	}

	return 0;
}
#endif

int main()
{
	int i;
	v3f a;
	int x, y, z;
	char c;
	FILE *f;
	v3f b;
	v3f d;
	v3f e;
	float dd;
	v3f g;

	i = 0;

asking:

	printf("load data? ");
	scanf("%c", &c);

	if (c == 'y')
		goto loading;
	else if (c == 'n')
		goto reading;
	else
		goto asking;

loading:

	f = fopen("data.dat", "rb");
	fread(&BROSZ, sizeof(int), 3, f);

	g_brode = (brode***)malloc(sizeof(brode**)*BROSZ[0]);
	for (x = 0; x < BROSZ[0]; ++x)
	{
		g_brode[x] = (brode**)malloc(sizeof(brode*)*BROSZ[1]);
		for (y = 0; y < BROSZ[1]; ++y)
		{
			g_brode[x][y] = (brode*)malloc(sizeof(brode)*BROSZ[2]);
			for (z = 0; z < BROSZ[2]; ++z)
			{
				fread(&g_brode[x][y][z].p, sizeof(float), 4, f);
				printf("%d,%d,%dx: %f", x, y, z, g_brode[x][y][z].p.x);
				printf("%d,%d,%dy: %f", x, y, z, g_brode[x][y][z].p.y);
				printf("%d,%d,%dz: %f", x, y, z, g_brode[x][y][z].p.z);
				printf("%d,%d,%dw: %f", x, y, z, g_brode[x][y][z].p.w);
			}
		}
	}

	fclose(f);

	goto starting;

reading:

	printf("xn: ");
	scanf("%d", &BROSZ[0]);
	printf("yn: ");
	scanf("%d", &BROSZ[1]);
	printf("zn: ");
	scanf("%d", &BROSZ[2]);

	g_brode = (brode***)malloc(sizeof(brode**)*BROSZ[0]);
	for (x = 0; x < BROSZ[0]; ++x)
	{
		g_brode[x] = (brode**)malloc(sizeof(brode*)*BROSZ[1]);
		for (y = 0; y < BROSZ[1]; ++y)
		{
			g_brode[x][y] = (brode*)malloc(sizeof(brode)*BROSZ[2]);
		}
	}

	for (z = 0; z < BROSZ[2]; ++z)
	{
		for (y = 0; y < BROSZ[1]; ++y)
		{
			for (x = 0; x < BROSZ[0]; ++x)
			{
				printf("%d,%d,%dx: ", x, y, z);
				scanf("%f", &g_brode[x][y][z].p.x);
				printf("%d,%d,%dy: ", x, y, z);
				scanf("%f", &g_brode[x][y][z].p.y);
				printf("%d,%d,%dz: ", x, y, z);
				scanf("%f", &g_brode[x][y][z].p.z);
				printf("%d,%d,%dw: ", x, y, z);
				scanf("%f", &g_brode[x][y][z].p.w);
			}
		}
	}

saving:

	printf("save data? ");
	scanf("%c", &c);

	if (c == 'y')
		goto save;
	else if (c == 'n')
		goto starting;
	else
		goto saving;

save:

	f = fopen("data.dat", "wb");
	fwrite(&BROSZ, sizeof(int), 3, f);

	for (x = 0; x < BROSZ[0]; ++x)
	{
		for (y = 0; y < BROSZ[1]; ++y)
		{
			for (z = 0; z < BROSZ[2]; ++z)
			{
				fwrite(&g_brode[x][y][z].p, sizeof(float), 4, f);
			}
		}
	}

	fclose(f);

starting:

	while (1)
	{
		printf("xf%d: ", i);
		scanf("%f", &a.x);
		printf("yf%d: ", i);
		scanf("%f", &a.y);
		printf("zf%d: ", i);
		scanf("%f", &a.z);
		printf("xt%d: ", i);
		scanf("%f", &b.x);
		printf("yt%d: ", i);
		scanf("%f", &b.y);
		printf("zt%d: ", i);
		scanf("%f", &b.z);
		printf("d%d: ", i);
		scanf("%f", &dd);

		v3fsub(&d, b, a);
		refract(a.x, a.y, a.z, d.x, d.y, d.z, dd, &e);
		v3fadd(&g, e, b);

		printf("refract: d(t,f)'=%f,%f,%f   t'=%f,%f,%f\r\n\r\n", e.x, e.y, e.z, g.x, g.y, g.z);

		i++;
	}


#if 0
	g_camf.pos.x = 0;
	g_camf.pos.y = 0;
	g_camf.pos.z = -1;
	g_camf.strafe.x = 1;
	g_camf.strafe.y = 0;
	g_camf.strafe.z = 0;
	g_camf.up.x = 0;
	g_camf.up.y = 1;
	g_camf.up.z = 0;
	g_camf.view.x = -10;
	g_camf.view.y = 0;
	g_camf.view.z = 0;

#define X	0
#define Y	1
#define Z	2
#define NEARP	0
#define FARP		1
#define POS		0
#define NEG		1

	vvvv[X][NEARP][POS] = g_camf.strafe;
	v3fsub(&vvvv[Z][NEARP][POS], g_camf.view, g_camf.pos);
	vvvv[Z][NEARP][POS] = norm3f(vvvv[Z][NEARP][POS]);
	vvvv[Y][NEARP][POS] = norm3f(cross3f(vvvv[X][NEARP][POS], vvvv[Z][NEARP][POS]));

	v3fmul(&vvvv[X][FARP][POS], vvvv[X][NEARP][POS], MAX_DISTANCE - MIN_DISTANCE);
	v3fmul(&vvvv[Y][FARP][POS], vvvv[Y][NEARP][POS], MAX_DISTANCE - MIN_DISTANCE);
	v3fmul(&vvvv[Z][FARP][POS], vvvv[Z][NEARP][POS], MAX_DISTANCE - MIN_DISTANCE);

	v3fmul(&vvvv[X][NEARP][NEG], vvvv[X][NEARP][POS], -1);
	v3fmul(&vvvv[Y][NEARP][NEG], vvvv[Y][NEARP][POS], -1);
	v3fmul(&vvvv[Z][NEARP][NEG], vvvv[Z][NEARP][POS], -1);
	v3fmul(&vvvv[X][FARP][NEG], vvvv[X][FARP][POS], -1);
	v3fmul(&vvvv[Y][FARP][NEG], vvvv[Y][FARP][POS], -1);
	v3fmul(&vvvv[Z][FARP][NEG], vvvv[Z][FARP][POS], -1);

	v3fadd(&pv[0], vvvv[X][FARP][NEG], vvvv[Y][FARP][POS]);
	v3fadd(&pv[0], pv[0], vvvv[Z][FARP][POS]);
	v3fadd(&pv[1], vvvv[X][FARP][POS], vvvv[Y][FARP][POS]);
	v3fadd(&pv[1], pv[1], vvvv[Z][FARP][POS]);
	v3fadd(&pv[2], vvvv[X][FARP][POS], vvvv[Y][FARP][NEG]);
	v3fadd(&pv[2], pv[2], vvvv[Z][FARP][POS]);
	v3fadd(&pv[3], vvvv[X][FARP][NEG], vvvv[Y][FARP][NEG]);
	v3fadd(&pv[3], pv[3], vvvv[Z][FARP][POS]);
	v3fadd(&pv[4], vvvv[X][NEARP][NEG], vvvv[Y][NEARP][POS]);
	v3fadd(&pv[4], pv[4], vvvv[Z][NEARP][POS]);
	v3fadd(&pv[5], vvvv[X][NEARP][POS], vvvv[Y][NEARP][POS]);
	v3fadd(&pv[5], pv[5], vvvv[Z][NEARP][POS]);
	v3fadd(&pv[6], vvvv[X][NEARP][POS], vvvv[Y][NEARP][NEG]);
	v3fadd(&pv[6], pv[6], vvvv[Z][NEARP][POS]);
	v3fadd(&pv[7], vvvv[X][NEARP][NEG], vvvv[Y][NEARP][NEG]);
	v3fadd(&pv[7], pv[7], vvvv[Z][NEARP][POS]);

#undef X
#undef Y
#undef Z
#undef NEARP
#undef FARP
#undef POS
#undef NEG

	v3fadd(&pv[0], pv[0], g_camf.pos);
	v3fadd(&pv[1], pv[1], g_camf.pos);
	v3fadd(&pv[2], pv[2], g_camf.pos);
	v3fadd(&pv[3], pv[3], g_camf.pos);
	v3fadd(&pv[4], pv[4], g_camf.pos);
	v3fadd(&pv[5], pv[5], g_camf.pos);
	v3fadd(&pv[6], pv[6], g_camf.pos);
	v3fadd(&pv[7], pv[7], g_camf.pos);

	//glUniform3fv(s->slot[SSLOT_P], 8, pv);

#if 0
	v[0].x = 2.0f;
	v[0].y = 5.0f;
	v[0].z = -50.0f;
	v[1].x = 7.0f;
	v[1].y = 1.0f;
	v[1].z = -50.0f;
	v[2].x = 2.0f;
	v[2].y = 1.0f;
	v[2].z = -50.0f;
#endif

	glUniform4f(s->slot[SSLOT_COLOR], 1, 1, 1, 1);

	for (x = 0; x < NPX; ++x)
	{
		for (y = 0; y < NPX; ++y)
		{
			for (z = 0; z < NPX; ++z)
			{
				v[0] = p[z*NPX*NPX + y*NPX + x];
				//xy = toxy(v[0], WX, WY, pv);
				//im[0] = 255;

				//xy.x += 1;

#if 0
				if (xy.x >= 0 && xy.x < WX && xy.y >= 0 && xy.y < WY)
				{
					//fprintf(g_applog, "xy.x%f %f\r\n", xy.x, xy.y);
					//fflush(g_applog);
					im[((int)xy.x%WX) * 3 + 3 * WX * ((int)xy.y%WY)] = 255;
					im[((int)xy.x%WX) * 3 + 3 * WX * ((int)xy.y%WY) + 1] = 255;
					im[((int)xy.x%WX) * 3 + 3 * WX * ((int)xy.y%WY) + 2] = 255;
				}
#endif
				v3fsub(&v[0], v[0], cv);

				v3fmul(&v[0], v[0], 10.0f);

				if (z + 1 < NPX)
				{
					v[1] = p[(z + 1)*NPX*NPX + y*NPX + x];
					//v3fsub(&v[1], v[1], cv);
					//v3fmul(&v[1], v[1], 10.0f);
					linexy(v[0], v[1], WX, WY, im, 1, 0, 0);
					//glVertexPointer(3, GL_FLOAT, 0, v);
					//glDrawArrays(GL_LINES, 0, 2);
				}
				if (x + 1 < NPX)
				{
					v[1] = p[z*NPX*NPX + y*NPX + x + 1];
					//v3fsub(&v[1], v[1], cv);
					//v3fmul(&v[1], v[1], 10.0f);
					linexy(v[0], v[1], WX, WY, im, 0, 1, 0);
					//glVertexPointer(3, GL_FLOAT, 0, v);
					//glDrawArrays(GL_LINES, 0, 2);
				}
				if (y + 1 < NPX)
				{
					v[1] = p[z*NPX*NPX + (y + 1)*NPX + x];
					//v3fsub(&v[1], v[1], cv);
					//v3fmul(&v[1], v[1], 10.0f);
					linexy(v[0], v[1], WX, WY, im, 0, 0, 1);
					//glVertexPointer(3, GL_FLOAT, 0, v);
					//glDrawArrays(GL_LINES, 0, 2);
				}
			}
		}
	}
#endif

	return 0;
}


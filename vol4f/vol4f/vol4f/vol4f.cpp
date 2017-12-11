// vol4f.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>

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
	v3f ab;
	v3f ac;
	float f;
	v3fsub(&ab, a, b);
	v3fsub(&ac, a, c);
	f = mag3f(ab);
	return 0.5 * f * (1.0 - dot3f(ab, ac) / f)*f;
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

v3f toxy(v3f vi, float wx, float wy, v3f p[8])
{
	float v[12];
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
	vi.x = wx * (v[8] + v[9]) / (v[8] + v[9] + v[4] + v[5]);
	vi.y = wy * (v[6] + v[7]) / (v[6] + v[7] + v[2] + v[3]);
	vi.z = ((v[10] + v[11]) / (v[0] + v[1] + v[10] + v[11]));
	return vi;
}

int main()
{
	int i;
	v4f a[4];

	i = 0;

	while (1)
	{
		for (i = 0; i < 4; i++)
		{
			printf("x%d: ", i);
			scanf("%f", &a[i].x);
			printf("y%d: ", i);
			scanf("%f", &a[i].y);
			printf("z%d: ", i);
			scanf("%f", &a[i].z);
			printf("w%d: ", i);
			scanf("%f", &a[i].w);
		}

		printf("vol4f: %f\r\n\r\n", vol4f(a[0], a[1], a[2], a[3]));
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


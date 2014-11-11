/* FAST code from http://svr-www.eng.cam.ac.uk/~er258/work/fast.html. See LICENCE */

#ifndef FAST_H
#define FAST_H

#include "../../../util/util/optimisation_attributes.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { int x, y; } xy;
typedef unsigned char byte;

int fast9_corner_score(const byte* p, const int pixel[], int bstart) HOT;
int fast10_corner_score(const byte* p, const int pixel[], int bstart) COLD;
int fast11_corner_score(const byte* p, const int pixel[], int bstart) COLD;
int fast12_corner_score(const byte* p, const int pixel[], int bstart) COLD;

xy* fast9_detect(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners) HOT;
xy* fast10_detect(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners) COLD;
xy* fast11_detect(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners) COLD;
xy* fast12_detect(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners) COLD;

int* fast9_score(const byte* i, int stride, xy* corners, int num_corners, int b) HOT;
int* fast10_score(const byte* i, int stride, xy* corners, int num_corners, int b) COLD;
int* fast11_score(const byte* i, int stride, xy* corners, int num_corners, int b) COLD;
int* fast12_score(const byte* i, int stride, xy* corners, int num_corners, int b) COLD;


xy* fast9_detect_nonmax(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners) HOT;
xy* fast10_detect_nonmax(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners) COLD;
xy* fast11_detect_nonmax(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners) COLD;
xy* fast12_detect_nonmax(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners) COLD;

xy* nonmax_suppression(const xy* corners, const int* scores, int num_corners, int* ret_num_nonmax) HOT;

#ifdef __cplusplus
}
#endif

#endif

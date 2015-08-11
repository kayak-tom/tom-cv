#ifndef RECONSTRUCT3D_H
#define RECONSTRUCT3D_H

#include "newCamera.h"
#include <geom/levMarNumerical.h>

optional<C3dWorldPoint> reconstruct(const CWorldCamera &P, const CWorldCamera &Pp, const C2dImagePointPx &p, const C2dImagePointPx &pp, const double dReprojectionErrorPx, int bVerbose);
optional<const C3dWorldPoint> reconstructLM(const CWorldCamera &P, const CWorldCamera &Pp, const C2dImagePointPx &p, const C2dImagePointPx &pp);
void test3DRecon(const CWorldCamera & P1, const CWorldCamera & P2);

#endif // RECONSTRUCT3D_H

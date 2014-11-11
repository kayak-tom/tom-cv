#pragma once

#include "params/param.h"

PARAMCLASS(Patch)
	PARAMB(MONO_DESCRIPTOR, true, "Greyscale rather than RGB descriptors")
	PARAMB(NORMALISE, false, "Normalise intensity in each descriptor (so each has the same mean; for illumination invariance). Normalising image intensity probably works better.")
	PARAMB(ORIENT, false, "Rotation invariance: orient with strongest change in intensity (a bit like SURF).")
	PARAM(PATCH_SCALE, 1, 10, 4, "Patches are subsampled. Essentially the same as downsizing input images (except corner localisation will be better). Don't know what is best")
	//PARAM(PATCH_BLUR, 0, 10, 0, "A Gaussian blur can be applied to whole image (if > 0) to get rid of noise effects when subsampling.")
	{}

	CNumParam<bool> MONO_DESCRIPTOR, NORMALISE, ORIENT;
	CNumParam<int> PATCH_SCALE;//, PATCH_BLUR;
	enum eColourTypeSize {GREY=1, RG=2, RGB=3};
};

PARAMCLASS(PatchDescriptor)
	PARAM(PATCH_RAD, 1, 6, 5, "Patch will be square with sides 2*PATCH_RAD+1")
	PARAME(PATCH_COMP_METHOD, PatchEuclidFast, "Choose norm for comparing patches. Euclidean seems best. The squared val is used.")
	PARAM(PX_RADIUS, 1, 255, 50, "If 2 descriptors have average differnce of more than this many grey levels they are different.")
	CHILDCLASS(Patch, "Simplest descriptor: Just an image patch centred around corner")
	{}

	CNumParam<int> PATCH_RAD;
	MAKEENUMPARAM8(PATCH_COMP_METHOD, PatchEuclidFast, PatchL1Fast, PatchMaxDist, PatchCorrel, PatchL1, PatchEuclid, PatchEuclidParallel, PatchL1Parallel); /*Todo: consider CPatchDescriptorParams::ePatchCos, CPatchDescriptorParams::ePatchChiSquared, CPatchDescriptorParams::ePatchJeffereys??*/
	CNumParam<int> PX_RADIUS;
	MAKECHILDCLASS(Patch);

	int margin() const
	{
		const int SQUARE_RAD=Patch.PATCH_SCALE * PATCH_RAD;
		if(Patch.ORIENT)
			return (SQUARE_RAD*15)/10;
		else
			return SQUARE_RAD+1; //Todo: SQUARE_RAD should be ok
	}

	//PX_RADIUS should be the borderline average absolute pixel dist between 2 pixels in patches that are the same/different
	//radius() converts this to the distance between these patches using the current norm
	int radius() const
	{
		int AREA = (int)(Patch.MONO_DESCRIPTOR ? 1 : 2) * (int)sqr((2*(int)(PATCH_RAD))+1);

		switch(PATCH_COMP_METHOD)
		{
		case CPatchDescriptorParams::ePatchEuclid:
		case CPatchDescriptorParams::ePatchEuclidFast:
			return sqr(PX_RADIUS) * AREA;

		case CPatchDescriptorParams::ePatchEuclidParallel:
			return (sqr(PX_RADIUS) * AREA)/2; //Because quantising to 0..128

		case CPatchDescriptorParams::ePatchL1:
			return PX_RADIUS * AREA;

		case CPatchDescriptorParams::ePatchL1Parallel:
			return (PX_RADIUS * AREA)/2; //Because quantising to 0..128

		case CPatchDescriptorParams::ePatchMaxDist:
			return PX_RADIUS;

		default:
			THROW("Unhandled patch type")
		}
	}
};


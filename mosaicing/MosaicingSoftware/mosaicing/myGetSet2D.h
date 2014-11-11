#ifndef _MY_GET_SET_2D_
#define _MY_GET_SET_2D_

#include "util/exception.h"
#include "util/opencv.h"

namespace grc {
	/*!
	 * \brief Fast version of cvGet2D
	 *
	 * Gets a pointer to the first colour channel of an image at a given location.
	 * There is no bounds checking or other safety checks, but this should be very fast.
	 * This only works for 8-bit images (type IPL_DEPTH_8U).
	 * For single channel images this gets a pointer to the grey value.
	 *
	 * \param img The image to get pixel values from
	 * \param x   The x-coordinate of the pixel to get
	 * \param y   The y-coordinate of the pixel to get
	 *
	 * \returns A pointer to the first element of the image data for the given pixel
	 */
	static inline unsigned char myGet2D(IplImage *img, int x, int y) HOT PURE_FN HARD_INLINE;
	static inline unsigned char myGet2D(IplImage *img, int x, int y) {
		return img->imageData[y*img->widthStep + x*img->nChannels];
	}
	static inline unsigned char * myGet2D_ref(IplImage *img, int x, int y) HOT PURE_FN HARD_INLINE;
	static inline unsigned char * myGet2D_ref(IplImage *img, int x, int y) {
		return reinterpret_cast<unsigned char *>(img->imageData + y*img->widthStep + x*img->nChannels);
	}

	/*!
	 * \brief Fast version of cvGet2D
	 *
	 * Gets a pointer to the specified colour channel of an image at a given location.
	 * There is no bounds checking or other safety checks, but this should be very fast.
	 * This only works for 8-bit images (type IPL_DEPTH_8U).
	 *
	 * \param img The image to get pixel values from
	 * \param x   The x-coordinate of the pixel to get
	 * \param y   The y-coordinate of the pixel to get
	 * \param c   The (0-based) index of the colour channel to get
	 *
	 * \returns A pointer to the (c+1)th element of the image data for the given pixel
	 */
	static inline unsigned char myGet2D(IplImage *img, int x, int y, int c) HOT PURE_FN HARD_INLINE;
	static inline unsigned char myGet2D(IplImage *img, int x, int y, int c) {
		return img->imageData[y*img->widthStep + x*img->nChannels + c];
	}
	static inline unsigned char * myGet2D_ref(IplImage *img, int x, int y, int c) HOT PURE_FN HARD_INLINE;
	static inline unsigned char * myGet2D_ref(IplImage *img, int x, int y, int c) {
		return reinterpret_cast<unsigned char *>(img->imageData + y*img->widthStep + x*img->nChannels + c);
	}
	static inline int * myGet2D_int(IplImage *img, int x, int y) HOT PURE_FN HARD_INLINE;
	static inline int * myGet2D_int(IplImage *img, int x, int y) {
		if(IS_DEBUG) CHECK(img->nChannels != 4, "Wrong num of channels, need 4 for int cast")
		return reinterpret_cast<int *>(img->imageData + y*img->widthStep + x*4);
	}

	/*!
	 * \brief Fast version of cvSet2D
	 *
	 * Sets the first colour channel of an image at a given location.
	 * There is no bounds checking or other safety checks, but this should be very fast.
	 * This only works for 8-bit images (type IPL_DEPTH_8U).
	 * For single channel images this sets the grey value.
	 *
	 * \param img The image to get pixel values from
	 * \param x   The x-coordinate of the pixel to get
	 * \param y   The y-coordinate of the pixel to get
	 * \param value  The new colour value
	 */
	static inline void mySet2D(IplImage *img, int x, int y, int c, unsigned char value) HOT HARD_INLINE;
	static inline void mySet2D(IplImage *img, int x, int y, int c, unsigned char value) {
		img->imageData[y*img->widthStep + x*img->nChannels + c] = value;
	}
	/*static inline void mySet2D_int(IplImage *img, int x, int y, int value) HOT HARD_INLINE;
	static inline void mySet2D_int(IplImage *img, int x, int y, int value) {
		reinterpret_cast<int *>(img->imageData[y*img->widthStep + x*4] = value;
	}*/

	/*!
	 * \brief Fast version of cvSet2D
	 *
	 * Sets the specified colour channel of an image at a given location.
	 * There is no bounds checking or other safety checks, but this should be very fast.
	 * This only works for 8-bit images (type IPL_DEPTH_8U).
	 *
	 * \param img The image to get pixel values from
	 * \param x   The x-coordinate of the pixel to get
	 * \param y   The y-coordinate of the pixel to get
	 * \param c   The (0-based) index of the colour channel to get
	 * \param value  The new colour value
	 */
	static inline void mySet2D(IplImage *img, int x, int y, unsigned char value) HOT HARD_INLINE;
	static inline void mySet2D(IplImage *img, int x, int y, unsigned char value) {
		img->imageData[y*img->widthStep + x*img->nChannels] = value;
	}

	static inline unsigned char *myGetRowStart(IplImage *img, int y) {
		return (unsigned char *)(img->imageData + y*img->widthStep);
	}

}

#endif

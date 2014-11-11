/* \file image.h \brief Wraps Intel IplImage to provide easy access to pixel values
 *
*/

#pragma once

#ifndef _grc_image_
#define _grc_image_

#include "util/opencv.h"

#include <iostream>
#include <string>
#include <limits>

pragma_warning (push)
pragma_warning (disable: 4100)

using namespace std;

namespace grc {

enum ImageDepth { GREY = 1, RGB = 3, RGBA = 4 };

template <ImageDepth t_depth = RGB, class t_DataType = unsigned char> class Image {

    bool imageIsCopy_;
public:

    Image() : imageIsCopy_(false) {
		m_img = NULL;
	}

	Image(const Image &I) : imageIsCopy_(false) {
		if (I.m_img)
			m_img = cvCloneImage(I.m_img);
		else
			m_img = NULL;
	}

	Image(unsigned int width, unsigned int height) : imageIsCopy_(false) {
		m_img = cvCreateImage(cvSize(width, height), getIplDepth(), int(t_depth));
		m_img->origin = 1;
	}

	Image(const IplImage *img, bool imageIsCopy = false) : imageIsCopy_(imageIsCopy) {
        if(imageIsCopy_)
        {
            m_img = const_cast<IplImage *>(img);
        }
        else
        {
		    if (img->nChannels != t_depth) {
			    cerr << "Warning, different number of channels in Image(const IplImage *img)" << endl;
		    }
		    if (img->depth != getIplDepth()) {
			    cerr << "Warning, different pixel types in Image(const IplImage *img)" << endl;
		    }
		    m_img = cvCloneImage(img);
		    cvFlip(m_img);
        }

	}

	~Image() {
		if (m_img && !imageIsCopy_) cvReleaseImage(&m_img);
	}

	Image &operator=(const Image &I) {
		if (this != &I) {
			// Not self-assignment
			if (I.m_img) {
				if (m_img && !imageIsCopy_) cvReleaseImage(&m_img);
				m_img = cvCloneImage(I.m_img);
			} else {
				m_img = NULL;
			}
		}
		return *this;
	}

	//! Red, Green, Blue, Alpha. NB: Use Blue channel if its a greyscale image or will be OOB!
	t_DataType &r(unsigned int x, unsigned int y) { return index(x,y,2); }
	t_DataType &g(unsigned int x, unsigned int y) { return index(x,y,1); }
	t_DataType &b(unsigned int x, unsigned int y) { return index(x,y,0); }
	t_DataType &a(unsigned int x, unsigned int y) { return index(x,y,3); }

	// Const correct versions
	const t_DataType &r(unsigned int x, unsigned int y) const { return index(x,y,2); }
	const t_DataType &g(unsigned int x, unsigned int y) const { return index(x,y,1); }
	const t_DataType &b(unsigned int x, unsigned int y) const { return index(x,y,0); }
	const t_DataType &a(unsigned int x, unsigned int y) const { return index(x,y,2); }

	// Dimensions
	unsigned int getWidth() const { return m_img->width; }
	unsigned int getHeight() const {return m_img->height; }
	ImageDepth getDepth() const {return t_depth; }

	operator IplImage*() { return m_img;}
	operator const IplImage*() const { return m_img;}

private:

	IplImage *m_img;

	int getIplDepth() const {
		t_DataType dummy = 0;
		return getIplDepth(dummy);
	}

	template <class T> int getIplDepth(T dummy) const {
		throw "This will throw an exception if someone tries to make an image with the wrong type.";
		//The explicit overrides that follow are used for the correct types, and avoid this 'code'.
	}
	int getIplDepth(char dummy) const { return IPL_DEPTH_8S; }
	int getIplDepth(unsigned char dummy) const { return IPL_DEPTH_8U; }
	int getIplDepth(short dummy) const { return IPL_DEPTH_16S; }
	int getIplDepth(unsigned short dummy) const { return IPL_DEPTH_16U; }
	int getIplDepth(int dummy) const { return IPL_DEPTH_32S; }
	int getIplDepth(float dummy) const { return IPL_DEPTH_32F; }
	int getIplDepth(double dummy) const { return IPL_DEPTH_64F; }

	t_DataType &index(int x, int y, int c) {
		return ( (t_DataType *)(m_img->imageData + m_img->widthStep*y) )[x*t_depth + c];
	}
	const t_DataType &index(int x, int y, int c) const {
		return ( (t_DataType *)(m_img->imageData + m_img->widthStep*y) )[x*t_depth + c];
	}
}; // class image


} // namespace grc

pragma_warning (pop)

#endif

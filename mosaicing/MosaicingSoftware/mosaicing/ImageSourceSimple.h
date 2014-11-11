#pragma once

namespace grc
{

//! Parent class for image source that avoids the need for circular header file includes (and only exposes a minimal interface to the renderer, etc.).
class ImageSourceSimple
{
public:
    //! Returns a pointer to a frame. Don't release this image.
	virtual const IplImage * getImage(size_t id) = 0;
};

}


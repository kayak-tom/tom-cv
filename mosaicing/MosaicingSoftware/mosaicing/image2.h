/* \file image2.h \brief Useful typedefs.
 *
 */
#pragma once

#include "image.h"
namespace grc 
{
	typedef Image <RGB, unsigned char> ImageRGB; //! RGB image type
	typedef Image <GREY, unsigned char> ImageGREY; //! Greyscale image type
}

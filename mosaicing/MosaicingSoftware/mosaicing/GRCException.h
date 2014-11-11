/*! \file grcException.h
 * \brief Exception handling
 *
 *  Defines exceptions and GRCCATCH/FINALLY exception handling macros.
 */
#pragma once
#ifndef _GRC_EXCEPTION
#define _GRC_EXCEPTION

#include "util/exception.h"

namespace grc
{
typedef CException GRCException;
}

/*
#include <exception>
#include <string>

namespace grc
{

//! Generic exception class.
/ *! Generally thrown when a consistency check fails, indicating an error in the program.
 * /
class GRCException
#ifndef GCC
	: public std::exception
#endif
{
protected:
    std::string * pstrErrorMsg_;
public:
    //! Constructor takes an error message that may be displayed by handler.
    GRCException(const char * szErrorMsg ///< Pass in an error message
        )
    {
        pstrErrorMsg_ = new std::string(szErrorMsg);
    };

    ~GRCException(){ delete pstrErrorMsg_; };

    const std::string * getErrorMessage() { return pstrErrorMsg_; };
};

//! Can be #defined to nothing for faster execution (generally unused) 
#define if(IS_DEBUG) CHECK(x, s) if(x) throw new GRCException(s);

}*/

#endif // _GRC_EXCEPTION

/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * imageSourceFromDir.h
 *
 *  Created on: 5/06/2009
 *      Author: tom
 */

#ifndef IMAGESOURCEFROMDIR_H_
#define IMAGESOURCEFROMDIR_H_

#include "imageSource.h"
#include <vector>
#include <string>
#include <boost/thread/mutex.hpp>

class CImageSourceFromDir: public CImageSource
{
    std::vector<std::string> aAllFileNames;
    enum eImTypes { eImNone, eImJPG, eImPNG, eImTIFF, eImBMP, eImPGM };
    eImTypes imType;
    eImTypes getImType(const std::string & str) const;
    boost::mutex imLoadMutex;
    const int nLoadSubset;
    int nLastId, nLastImageId;
    std::vector<int> aJumps;

    void init(const char * szFolder, CImParams & IM_PARAMS_in);
public:
    CImageSourceFromDir(CImParams & IM_PARAMS_in);
    //only use this one when multiple image loaders are needed
    CImageSourceFromDir(const char *, CImParams & IM_PARAMS_in);
    virtual ~CImageSourceFromDir();

    virtual bool loadImage_int(int nId, IplImage * pIm) HOT;
    virtual bool loadNextImage_int(int & nId, IplImage * pIm) ;
    virtual void setGlobalDims(int & x, int & y, int & channels) ;
};

#endif /* IMAGESOURCEFROMDIR_H_ */

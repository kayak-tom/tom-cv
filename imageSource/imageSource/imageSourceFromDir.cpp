/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * **** if not building with GCC4.4 you may need this patch https://svn.boost.org/trac/boost/ticket/3041 *******
 *
 *
 * imageSourceFromDir.cpp
 *
 *  Created on: 5/06/2009
 *      Author: tom
 */

#include "imageSourceFromDir.h"
#include <boost/filesystem.hpp>
#include "util/opencv.h"
#include "util/opencv_highgui.h"
#ifdef __GNUC__
#include <boost/gil/extension/io/jpeg_io.hpp>
#include <boost/gil/extension/io/png_io.hpp>
#ifdef HAVE_TIFF
#include <boost/gil/extension/io/tiff_io.hpp>
#endif
#endif
#include <boost/gil/gil_all.hpp>
#include <boost/thread/mutex.hpp>
#include <fstream>
#include "image/convert_OpenCV.h"

using namespace std;
using namespace boost::filesystem;

void greyImCopy(const boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::gray8_pixel_t*> > > & imView, IplImage *pIm) HOT;
void greyImCopy(const boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::gray8_pixel_t*> > > & imView, IplImage *pIm)
{
    CIplImIt<uchar> imOutIt(pIm);
    for(int y = 0;y < imView.height();++y){
        boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::gray8_pixel_t*> > >::x_iterator src_it = imView.row_begin(y);
        for(int x = 0;x < imView.width();++x, src_it++, imOutIt++)
            *imOutIt = *src_it;
    }
}
void rgbImCopy(const boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::rgb8_pixel_t*> > > & imView, IplImage *pIm) HOT;
void rgbImCopy(const boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::rgb8_pixel_t*> > > & imView, IplImage *pIm)
{
    CIplImIt<uchar> imOutIt(pIm);
    for(int y = 0;y < imView.height();++y)
    {
        //* **** if not building with GCC4.4 you may need this patch https://svn.boost.org/trac/boost/ticket/3041 *******
        boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::rgb8_pixel_t*> > >::x_iterator src_it = imView.row_begin(y);
        for(int x = 0;x < imView.width(); ++x, imOutIt++, src_it++)
        {
            const boost::gil::rgb8_pixel_t & col = *src_it;
            /*CIplPx<unsigned char>::setRed(pIm, x, y, boost::gil::get_color(col, boost::gil::red_t()));
            CIplPx<unsigned char>::setGreen(pIm, x, y, boost::gil::get_color(col, boost::gil::green_t()));
            CIplPx<unsigned char>::setBlue(pIm, x, y, boost::gil::get_color(col, boost::gil::blue_t()));*/
            imOutIt.setBGR(boost::gil::get_color(col, boost::gil::blue_t()), boost::gil::get_color(col, boost::gil::green_t()), boost::gil::get_color(col, boost::gil::red_t()));
        }
    }
}

void greyImCopyOld(const boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::gray8_pixel_t*> > > & imView, IplImage *pIm)
{
    for(int y = 0;y < imView.height();++y){
        boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::gray8_pixel_t*> > >::x_iterator src_it = imView.row_begin(y);
        for(int x = 0;x < imView.width();++x)
            CIplPx<unsigned char>::setGrey(pIm, x, y, src_it[x]);
    }
}

void rgbImCopyOld(const boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::rgb8_pixel_t*> > > & imView, IplImage *pIm)
{
    for(int y = 0;y < imView.height();++y)
    {
        boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::rgb8_pixel_t*> > >::x_iterator src_it = imView.row_begin(y);
        for(int x = 0;x < imView.width();++x)
        {
            CIplPx<unsigned char>::setRed(pIm, x, y, boost::gil::get_color(src_it[x], boost::gil::red_t()));
            CIplPx<unsigned char>::setGreen(pIm, x, y, boost::gil::get_color(src_it[x], boost::gil::green_t()));
            CIplPx<unsigned char>::setBlue(pIm, x, y, boost::gil::get_color(src_it[x], boost::gil::blue_t()));
        }
    }
}

CImageSourceFromDir::eImTypes CImageSourceFromDir::getImType(const string & str_in) const
{
    string str = str_in;
    for(int i=0; i<(int)str.size();++i)
       str[i] = tolower(str[i]);

    size_t extPos = str.length()-4;

    if(extPos == str.find(".jpg")) return eImJPG;
    if(extPos == str.find(".png")) return eImPNG;
    if(extPos == str.find(".tif")) return eImTIFF;
    if(extPos == str.find(".bmp")) return eImBMP; //todo
    if(extPos == str.find(".pgm")) return eImPGM; //todo

    /*if(string::npos != str.find(".jpg")) return eImJPG;
    if(string::npos != str.find(".png")) return eImPNG;
    if(string::npos != str.find(".tif")) return eImTIFF;
    if(string::npos != str.find(".bmp")) return eImBMP; //todo*/
    return eImNone;
}

CImageSourceFromDir::CImageSourceFromDir(const char * szFolder, CImParams & IM_PARAMS_in) : CImageSource(IM_PARAMS_in), nLoadSubset(IM_PARAMS.ImageDir.LOAD_SUBSET)
{
    if(IM_PARAMS.ImageDir.IMAGE_DIR.isInit())
    {
        cout << "WARNING: Overriding image source dir given in config file (" << (const string &)IM_PARAMS.ImageDir.IMAGE_DIR << ") with " << szFolder << endl;
    }
    init(szFolder, IM_PARAMS_in);
}

CImageSourceFromDir::CImageSourceFromDir(CImParams & IM_PARAMS_in) : CImageSource(IM_PARAMS_in), nLoadSubset(IM_PARAMS.ImageDir.LOAD_SUBSET)
{
    cout << IM_PARAMS.ImageDir.IMAGE_DIR;
    CHECK(!IM_PARAMS.ImageDir.IMAGE_DIR.isInit() || ((const string &)IM_PARAMS.ImageDir.IMAGE_DIR).length()==0, "Image source folder is required, specify Im.ImageDir.IMAGE_DIR=\"/path/to/image/dir\" in config file")
    init(IM_PARAMS.ImageDir.IMAGE_DIR.asSz(), IM_PARAMS_in);
}

void CImageSourceFromDir::init(const char * szFolder, CImParams & IM_PARAMS_in)
{
    nLastId = -1; nLastImageId = -1;

    CHECK(!szFolder, "Image source folder is null; should be specified in config file (Im.ImageDir.IMAGE_DIR=\"/path/to/image/dir\")")
    if(!boost::filesystem::exists(szFolder))
    {
        cout << "Folder \"" << szFolder << "\" does not exist\n";
        THROW("Image source folder does not exist")
    }
    else if(!boost::filesystem::is_directory(szFolder))
    {
        cout << "Folder \"" << szFolder << "\" is not a directory\n";
        THROW("Image source folder is not a directory")
    }

    imType = eImNone;
    eImTypes possType = eImNone;
    directory_iterator end_itr;
    for ( directory_iterator itr( szFolder ); itr != end_itr; ++itr )
    {
        if ( is_regular(itr->status()) )
        {
            const path p = itr->path();
//            string strImFilename = p.relative_path().file_string();
            string strImFilename = p.string();
            eImTypes thisType = getImType(strImFilename);
            if(eImJPG == thisType)
            {
                imType = eImJPG;
                //cout << "Found " << strImFilename << endl;
                aAllFileNames.push_back(strImFilename);
            }
            else if (thisType != eImNone)
                possType = thisType;
            else
                if(strImFilename.find("JUMPS") == strImFilename.length() - 5)
                {
                    ifstream jumps(strImFilename.c_str());
                    int n=-1, nLast = -1;
                    while(!jumps.eof())
                    {
                        jumps >> n;

                        CHECK(n<=0, "Bad val for jump val");

                        if(n != nLast)
                        {
                            nLast = n;
                            cout << "Will jump at id " << n << endl;
                            aJumps.push_back(n);
                        }
                    }
                }
        }
    }
    if(aAllFileNames.size()==0 && possType != eImNone) //try pngs or whatever:
        for ( directory_iterator itr( szFolder ); itr != end_itr; ++itr )
        {
            if ( is_regular(itr->status()) )
            {
                const path p = itr->path();
                string strImFilename = p.string();
                if(possType == getImType(strImFilename))
                {
                    imType = getImType(strImFilename);
                    aAllFileNames.push_back(strImFilename);
                }
            }
        }

    //Now sort by name so that ids represent times
    std::sort(aAllFileNames.begin(), aAllFileNames.end());

    CHECK(aAllFileNames.size()==0, "No images found in image source directory");

    if(!IM_PARAMS.IM_WIDTH.isInit())
    {
        cout << "Dynamically initialising image height and width and colour channels...\n";
        setImageParams(IM_PARAMS_in);

        IM_PARAMS_in.initCalibration(szFolder);
    }
    else
    {
        cout << "Image height and width and colour channels initialised already by a different image loader\n";
        setImageParams(IM_PARAMS_in);
    }
}

CImageSourceFromDir::~CImageSourceFromDir()
{
}

void CImageSourceFromDir::setGlobalDims(int & x, int & y, int & channels)
{
    boost::gil::point2<std::ptrdiff_t> dims;

#ifdef __GNUC__
    if(imType == eImJPG)
    {
        boost::gil::detail::jpeg_reader m(aAllFileNames[0].c_str());
        dims = m.get_dimensions();    //    dims = boost::gil::jpeg_read_dimensions(aAllFileNames[0].c_str());
        try
        {
            boost::gil::rgb8_image_t img( dims.x, dims.y );
            m.apply(boost::gil::view(img));
            channels = 3;
        }
        catch(std::ios_base::failure)
        {
            channels = 1;
        }
    }
    else if(imType == eImPNG)
    {
        boost::gil::detail::png_reader m(aAllFileNames[0].c_str());
        dims = m.get_dimensions();    //    dims = boost::gil::jpeg_read_dimensions(aAllFileNames[0].c_str());
        try
        {
            boost::gil::rgb8_image_t img( dims.x, dims.y );
            m.apply(boost::gil::view(img));
            channels = 3;
        }
        catch(std::ios_base::failure)
        {
            channels = 1;
        }
    }
    else if(imType == eImTIFF)
    {
#ifdef HAVE_TIFF
        boost::gil::detail::tiff_reader m(aAllFileNames[0].c_str());
        dims = m.get_dimensions();    //    dims = boost::gil::jpeg_read_dimensions(aAllFileNames[0].c_str());
        try
        {
            boost::gil::rgb8_image_t img( dims.x, dims.y );
            m.apply(boost::gil::view(img));
        }
        catch(std::ios_base::failure)
        {
            nDepth = 1;
        }
#endif
    }
    if(dims.x>0&&dims.y>0)
    {
        x = dims.x;
        y = dims.y;
    }
#endif

#ifdef __GNUC__
    if(imType == eImBMP || imType == eImPGM) //or pgm
#endif
    {
        IplImage * pIm = cvLoadImage(aAllFileNames[0].c_str());
        if(!pIm)
        {
            cout << "ERROR loading " << aAllFileNames[0] << " to determine file dims\n";
            THROW( "Error getting dimensions");
        } 
#ifdef __GNUC__
        if(imType == eImBMP)//only for bitmaps
        {
            ifstream fileGetDepth(aAllFileNames[0].c_str(), ios::binary);
            fileGetDepth.seekg (28, ios::beg);
            int nBitDepth = 0;
            fileGetDepth.read ((char*)&nBitDepth, sizeof(int));
            fileGetDepth.close();
            channels = nBitDepth==24 ? 3 : 1; // doesn't work pIm->nChannels;
        } else
            channels = pIm->nChannels; //Assume greyscale (??)
#else
            cout << "COLOUR IMAGES assumed under Windows (I think? todo: fix)\n";
            channels = pIm->nChannels;
#endif


        x = pIm->width;
        y = pIm->height;
        cvReleaseImage(&pIm);
    }

    if(channels == 1)
        cout << "GREYSCALE images found\n";
    else
        cout << "RGB images found\n";
}

//Also set id/
bool CImageSourceFromDir::loadNextImage_int(int & nId, IplImage * pIm)
{
    if(nLastId < 0)
    {
        nLastId = 0;
        nLastImageId = 0;
    }
    else
    {
        nLastId++;
        nLastImageId++;
    }
    for(std::vector<int>::const_iterator pJump = aJumps.begin(); pJump != aJumps.end(); pJump++)
        if(nLastImageId == *pJump)
            nLastId+=1000;

    nId = nLastId;

    cout << "Returning image " << nLastImageId << " with id " << nId << endl;
    return loadImage_int(nLastImageId, pIm);
}

bool CImageSourceFromDir::loadImage_int(int nId, IplImage * pIm)
{
    nId *= nLoadSubset;
    if(nId >= (int)aAllFileNames.size()) return false;

    const char * szFileName = aAllFileNames[nId].c_str();
    CHECK(SOURCE_WIDTH != (int)pIm->width || SOURCE_HEIGHT != (int)pIm->height, "CImageSourceFromDir::loadImage: Supplied image dims do not match SOURCE_WIDTH, SOURCE_HEIGHT");

    boost::unique_lock<boost::mutex> write_lock(imLoadMutex);

#ifdef __GNUC__
    if(imType == eImBMP || imType == eImPGM)
#endif
    {
        IplImage * pImIn = cvLoadImage(szFileName);
        cvCopy(pImIn, pIm);
        cvReleaseImage(&pImIn);
    }
#ifdef __GNUC__
    else
    {
        try
        {
            if(SOURCE_CHANNELS == 3)
            {
                static boost::gil::rgb8_image_t img( SOURCE_WIDTH, SOURCE_HEIGHT );
                const boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::rgb8_pixel_t*> > > & imView = boost::gil::view( img );
                if(    imType == eImJPG)
                {
                    boost::gil::jpeg_read_view( szFileName,  imView  );
                }
                else if(imType == eImPNG)
                {
                    boost::gil::png_read_view( szFileName,  imView  );
                }
                else if(imType == eImTIFF)
                {
#ifdef HAVE_TIFF
                    boost::gil::tiff_read_view( szFileName,  imView  );
#endif
                }
                else
                    THROW( "Unrecognised im type")

                rgbImCopy(imView, pIm);
            }
            else if(SOURCE_CHANNELS == 1)
            {
                static boost::gil::gray8_image_t img( SOURCE_WIDTH, SOURCE_HEIGHT );
                const boost::gil::image_view<boost::gil::memory_based_2d_locator<boost::gil::memory_based_step_iterator<boost::gil::gray8_pixel_t*> > > & imView = boost::gil::view( img );

                if(    imType == eImJPG)
                {
                    boost::gil::jpeg_read_view( szFileName,  imView  );
                }
                else if(imType == eImPNG)
                {
                    boost::gil::png_read_view( szFileName,  imView  );
                }
                else if(imType == eImTIFF)
                {
#ifdef HAVE_TIFF
                    boost::gil::tiff_read_view( szFileName,  imView  );
#endif
                }
                else
                    THROW( "Unrecognised im type")

                greyImCopy(imView, pIm);
            }
        }
        catch(std::ios_base::failure)
        {
            THROW( "CImageSourceFromDir::loadImage: boost::gil::jpeg_read_view failed--wrong size/depth. Has the image size changed?");
        }
        catch(...)
        {
            THROW( "CImageSourceFromDir::loadImage: boost::gil::jpeg_read_view failed--unknown exception. Define HAVE_TIFF and link to libtiff to support TIFF images");
        }
    }
#endif
    return true;
}

/* \file This file defines enums used in both C# and C++. The GUI sets them from 
 * its combo box selected-indices (so the order should not be changed).
 * 
 */

namespace grc
{
    namespace Enums
    {
        enum eVideoSource { eAttachedCam, eVideoFile, eImageDirectory }; //!< Video or image source
        enum eDisplayMode { eDisplayOnly, eSaveImages }; //!< Controls whether to save mosaics as well as displaying them (save video not implemented)
        enum eTransformType { ePerspective, eAffine, eSimilarity }; //!< Transform type to use. Order MUST MATCH combo box in gui
        enum eFeatureDetector { eShiTomasi, eShiTomasiSubpix, eDoGBlob };//!< Salient feature detector to use. Order MUST MATCH combo box in gui
        enum eDescriptorInvarianceType { eSimplePatch, eOriented, eNormalised, eOrientedNormalised };//!< Patch descriptor to use. Order MUST MATCH combo box in gui
        enum eDescriptorSize { e9x9, e11x11, e13x13, e15x15 };//!< Size of patch descriptor to use. Order MUST MATCH combo box in gui
        enum eRenderer { eBasicRenderer, eDijkstraRenderer, eFeatheredRenderer, eMultiScaleRenderer };//!< renderer to use. Order MUST MATCH combo box in gui
        enum eRendererEvalFunction { eNoEvaluation, eSSDEvaluation };//!< Mosaic evaluation function to use. Order MUST MATCH combo box in gui
        enum eWarpMethod { eWarpNN, eWarpBilinear };//!< Interpolation method to use when warping images into mosaic. Nearest-Neighbour is faster, Bilinear is slightly better. Order MUST MATCH combo box in gui
        enum ePositionMethod { ePosCentre, ePosLatest }; //!< Method used to position frames within mosaic. ePosLatest will fit the latest frame and as much as possible of the rest of the mosaic.
        enum eFrameChoiceMethod { eChooseSequential, eChooseSequentialSkip, eChooseBoWRecurse }; //!< Choose whether to render frames sequentially, or to select some to give a larger mosaic.
    }
}

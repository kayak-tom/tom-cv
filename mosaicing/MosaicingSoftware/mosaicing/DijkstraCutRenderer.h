#ifndef _GRC_DIJKSTRA_CUT_RENDERER_
#define _GRC_DIJKSTRA_CUT_RENDERER_

#include "Renderer.h"
#include "util/set2.h"

namespace grc {

	/*!
	 * \brief Reads the latest TransformSet from the TransformEngine and renders an image from it using an optimal seam method.
	 *
	 * Renderers are responsible for creating a mosaic from a set of images and Transforms that align 
	 * them to a common co-ordinate frame. The DijkstraCutRenderer creates a mosaic by using a graph search
	 * to find seams along which to cut the images in order to create a visually pleasing result.
	 */	
	class DijkstraCutRenderer : public Renderer {

	public:
	    DijkstraCutRenderer(ImageSource & imageSource, TransformEngine & transformEngine, CvSize mosaicSize, unsigned int scale, EvaluationFunction *evaluationFunction);
	    ~DijkstraCutRenderer();
	protected:
	    IplImage *renderImagesToMosaic(TransformSet *TS);
	private:
	    template<unsigned int SCALE, unsigned int CHANNELS>
	    void findAndApplyCut_template(const IplImage *frame, const Transform *transform);

	    template<unsigned int CHANNELS>
	    void findAndApplyCut2(const IplImage *frame, const Transform *transform);

	    void findAndApplyCut(const IplImage *frame, const Transform *transform);
	    void initialiseStorage(TransformSet *TS);

	    template<unsigned int SCALE, unsigned int CHANNELS>
	    void resizeAndThresh();

	    bool initialised_;
	    unsigned int scale_;
	    CvSize frameSize_;
	    CvSize scaledFrameSize_;
	    CvSize scaledMosaicSize_;
	    IplImage *mosaicMask_;
	    IplImage *warpedMask_;
	    //IplImage *mosaicEdge_;
	    //IplImage *warpedEdge_;
	    //IplImage *crossingMask_;
	    IplImage *totalMask_;
	    IplImage *largeMask_;
	    //IplImage *checkMask_;
	    //IplImage *checkFrame_;
	    IplImage *largeMosaic_;
	    IplImage *largeWarped_;
	    CvMat *scaleTransform_;
	    CvMat *corners_;
	    CvMat *warpedCorners_;
	    double **distance_;
	    double **cost_;
	    bool **expanded_;
	    signed char **dx_;
	    signed char **dy_;
	    class PriorityQueue
	    {
	    public:
	        PriorityQueue();
	        ~PriorityQueue();
	        void push(int x, int y, double oldDistance, double newDistance);
	        bool isEmpty() const;
	        void pop(int *x, int *y);
	        void clear();
	    private:
	        //Probably should create a new queue each time
	        typedef multimap<double, pair<int,int> ,less<double>, CIndividualPoolAllocator<std::pair<const int, int>, 48 > > TQueue;
			TQueue Q_;
	    };

	};

}

#endif

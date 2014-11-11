/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * plot.h
 *
 *  Created on: 28/05/2009
 *      Author: tom
 */

#ifndef PLOT_H_
#define PLOT_H_

#include<map>
#include "util/opencv.h"

class CPlot
{
public:
	static void plot(IplImage * image, const std::map<int, double> & aData, bool bPlot0);
	static void setWhite(IplImage * image);
};

#endif /* PLOT_H_ */

/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * bowSpeedoParams.h
 *
 *  Created on: 30/10/2009
 *      Author: tom
 */

#ifndef BOWSPEEDOPARAMS_H_
#define BOWSPEEDOPARAMS_H_

#include "params/param.h"

PARAMCLASS(BOWSpeedo)
	PARAM(NUM_OBJECTS, 0, 1000000, 250, "Detect this many objects. 0 turns off speedo.")
	PARAM(MAX_PROP_COOCCURANCES, 0.000001, 1, 0.01, "Use at most this proportion of objects (only relevent for first few frames)") //Top 1%
	PARAMB(DRAW_OBJECTS_AFTER_RETRAIN, false, "Output all frames with best objects marked. Subject to random-access image loader (=from directory).")
	{}

	CNumParam<int> NUM_OBJECTS;
	CNumParam<double> MAX_PROP_COOCCURANCES;
	CNumParam<bool> DRAW_OBJECTS_AFTER_RETRAIN;
};


#endif /* BOWSPEEDOPARAMS_H_ */

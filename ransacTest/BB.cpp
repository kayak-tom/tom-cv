#include "util/random.h"
#include "util/dynArray.h"
#include <fstream>
#include <algorithm>

class CBBSimulator
{
	CDynArray<double> aValues;
	const int nSteps;
public:
	CBBSimulator(int nSteps = 1000000) : aValues(nSteps+1), nSteps(nSteps) {}

	void simulate()
	{
		aValues[0]=0;
		double dSD = sqrt(1.0/nSteps);
		for(int i=1; i<=nSteps; i++)
			aValues[i]=aValues[i-1]+CRandom::Normal(0, dSD);

		double dEnd = aValues[nSteps];
		double dStep = dEnd/nSteps;
		//cout << dEnd << endl;

		for(int i=1; i<=nSteps; i++)
			aValues[i] -= i*dStep;
	}

	//Return the val of k where the BB will lie entirely above k*t(1-t)
	double getk()
	{
		double k = 0, t=0;
		double dtStep = 1.0/nSteps;
		for(int i=0; i<=nSteps; i++, t += dtStep)
		{
			if(aValues[i] < 0)
			{
				double dkNew = aValues[i]/sqrt((t*(1-t)));
				if(dkNew < k)
					k=dkNew;
			}
		}
		
		return k;
	}
};

void findK()
{
	const int iters=100000, nSteps=100000;
	CDynArray<double> ak(iters);
	CBBSimulator bb(nSteps);
	for(int i=0;i<iters;i++)
	{
		bb.simulate();
		ak[i] = bb.getk();
	}
	std::ofstream kfile("kfile.tsv");
	
	const bool bOutputAll = false;
	if(bOutputAll)
	{
		std::sort(ak.begin(), ak.end());
		ak.pp("\n", kfile);
	}
	else
	{
		double dMean, dVar;
		ak.getMeanVar(dMean, dVar);
		kfile << dMean << endl;
		kfile << dVar << endl;
	}
	
	kfile.close();



}
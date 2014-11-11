/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#define templateStat template<class TNumber, bool QUAL_SCORE, eVarMeasure eVarianceMethod>
#define CCBoWCorrespondences CBoWCorrespondences<TNumber, QUAL_SCORE, eVarianceMethod>

enum eVarMeasure {eStatVar, eStatRange};
//Point correspondences between 2 images
templateStat
class CBoWCorrespondences : public std::vector<CCorrespondence>
{
//	std::vector<CBoWCorrGroup *> vCorrGroups;
//	typedef std::vector<CBoWCorrGroup *> corrGrpIter;

    static const int INT_SCALE = 1;//256;
    static const int NO_STAT = -1;

    TNumber nStat;
    static TNumber getVar(const TNumber * v, int N);
    static TNumber getWeightedVar(const TNumber * v, const TNumber * weight, int N);
    static inline void scaleProjection(TNumber * axProj, int N);
    static void scaleProjection(TNumber * axProj, TNumber SCALE, int N);
    static TNumber getProjectionStatInt(const TNumber * x_proj_x, const TNumber * xp_proj_x, const TNumber * xp_proj_y, const TNumber * matchQual, int N);
    static void normaliseProjection(TNumber * axProj, TNumber SCALE, int N);

    static void GetSums(const TNumber *p1, const TNumber *p2, TNumber *pSum, TNumber *pDif, int N);

    static TNumber scaleRatio(TNumber a, TNumber b);
    //TNumber * operator+(const TNumber *p1, const TNumber *p2)
    //TNumber * operator-(const TNumber *p1, const TNumber *p2)
public:
    static const int MIN_CORRESPONDENCES = 3;
    TNumber CompStat() ;

    CBoWCorrespondences() : nStat(NO_STAT) { reserve(450); };

    inline TNumber Stat() const
    {
        if(IS_DEBUG) CHECK(nStat == NO_STAT, "Stat not setup");
        return nStat;
    };

    ~CBoWCorrespondences()
    {
//        iterator ppEnd = end();
//        for(iterator ppCorrespondence = begin(); ppCorrespondence<ppEnd; ppCorrespondence++)
//            delete *ppCorrespondence;
//
//        clear();
    };

    /*void push_back(CBoWCorrGroup * pCorrGroup)
    {
    	vCorrGroups.push_back(pCorrGroup);
    };

    int numSelectable() const
    {
    	int nSelectable = 0;
    	corrGrpIter ppCGend = vCorrGroups.end();
    	for(corrGrpIter ppCG = vCorrGroups.begin(); ppCG < ppCGend; ppCG++)
    		nSelectable += (*ppCG)->maxPossCorrespondences();
    };*/


    typedef TNumber TNumberType;
};

//Normalise distn to mean 0, S.D. SCALE
templateStat
void CCBoWCorrespondences::normaliseProjection(TNumber * axProj, TNumber SCALE, int N)
{

    //Find mean
    TNumber * pxProj = axProj;
    TNumber Mean = 0;
    for(int i=N; i>0; i--)
    {
        Mean += *pxProj;
        pxProj++;
    }
    Mean /= N;

    //Find Var
    TNumber Var = 0;
    pxProj = axProj;
    for(int i=N; i>0; i--)
    {
        TNumber dif = *pxProj - Mean;
        Var += dif*dif;
        pxProj++;
    }
    Var /= N;
    TNumber SD = sqrt(Var);

    pxProj = axProj;
    for(int i=N; i>0; i--)
    {
        TNumber x_elem = *pxProj;
        TNumber newElem = (SCALE*(x_elem-Mean)) / SD;
        *pxProj = newElem;
        pxProj++;
    }
};

//Scale to -k..k
templateStat
void CCBoWCorrespondences::scaleProjection(TNumber * axProj, TNumber SCALE, int N)
{
    TNumber * pxProj = axProj;

    TNumber Min = MAX_INT, Max=-Min;

    //Find min+max
    for(int i=N; i>0; i--)
    {
        TNumber x_elem = *pxProj;
        if(Min > x_elem)
            Min = x_elem;
        if(Max < x_elem)
            Max = x_elem;

        pxProj++;
    }

    TNumber HalfRange = (Max-Min)/2;
    TNumber Sub = Min + HalfRange ;

    if(HalfRange==0)
    {
        HalfRange = 1;
    }

    pxProj = axProj;
    for(int i=N; i>0; i--)
    {
        TNumber x_elem = *pxProj;
        TNumber newElem = (SCALE*(x_elem-Sub))/HalfRange;
        *pxProj = newElem;
        pxProj++;
    }
};

/* calc p1+p2 and p1-p2 (vectors), put into pSum and pDif */
templateStat
void CCBoWCorrespondences::GetSums(const TNumber *p1, const TNumber *p2, TNumber *pSum, TNumber *pDif, int N)
{
    for(int i=N; i>0; i--)
    {
        TNumber n1=*p1;
        TNumber n2=*p2;
        *pSum = n1 + n2;
        *pDif = n1 - n2;
        p1++;
        p2++;
        pSum++;
        pDif++;
    }
}
extern int statQS;
templateStat
TNumber CCBoWCorrespondences::getProjectionStatInt(const TNumber * x_proj_x, const TNumber * xp_proj_x, const TNumber * xp_proj_y, const TNumber * matchQual, int N)
{
    TNumber * aSum = new TNumber[2*N];
    TNumber * aDif = aSum+N;
    TNumber xx_var, x_x_var, xy_var, x_y_var;

    GetSums(x_proj_x, xp_proj_x, aSum, aDif, N);
    if(statQS)//QUAL_SCORE
    {
        xx_var=getWeightedVar(aSum,matchQual, N);
        x_x_var=getWeightedVar(aDif,matchQual, N);
    } else {
        xx_var=getVar(aSum, N);
        x_x_var=getVar(aDif, N);
    }
    GetSums(x_proj_x, xp_proj_y, aSum, aDif, N);
    if(QUAL_SCORE)
    {
        xy_var=getWeightedVar(aSum,matchQual, N);
        x_y_var=getWeightedVar(aDif,matchQual, N);
    } else {
        xy_var=getVar(aSum, N);
        x_y_var=getVar(aDif, N);
    }

    delete [] aSum; aSum=aDif=0;

    TNumber min1 = min(xx_var,x_x_var);
    TNumber min2 = min(min1,xy_var);
    TNumber min_var=min(min2,x_y_var);
    //this is recursive and horrible when macros expand but more obvious double min_var=min(min(min(xx_var,x_x_var),xy_var),x_y_var);

    if (min_var==xx_var)
        return scaleRatio(xx_var,x_x_var);
    if (min_var==x_x_var)
        return scaleRatio(x_x_var,xx_var);
    if (min_var==xy_var)
        return scaleRatio(xy_var,x_y_var);

    return scaleRatio(x_y_var, xy_var);
};

templateStat
inline TNumber CCBoWCorrespondences::scaleRatio(TNumber a, TNumber b)
{
    return (INT_SCALE*a)/b;//a/b;
};

//template<bool QUAL_SCORE>
//inline int CBoWCorrespondences<int, QUAL_SCORE>::scaleRatio(int a, int b)
//{
//    return (INT_SCALE*a)/b;
//};

//Get some measure of variance
templateStat
TNumber CCBoWCorrespondences::getWeightedVar(const TNumber * v, const TNumber * weight, int N)
{
    TNumber SS = 0;
    const TNumber * pv = v;

    for(int i=N; i>0; i--)
    {
        TNumber t=*pv;
        TNumber w=*weight;
        SS += w*t*t;
        pv++;
        weight++;
    }

    return SS/N;
};

//Get some measure of variance
templateStat
TNumber CCBoWCorrespondences::getVar(const TNumber * v, int N)
{
    TNumber SS = 0;
    const TNumber * pv = v;

    for(int i=N; i>0; i--)
    {
        TNumber t=*pv;
        SS += t*t;
        pv++;
    }

    return SS/N;
};
extern int b45Deg;
templateStat
TNumber CCBoWCorrespondences::CompStat()
{
    if(IS_DEBUG) CHECK(nStat != NO_STAT, "Stat already setup");

    int N=(int)size();
    if (N < MIN_CORRESPONDENCES) return 0;

    TNumber * x_proj_x = new TNumber[N];
    TNumber * x_proj_y = new TNumber[N];
    TNumber * xp_proj_x = new TNumber[N];
    TNumber * xp_proj_y = new TNumber[N];
    TNumber * xp_proj_x45 = 0;
    TNumber * xp_proj_y45 = 0;
    if(b45Deg)
    {
        xp_proj_x45 = new TNumber[N];
        xp_proj_y45 = new TNumber[N];
    }

    TNumber * matchQual = 0;
    if(QUAL_SCORE)
        matchQual = new TNumber[N];

    int nMatch=0;
    iterator ppEnd = end();
    for(iterator pCorrespondence = begin(); pCorrespondence<ppEnd; pCorrespondence++)
    {
        //const CCorrespondence * pCorrespondence = *ppCorrespondence;
        x_proj_x[nMatch]=pCorrespondence->Location1().x();
        x_proj_y[nMatch]=pCorrespondence->Location1().y();
        xp_proj_x[nMatch]=pCorrespondence->Location2().x();
        xp_proj_y[nMatch]=pCorrespondence->Location2().y();
        if (QUAL_SCORE)
            matchQual[nMatch]= 0;//pCorrespondence->Strength();
        THROW( "deprecated")

        if(b45Deg)
        {
        	const double sin45 = 0.707106781;
        	const double cos45 = sin45;
            xp_proj_x45[nMatch] = cos45*xp_proj_x[nMatch] - sin45*xp_proj_y[nMatch];
            xp_proj_y45[nMatch] = sin45*xp_proj_x[nMatch] + cos45*xp_proj_y[nMatch];
        }

        nMatch++;
    }

    scaleProjection(x_proj_x, N);
    scaleProjection(x_proj_y, N);
    scaleProjection(xp_proj_x, N);
    scaleProjection(xp_proj_y, N);
    if(b45Deg)
    {
        scaleProjection(xp_proj_x45, N);
        scaleProjection(xp_proj_y45, N);
    }

    TNumber s1=getProjectionStatInt(x_proj_x, xp_proj_x, xp_proj_y, matchQual, N);
    TNumber s2=getProjectionStatInt(x_proj_y, xp_proj_x, xp_proj_y, matchQual, N);

    if(b45Deg)
    {
		TNumber s3=getProjectionStatInt(x_proj_x, xp_proj_x45, xp_proj_y45, matchQual, N);
		TNumber s4=getProjectionStatInt(x_proj_y, xp_proj_x45, xp_proj_y45, matchQual, N);
		nStat = std::min<TNumber>(std::min<TNumber>(s1,s2), std::min<TNumber>(s3,s4));
    }
    else
    	nStat = std::min<TNumber>(s1,s2);

    delete [] matchQual;
    delete [] x_proj_x;
    delete [] x_proj_y;
    delete [] xp_proj_x;
    delete [] xp_proj_y;

    return nStat;
};

//template<bool QUAL_SCORE>
//inline void CBoWCorrespondences<int, QUAL_SCORE>::scaleProjection(int * axProj, int N)
//{
//    scaleProjection(axProj, INT_SCALE, N);
//};
extern int statVMtune;

templateStat
inline void CCBoWCorrespondences::scaleProjection(TNumber * axProj, int N)
{
/*    if(eVarianceMethod == eStatVar)
        normaliseProjection(axProj, INT_SCALE, N);
    else if(eVarianceMethod == eStatRange)
        scaleProjection(axProj, INT_SCALE, N);
*/
    if(statVMtune)
        normaliseProjection(axProj, INT_SCALE, N);
    else
        scaleProjection(axProj, INT_SCALE, N);

};


/*/Scale to -1..1
void scaleProjection(ColumnVector &v)
{
    v -= v.minimum();
    v /= (v.maximum()*0.5);
    v -= 1;
}

//Get some measure of variance
double getVar(const ColumnVector &v)
{
    return v.SumSquare()/v.nrows();
}

double getProjectionStatInt(const ColumnVector &x_proj_x, const ColumnVector &xp_proj_x, const ColumnVector &xp_proj_y)
{
    double xx_var=getVar(x_proj_x+xp_proj_x);
    double x_x_var=getVar(x_proj_x-xp_proj_x);
    double xy_var=getVar(x_proj_x+xp_proj_y);
    double x_y_var=getVar(x_proj_x-xp_proj_y);

    double min1 = min(xx_var,x_x_var);
    double min2 = min(min1,xy_var);
    double min_var=min(min2,x_y_var);
    //this is recursive and horrible when macros expand but more obvious double min_var=min(min(min(xx_var,x_x_var),xy_var),x_y_var);

    if (min_var==xx_var)
        return xx_var/x_x_var;
    if (min_var==x_x_var)
        return x_x_var/xx_var;
    if (min_var==xy_var)
        return xy_var/x_y_var;

    return x_y_var/xy_var;
}

double getProjectionStat(const TMatchVec &vPossibleMatches)
{
    int nMatches=(int)vPossibleMatches.size();
    if (nMatches<3) return 0;

    ColumnVector x_proj_x(nMatches);
    ColumnVector x_proj_y(nMatches);
    ColumnVector xp_proj_x(nMatches);
    ColumnVector xp_proj_y(nMatches);

    for(int nMatch=0;nMatch<nMatches;nMatch++)
    {
        ColumnVector v = vPossibleMatches[nMatch]->vPoint(1);
        ColumnVector vp = vPossibleMatches[nMatch]->vPoint(2);
        x_proj_x(nMatch+1)=v(1);
        x_proj_y(nMatch+1)=v(2);
        xp_proj_x(nMatch+1)=vp(1);
        xp_proj_y(nMatch+1)=vp(2);
    }
    scaleProjection(x_proj_x);
    scaleProjection(x_proj_y);
    scaleProjection(xp_proj_x);
    scaleProjection(xp_proj_y);

    double s1=getProjectionStatInt(x_proj_x, xp_proj_x, xp_proj_y);
    double s2=getProjectionStatInt(x_proj_y, xp_proj_x, xp_proj_y);

    return min(s1,s2);
}
*/

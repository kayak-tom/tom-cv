/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#ifdef __GNUC__
#define BOW_DLL
#else
#ifdef BOWDLL_EXPORTS
#define BOW_DLL __declspec(dllexport)
#else
#define BOW_DLL __declspec(dllimport)
#endif
#endif

#include "util/location.h"
#include <vector>
#include "util/set2.h"


template<class strengthType = double>
class CBoWMatch
{
protected:
    int nId;
    strengthType dMatchStrength;
public:
    CBoWMatch(const CBoWMatch<strengthType> &m) : nId(m.id()), dMatchStrength(m.MatchStrength()) {};

    CBoWMatch(int nId_in, strengthType dMatchStrength_in) : nId(nId_in), dMatchStrength(dMatchStrength_in) {};

    inline int id() const { return nId; };
    inline strengthType MatchStrength() const { return dMatchStrength; };

    class MatchSort 
    {
    public:
       bool operator( ) ( const CBoWMatch<strengthType> & a,
                     const CBoWMatch<strengthType> & b )
       {
          return a.MatchStrength() > b.MatchStrength() || (a.MatchStrength() == b.MatchStrength() && a.id() > b.id());
       }
    };
};

//typedef std::multiset<CBoWMatch<int>, CBoWMatch<int>::MatchSort> TBoWMatchVector;
typedef set2<CBoWMatch<int>, CBoWMatch<int>::MatchSort/*, CIndividualPool_NoFree_Allocator< CBoWMatch<int>, 128 >*/ > TBoWMatchVector;

class CBoW;
class CDescriptorSet;

class BOW_DLL CBagOfWords
{
    CBoW * pBoW;
public:
    CBagOfWords(bool bAutoCluster, bool bIOwnAddedDSs = true);
    ~CBagOfWords();

    static const int DONT_ADD = -1;

    TBoWMatchVector * getMatches(CDescriptorSet * pDescriptors_in, int nId_in = DONT_ADD);
    TBoWMatchVector * getMatches(int nId_in);
    void addImage(CDescriptorSet * pDescriptors, int nNewImageId);

    void recreateDictionary();

    // Get correspondences between 2 images. Returns all possibilities if there are N matches in one and M in the other
    // nN <= nM
    std::vector<CCorrespondence *> * getCorrespondences(int nId1, int nId2, int nN=1, int nM=1);
};

BOW_DLL CDescriptorSet * getImageDescriptors(const char * szFileName);
BOW_DLL CDescriptorSet * getImageDescriptors(const void * pIplImage_RGB_8U); //Pass in a 'const IplImage * ' to a normal RGB image (3-channel, 8 bit)

//Must use these functions instead of 'delete', as dll has its own address space
BOW_DLL void releaseMatches(TBoWMatchVector ** ppMatches);
BOW_DLL void releaseCorrespondences(std::vector<CCorrespondence *> ** ppCorrespondences);
BOW_DLL void releaseDescriptorSet(CDescriptorSet ** ppDS);
BOW_DLL void getImageDims(int & nWidth, int & nHeight);

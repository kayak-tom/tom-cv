/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * CImageSimulator.cpp
 *
 *  Created on: 10/05/2009
 *      Author: tom
 */

#include "cImageSimulator.h"
#include "geom/geom.h"
#include "util/random.h"
#include "description/vectorDescriptor.h"
#include <set>

using namespace std;

// Simulated descriptor (just an integer, with abs distances between descriptors). Used for simulated data.
class CSymDescriptor : public CDescriptor {
    int id_;
    CLocation loc_;
public:
    CSymDescriptor(CLocation loc, int id) : id_(id), loc_(loc) {
    }
    virtual TDist distance(const CDescriptor * pd) const {
        return abs(id_ - CAST<const CSymDescriptor *>(pd)->id_);
    }
    virtual CLocation location() const {
        return loc_;
    }
    virtual int size() const {
        return sizeof (CSymDescriptor);
    }
    virtual int length() const {
        return 1;
    };
};
class CWorldPoint {
    C3dPoint location3d_;
    int descriptor_;
public:
    CWorldPoint(const C3dPoint & location3d, int desc) : location3d_(location3d), descriptor_(desc) {
    }
    const C3dPoint & loc() const {
        return location3d_;
    }
    int descriptor() const {
        return descriptor_;
    }
};
class CSymCamPos {
    C3dPoint location3d_;
    C3dRotation orientation_;
public:
    CSymCamPos(const C3dPoint & location3d, const C3dRotation & orientation) : location3d_(location3d), orientation_(orientation) {
    }
};
void getCorner(double & nPos, C3dRotation & rot, C3dPoint & pos) {
    if (nPos < 10) {
        rot = C3dRotation();
        pos = C3dPoint(0, 0, nPos);
    } else if (nPos < 15) {
        C3dPoint rotAxis(0, 1, 0);
        C3dPoint firstAxis(1, 0, 0);
        double theta = (nPos - 10) * M_PI / 10; //2*M_PI / 250=0.0251
        double rad = 1;
        pos = C3dPoint(-rad, 0, 10);
        pos += C3dPoint(rad * cos(theta), 0, rad * sin(theta));
        rot = C3dRotation(rotAxis, theta);
    } else {
        C3dPoint rotAxis(0, 1, 0);
        double theta = M_PI / 2;
        pos = C3dPoint(14 - nPos, 0, 11);
        rot = C3dRotation(rotAxis, theta);
    }
    cout << pos << endl;
}
void getCircle(double & nPos, C3dRotation & rot, C3dPoint & pos) {
    C3dPoint rotAxis(0, 1, 0);
    C3dPoint firstAxis(1, 0, 0);
    double theta = nPos * 2 * M_PI / 300.123; //2*M_PI / 250=0.0251
    double rad = 25;
    pos = C3dPoint(rad * cos(theta), 0, rad * sin(theta)); //This is unclear as we start with z forewards
    C3dRotation rot1(rotAxis, theta);
    cout << rot1 << endl;
    C3dPoint rotAxisInPlane = rot1.t() * firstAxis;
    C3dRotation rot2; //(rotAxisInPlane, 0.5*(theta-floor(theta+.25))); //todo why isn't this working?
    //cout << rot2 << endl;
    rot = rot2; //*rot1;
    cout << rot << endl;
}
void getCircle2(double & nPos, C3dRotation & rot, C3dPoint & pos) {
    C3dPoint rotAxis(0, 1, 0);
    C3dPoint firstAxis(1, 0, 0);
    double theta = nPos * 2 * M_PI / 300.123; //2*M_PI / 250=0.0251
    double rad = 25;
    pos = C3dPoint(rad * cos(theta), 0, rad * sin(theta)); //This is unclear as we start with z forewards

    theta -= floor(theta / (2 * M_PI))*2 * M_PI;
    if (theta > M_PI)
        theta = 2 * M_PI - theta;

    C3dRotation rot1(firstAxis, 0.4 * theta);
    cout << rot1 << endl;
    C3dPoint rotAxisInPlane = rot1.t() * firstAxis;
    C3dRotation rot2; //(rotAxisInPlane, 0.5*(theta-floor(theta+.25))); //todo why isn't this working?
    cout << rot2 << endl;
    rot = rot2*rot1;
    cout << rot << endl;
}
void getStationary(const double & nPos, C3dRotation & rot, C3dPoint & pos) {
    double dStatPos = nPos < 10 ? nPos : 10;
    if (nPos > 20)
        dStatPos = nPos - 10;
    if (nPos > 30)
        dStatPos = 20;
    if (nPos > 40)
        dStatPos = nPos - 20;
    pos = C3dPoint(0, 0, dStatPos);
    C3dPoint vertical(0, 1, 0);
    if (nPos < 25)
        rot = C3dRotation();
    else
        rot = C3dRotation(vertical, (nPos - 25)*0.1);
}
void getSquare(double & nPos, C3dRotation & rot, C3dPoint & pos) {
    double squarePos = nPos - 400.123 * floor(nPos / 400);
    /*if(squarePos < 80)
            pos = C3dPoint(0, 0, squarePos);
    else if(squarePos < 200)
    {
            double sidePos=squarePos-80;
            pos = C3dPoint(0.666*sidePos, 0, 80);
    }
    else if(squarePos < 260)
    {
            double sidePos=squarePos-200;
            pos = C3dPoint(80, 0, 80-1.333*sidePos);
    }
    else
    {
            double sidePos=squarePos-260;
            pos = C3dPoint(80-sidePos*(80./140.), 0, 0);
    }*/
    if (squarePos < 100)
        pos = C3dPoint(0, 0, 0.8 * squarePos);
    else if (squarePos < 200) {
        double sidePos = squarePos - 100;
        pos = C3dPoint(0.8 * sidePos, 0, 80);
    } else if (squarePos < 300) {
        double sidePos = squarePos - 200;
        pos = C3dPoint(80, 0, 80 - 0.8 * sidePos);
    } else {
        double sidePos = squarePos - 300;
        pos = C3dPoint(80 - sidePos * 0.8, 0, 0);
    }

    double theta = ((int) (nPos) % 100) * 2 * M_PI / 119.123; //2*M_PI / 250=0.0251
    if (theta > M_PI)theta = 2 * M_PI - theta;
    C3dPoint rotAxis(0, 1, 0);
    C3dPoint firstAxis(sin(theta * 0.7 + 0.2), sin(theta * 0.4 + 0.3), 1);

    C3dRotation rot1(firstAxis, .5 * cos(1.69 * theta) * sin(theta)); //Make sure both 0 at 0
    //rot = rot1 * C3dRotation(rotAxis, .9*sin(theta)*cos(.69*(theta+0.4)));
    //rot = rot1 * C3dRotation(rotAxis, theta);
    rot = rot1 * C3dRotation(rotAxis, 0.5 * theta);
    //rot = C3dRotation();
    cout << rot << endl;
}
CImageSimulator::CImageSimulator(const CCamCalibMatrix & K, const CDescriptorSetClusteringParams & DSC_PARAMS) {
    const int IM_HEIGHT = 600, IM_WIDTH = 600;
    // Generate the world
    double xMin = -50, xMax = 130, yMin = -50, yMax = 50, zMin = -50, zMax = 130, maxDist = sqr(25);

    vector<CWorldPoint> aWorldPoints;

    for (int nPoint = 0; nPoint < 1000000; nPoint++) //100k is enough
    {
        C3dPoint p(CRandom::Uniform(xMin, xMax), CRandom::Uniform(yMin, yMax), CRandom::Uniform(zMin, zMax));
        //CWorldPoint wp(p, nPoint % 5000);
        int nPointId = nPoint / 10 + CRandom::Uniform(5000); //Introduces 'self similar' features
        CWorldPoint wp(p, nPoint);
        aWorldPoints.push_back(wp);
    }

    // Generate the camera trajectory
    vector <CSymCamPos> aCamPositions;

    // Generate descriptor sets

    //double * adRDcoeffs = 0;
    //int nRDcoeffs = 0;
    //Matrix K(3,3);geom
    //getCalibrationMat(K, &adRDcoeffs, nRDcoeffs);

    double dNoise = 0 * 1.; // 3 works ok

    double speed = 5;
    const int nStepsAtEachSpeed = 100;
    int nStepsAtThisSpeed = 0;
    double aSpeeds[8] = {1.3, 1.7, 3.3, 0.4, 5.1, 2.2, 1, 3.1};
    for (double nPos = 0; nPos < 500; nPos += speed, nStepsAtThisSpeed++) {
        if (nStepsAtThisSpeed >= nStepsAtEachSpeed) {
            nStepsAtThisSpeed = 0;
            speed = aSpeeds[(int) nPos % 8];
        }

        const bool bIntroduceGap = false;
        if (bIntroduceGap) {
            if (nPos > 100 && nPos < 200) {
                CDescriptorSet * pDS = new CMetricSpaceDescriptorSet(DSC_PARAMS, 0);
                aSimPhotos.push_back(pDS);
                int i = 0;
                for (vector<CWorldPoint>::const_iterator pWP = aWorldPoints.begin(); pWP != aWorldPoints.end() && pDS->Count() < 100; pWP++) {
                    i++; //Don't want 2 different descriptors in the same place
                    pDS->Push(new CSymDescriptor(CLocation(i + 10, i/*CRandom::Uniform(1000) + 10*/), pWP->descriptor()));
                }
                continue;
                //nPos = 200;
            }
        }

        C3dPoint pos;
        C3dRotation rot;
        //getCircle(nPos, rot, pos);
        //getCorner(nPos, rot, pos);
        getSquare(nPos, rot, pos);
        //getStationary(nPos, rot, pos);

        CSymCamPos cam(pos, rot);
        aCamPositions.push_back(cam); //Later these will be used for validation

        CCamera P = (rot | -(rot * pos));

        //cout << P << "=P--Camera\n";
        P.calibrate(K);
        //cout << P << "=P calibrated\n";
        //cout << -pos * (1.0/pos.length()) << "=T_ba_b_dir";

        set<CLocation> locationsUsed;
        CDescriptorSet * pDS = new CMetricSpaceDescriptorSet(DSC_PARAMS, 0);
        aSimPhotos.push_back(pDS);

        for (vector<CWorldPoint>::const_iterator pWP = aWorldPoints.begin(); pWP != aWorldPoints.end(); pWP++) {
            C3dPoint worldPointLoc = pWP->loc();

            if ((worldPointLoc - pos).sum_square() < maxDist) {
                if (worldPointLoc.testInFront(P)) {
                    //point is in front of cam
                    C2dPoint camPointExact = worldPointLoc.photo(P);
                    double x = camPointExact.getX() + CRandom::Uniform(-dNoise, dNoise);
                    if (CRandom::Uniform(1.0) > 0.5)
                        x += 25;

                    double y = camPointExact.getY() + CRandom::Uniform(-dNoise, dNoise);
                    CLocation camLoc(x, y);
                    if (camLoc.x() > 0 && camLoc.y() > 0 && camLoc.x() < IM_WIDTH && camLoc.y() < IM_HEIGHT) {
                        //Point lies within image.
                        if (locationsUsed.find(camLoc) == locationsUsed.end()) {
                            locationsUsed.insert(camLoc);
                            pDS->Push(new CSymDescriptor(camLoc, pWP->descriptor()));
                        }
                    }
                }
            }
        }
        cout << pDS->Count() << " descriptors in image " << nPos << endl;
    }
}
CImageSimulator::~CImageSimulator() {
    for (int i = 0; i < (int) aSimPhotos.size(); i++) {
        if (aSimPhotos[i]) {
            for (int j = 0; j < (int) aSimPhotos[i]->Count(); j++)
                delete (*(aSimPhotos[i]))[j];
            delete aSimPhotos[i];
        }
    }
}
CDescriptorSet * CImageSimulator::getSimDescriptors(int id) {
    if (id >= (int) aSimPhotos.size()) return 0;
    CDescriptorSet * pDS = aSimPhotos[id];
    aSimPhotos[id] = 0;
    return pDS;
}

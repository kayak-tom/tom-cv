/* 
 * File:   makeBasis_GramSchmidt.h
 * Author: tom
 *
 * Created on 19 April 2011, 2:14 PM
 */

#ifndef MAKEBASIS_GRAMSCHMIDT_H
#define	MAKEBASIS_GRAMSCHMIDT_H

template<typename M>
void makeClosestE(M & E, const bool bMakeF = false) {
#if EIGEN_VERSION_AT_LEAST(2,90,0)
    Eigen::JacobiSVD<M> svdE(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
#else
    Eigen::SVD<M> svdE(E);
#endif

    M D;
    D.setIdentity();
    if(bMakeF)
    {
        D(0, 0) = svdE.singularValues()(0);
        D(1, 1) = svdE.singularValues()(1);
    }
    D(2, 2) = 0;

    E = svdE.matrixU() * D * svdE.matrixV().transpose();

    CHECKNAN(E.sum());
}

template<typename TFloat>
TFloat getETEResid(const Eigen::Matrix<TFloat, 3, 3> & E)
{
    Eigen::Matrix<TFloat, 3, 3> EET = E*E.transpose();
    EET.diagonal().array() -= 0.5*EET.trace();
    EET *= E;
    return EET.squaredNorm();
}

template<typename TFloat>
TFloat getResids(const Eigen::Matrix<TFloat, 3, 3> & E, const Eigen::Matrix<TFloat, 3, 1> * am0, const Eigen::Matrix<TFloat, 1, 3> * am1)
{
    TFloat dErrSq = 0;
    for (int i = 0; i < 5; i++) {
        dErrSq += am1[i] * E * am0[i];
    }
    return dErrSq;
}


#define makeOrth(m1, m2) \
{    m1 -= m1.dot(m2) * m2;    if(IS_DEBUG) CHECK(!zero(m1.dot(m2)), "makeOrth failed"); }


//template<typename TMat> inline void makeNormal(TMat & M)
#define makeNormal(M) \
{  M /= M.stableNorm(); }

//Similar to the method recommended in Nister (04)
template<int NUM_POINTS>
void makeBasisForE_GramSchmidt(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, Eigen::Matrix<double, 9, 9-NUM_POINTS > & EE) {
    const int BASIS_SIZE = 9-NUM_POINTS;
    Eigen::Matrix<double, 1, NUM_POINTS> Q1_col1;
    Eigen::Matrix<double, 1, NUM_POINTS> Q2_col1;
    Eigen::Matrix<double, 1, NUM_POINTS> Q1_col2;
    Eigen::Matrix<double, 1, NUM_POINTS> Q2_col2;

    Eigen::Matrix< double, 9, NUM_POINTS > Q; 

    for (int i = 0; i < NUM_POINTS; i++) {
        Q1_col1(i) = m0[anHypSet[i]].getX();
        Q1_col2(i) = m0[anHypSet[i]].getY();
        Q2_col1(i) = m1[anHypSet[i]].getX();
        Q2_col2(i) = m1[anHypSet[i]].getY();

        //And put 1s in the last row
        Q(8, i) = 1;
    }

    Q.row(0) = Q1_col1.array() * Q2_col1.array();
    Q.row(1) = Q1_col2.array() * Q2_col1.array();
    Q.row(2) = Q2_col1;
    Q.row(3) = Q1_col1.array() * Q2_col2.array();
    Q.row(4) = Q1_col2.array() * Q2_col2.array();
    Q.row(5) = Q2_col2;
    Q.row(6) = Q1_col1;
    Q.row(7) = Q1_col2;
    //already Q.row(8) = 1;
    
    //First make cols of Q orthogonal:
    for(int nEq = 0; nEq<NUM_POINTS; nEq++)
    {
        for(int nEqAlreadyOrth = 0; nEqAlreadyOrth<nEq; nEqAlreadyOrth++)
             makeOrth(Q.col(nEq), Q.col(nEqAlreadyOrth));

        makeNormal(Q.col(nEq));
    }    
    
    //cout << Q.transpose() * Q << endl;
    
    //4 cols of E should be 4 vectors orth. to the 5 cols of Q
    //Eigen is col-major
    //EE.setRandom();
    EE.setIdentity();
    for(int nNewBasis = 0; nNewBasis<BASIS_SIZE; nNewBasis++)
    {
        for(int nEq = 0; nEq<NUM_POINTS; nEq++)
            makeOrth(EE.col(nNewBasis), Q.col(nEq));
        for(int nEqAlreadyOrth = 0; nEqAlreadyOrth<nNewBasis; nEqAlreadyOrth++)
            makeOrth(EE.col(nNewBasis), EE.col(nEqAlreadyOrth));
        
        makeNormal(EE.col(nNewBasis));
    }  
    //cout << EE.transpose() * EE << endl;
}

//Use 10D vectors which will vectorise well in Eigen:
template<int NUM_POINTS>
void makeBasisForE_GramSchmidt_Vectorise(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, Eigen::Matrix<double, 9, 9-NUM_POINTS > & EE) {
    const int BASIS_SIZE = 9-NUM_POINTS;
    Eigen::Matrix<double, 1, NUM_POINTS> Q1_col1;
    Eigen::Matrix<double, 1, NUM_POINTS> Q2_col1;
    Eigen::Matrix<double, 1, NUM_POINTS> Q1_col2;
    Eigen::Matrix<double, 1, NUM_POINTS> Q2_col2;

    Eigen::Matrix< double, 10, NUM_POINTS > Q;

    for (int i = 0; i < NUM_POINTS; i++) {
        Q1_col1(i) = m0[anHypSet[i]].getX();
        Q1_col2(i) = m0[anHypSet[i]].getY();
        Q2_col1(i) = m1[anHypSet[i]].getX();
        Q2_col2(i) = m1[anHypSet[i]].getY();

        //And put 1s in the last row
        Q(8, i) = 1;
        
        Q(9, i) = 0;
    }

    Q.row(0) = Q1_col1.array() * Q2_col1.array();
    Q.row(1) = Q1_col2.array() * Q2_col1.array();
    Q.row(2) = Q2_col1;
    Q.row(3) = Q1_col1.array() * Q2_col2.array();
    Q.row(4) = Q1_col2.array() * Q2_col2.array();
    Q.row(5) = Q2_col2;
    Q.row(6) = Q1_col1;
    Q.row(7) = Q1_col2;
    //already Q.row(8) = 1;
    
    //First make cols of Q orthogonal:
    for(int nEq = 0; nEq<NUM_POINTS; nEq++)
    {
        for(int nEqAlreadyOrth = 0; nEqAlreadyOrth<nEq; nEqAlreadyOrth++)
             makeOrth(Q.col(nEq), Q.col(nEqAlreadyOrth));

        makeNormal(Q.col(nEq));
    }    
    
    //cout << Q.transpose() * Q << endl;
    
    //4 cols of E should be 4 vectors orth. to the 5 cols of Q
    //Eigen is col-major
    Eigen::Matrix<double, 10, 9-NUM_POINTS > EE_aligned;
    
    EE_aligned.setIdentity();
    for(int nNewBasis = 0; nNewBasis<BASIS_SIZE; nNewBasis++)
    {
        for(int nEq = 0; nEq<NUM_POINTS; nEq++)
            makeOrth(EE_aligned.col(nNewBasis), Q.col(nEq));
        for(int nEqAlreadyOrth = 0; nEqAlreadyOrth<nNewBasis; nEqAlreadyOrth++)
            makeOrth(EE_aligned.col(nNewBasis), EE_aligned.col(nEqAlreadyOrth));
        
        makeNormal(EE_aligned.col(nNewBasis));
    }  
    
    //EE = EE_aligned.block<9, 9-NUM_POINTS>(0,0);
    for(int c=0;c<9-NUM_POINTS;c++)
        for(int r=0;r<9;r++)
            EE(r,c) = EE_aligned(r,c);
    
    //cout << EE.transpose() * EE << endl;
}

#endif	/* MAKEBASIS_GRAMSCHMIDT_H */


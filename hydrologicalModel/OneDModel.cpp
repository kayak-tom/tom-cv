#include <iostream>
#include <Eigen/Core>

class CTotalFlowConstraint;
std::vector<CTotalFlowConstraint*> g_aTotalFlowConstraints;

#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/block_solver.h> 
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/dense/linear_solver_dense.h>

#include "subcatchment.h"


g2o::SparseOptimizer optimizer; //Make optimiser global for now (todo: singleton)

void testOptimiser()
{
    CFlowVertex* pRainfallEstimateOneTimestep = new CFlowVertex(0); //pRainfallEstimateOneTimestep->setId(0);
    
    /*CFlowMeasurement * pRainGaugeMeasurement = */new CFlowMeasurement(pRainfallEstimateOneTimestep, 10, 100);
    /*CFlowMeasurement * pRainGaugeMeasurement2 = */new CFlowMeasurement(pRainfallEstimateOneTimestep, 11, 100);

    CStoreElement* pSurfaceWaterOneTimestep = new CStoreElement(0); //pSurfaceWaterOneTimestep->setId(1);
    CStoreElement* pSurfaceWaterNextTimestep = new CStoreElement(1); //pSurfaceWaterNextTimestep->setId(2);
    CStoreElement* pSurfaceWaterNextNextTimestep = new CStoreElement(2); //pSurfaceWaterNextTimestep->setId(3);
    
    CTotalFlowConstraint * pConserveFlow = new CTotalFlowConstraint(pSurfaceWaterOneTimestep, pSurfaceWaterNextTimestep);
    pConserveFlow->addInflow(pRainfallEstimateOneTimestep);
    
    //Initialise it to 'near 0'
    /*CStoreInitialisingAssumption * pInitAssumption =*/ new CStoreInitialisingAssumption(pSurfaceWaterOneTimestep);
    
    
    
    
    CDrainageFunction * pLinearSWOutflow= new CLinearDrainageFunction(0.3);
    CFlowVertex* pSWOutFlow = new CFlowVertex(1);//pSWOutFlow->setId(4);
    /*COutFlowConstraint * pSWOutFlowConstr =*/ new COutFlowConstraint(pSurfaceWaterNextTimestep, pSWOutFlow, pLinearSWOutflow);
    
    CTotalFlowConstraint * pConserveFlow2 = new CTotalFlowConstraint(pSurfaceWaterNextTimestep, pSurfaceWaterNextNextTimestep);
    pConserveFlow2->addOutflow(pSWOutFlow);
    
    /* vertex with temp estimate */CPETVertex* pTemperatureEstimateOneTimestep = new CPETVertex(0);
    
    const double dTemperatureSD = 0.5;
    //MUST have a measurement for every estimate, otherwise MLE will be unconstrained (could add forward-backward constraints instead?). Each CTemperatureMeasurement constrains a temperature estimate
    /* unary edge constraining temp estimate CTemperatureMeasurement * pTemperatureMeasurement =*/  new CPETMeasurement(pTemperatureEstimateOneTimestep, 20, sqr(dTemperatureSD));

    //For each flow vertex there is a ternary edge constraint on the flow vertex, the SW store and the temperature
    /* flow out vertex, with associated 3-edge  CFlowVertex* pETEstimateOneTimestep =* / new CFlowVertex(0); 
	CETDrainageFn etFunction(0.02);
	/ *CTemperatureDependentOutFlowElement * pElement =* / new CTemperatureDependentOutFlowElement(pSurfaceWaterOneTimestep, pTemperatureEstimateOneTimestep, &etFunction);*/
    throw "Not likely to work with line above removed";

    cout << "Optimizing" << endl;
    optimizer.initializeOptimization();
    optimizer.optimize(20);
    cout << "done." << endl; 

}

CCatchmentBase * setupSimCatchment() //Simple simulation with 3 subcatchment structure, 1000 timesteps, temperature 15 for 500 then 30 for 500.
{
    
    double C_c_m3 = 100; double c_t = 0.01; //parameters for ET from canopy
    double theta_pa_m3 = 10;
    
    const std::vector<CSubcatchmentParams> aSubcatchmentsParams =
    { CSubcatchmentParams("northBranch", 5, C_c_m3, c_t, theta_pa_m3), CSubcatchmentParams("southBranch", 3, C_c_m3, c_t, theta_pa_m3), CSubcatchmentParams("Downstream", 4,  C_c_m3, c_t, theta_pa_m3) };
    
    CCatchment * pCatchment = new CCatchment(5*60, aSubcatchmentsParams);
    
    for(int i=0; i<1000; i++)
    {
        double dRainfall = 0;
        if(i % 20 == 0)
            dRainfall = 5;
            
        double PET = 15;
        if( i > 500)
            PET = 30;
        
        pCatchment->addTimestep(i, i, dRainfall, PET);
    }
    //pCatchment->addFlow(298, 1);
    //pCatchment->addFlow(2, 3);
    
    return pCatchment;
}

int main(int argc, char** argv)
{
    //testOptimiser();

    //typedef BlockSolver< BlockSolverTraits<-1, -1> >  HydroBlockSolver;
    //typedef LinearSolverCholmod<HydroBlockSolver::PoseMatrixType> HydroLinearSolver;

    //SparseOptimizer optimizer;
    /*HydroLinearSolver* linearSolver = new HydroLinearSolver();
    linearSolver->setBlockOrdering(false);
    HydroBlockSolver* solver = new HydroBlockSolver(&optimizer, linearSolver);
    optimizer.setAlgorithm(solver);
    optimizer.setVerbose(true);*/
    g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType> * pLinearSolverCholmod = new g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>();
    pLinearSolverCholmod->setBlockOrdering(false); 
    g2o::BlockSolverX * pSolver ( new g2o::BlockSolverX(pLinearSolverCholmod) );
    g2o::OptimizationAlgorithmLevenberg * solver(new g2o::OptimizationAlgorithmLevenberg(pSolver));
    optimizer.setAlgorithm(solver);
    
    
    //Add all the nodes/edges to optimiser:
    CCatchmentBase * pCatchment = setupSimCatchment();
    //CCatchmentBase * pCatchment = setupSimTopnetCatchment();


    cout << "Optimizing" << endl;
    optimizer.initializeOptimization();
    optimizer.optimize(10);
        
    cout << "done." << endl; 
    
    std::ofstream outputFile("results.tsv");
    pCatchment->pp(outputFile);
    outputFile.close();
    
    return 0;
}

#ifndef _SIMONESEED_H
#define _SIMONESEED_H
#include <random>
#include <fstream>
#include <time.h>
#include <vector>
#include "Eigen/Core"
#include "Eigen/SparseCore" // Eigen; added by right click on project -> Properties -> C/C++ -> Additional Include Directories
#include "Elements.h"
#include "Couplings.h"
#include "Networks.h"
#include "Misc.h"
#include "Input.h"
#include "Distribution_c.h"

namespace Distribution
{
    class IDistribution;
}

namespace Elements
{
    struct ElementParam_t;
	class IElement;
}

namespace Couplings
{
    struct CouplingParam_t;
    class ICoupling;
}

namespace Networks
{
    struct AdjacencyMatrixParam_t;
}

namespace Misc
{
    struct AvgVar_t;
}

namespace Input
{
    struct InputParam_t;
    enum class SimulationType;
    enum class SeedType;
}

namespace Simulations
{
    /*struct OneSeedInputParam_t
    {
        int N;
        int nSim;
        Elements::ElementParam_t nodeParam;
        Couplings::CouplingParam_t edgeParam;
        Networks::AdjacencyMatrixParam_t netParam;
    };

    extern OneSeedInputParam_t param;//*/
    class OneSeedSpreading
    {
    private:
        typedef void(OneSeedSpreading::*TimeStepAndRecordNet_ptr)();
        typedef void(OneSeedSpreading::*TimeStepAndRecordNodes_ptr)();
        typedef void(OneSeedSpreading::*SpreadAllToAll_ptr)();

        std::mt19937_64& rand;

        std::unique_ptr<Distribution::IDistribution> dist;

        Input::InputParam_t param;

        std::clock_t initSimTime;
        std::clock_t endSimTime;

        std::vector<Elements::IElement*> nodes; // collection of network nodes
        std::vector<Couplings::ICoupling*> edges; // collection of network edges

        int timeStep; // current time step of the simulation
        int tTotal; // total time of the simulation (for EvolveAndRecord function)
        bool isNetFullyActive; // if all the nodes are active (during dynamics), this flag changes to true

        Eigen::ArrayXXd nodesData; // each column is the time evolution of each node
        Eigen::ArrayXXd spreadingTimeAvg; // spreading time from each node (first vector) to every other node (internal vector)
        Eigen::ArrayXXd spreadingTimeVar; // spreading time from each node (first vector) to every other node (internal vector)
        Eigen::ArrayXd netDynData; // stores the number of firings for each timestep of the simulation
        Eigen::ArrayXd totalNetSimData; // stores the total number of firings for each simulation
        Eigen::ArrayXXd thetaRhoChiData; // 3 cols Vs ntheta rows to store theta, rho and chi (var of rho)
        std::vector<double> totalNetSimDataPT; // stores the total number of firings for each simulation, in order to calculate average and variance

        std::string outputFileNodesName;
        std::string outputFileSpAvgName;
        std::string outputFileSpVarName;
        std::string outputFileNetDynName;
        std::string outputFileTotNetSimName;
        std::string outputFilePhaseTransName;
        std::fstream outputFileNodes; // output file stream if input.writeOnRun is set
        std::fstream outputFileSpAvg; // output file stream if input.writeOnRun is set
        std::fstream outputFileSpVar; // output file stream if input.writeOnRun is set
        std::fstream outputFileNetDyn;
        std::fstream outputFileTotNetSim;
        std::fstream outputFilePhaseTrans;

        void InitializeSimulation(Input::InputParam_t); // parses input params
        void CreateNodes(); // creates nodes
        Eigen::SparseMatrix<double> CreateAdjacencyMatrix(); // creates adjacency matrix
        void SetupNetwork(); // creates edges and connect nodes
        void SetInitialCondition(int);
        void ClearNetworkData();
        void ResetNodesAndEdges();
        std::string GetOutputFileName(std::string);

        // evolve and record nodes activity
        void PrepareToRecordNodes(int); // initializes nodes data variable to record the desired amount of time steps
        void TimeStepAndRecordNodesInMemory(); // time step of the simulation and record nodes data in RAM memory
        void TimeStepAndRecordNodesInFile(); // time step of the simulation and record nodes data in output file
        TimeStepAndRecordNodes_ptr TimeStepAndRecordNodes;
        void EvolveAndRecordNodesOnce(int, int, int); // seed and tTotal // evaluate the specified number of time steps of the network and records nodes data
        void WriteNodesData(int);
        void SetNodesThreshold(double);

        // evolve and record network activity
        void PrepareToRecordNet(int); // initializes nodes data variable to record the desired amount of time steps
        void TimeStepAndRecordNetInMemory(); // time step of the simulation and record nodes data in RAM memory
        void TimeStepAndRecordNetInFile(); // time step of the simulation and record nodes data in output file
        TimeStepAndRecordNet_ptr TimeStepAndRecordNet;
        void EvolveAndRecordNetworkOnce(int, int, int); // runs the simulation
        void WriteNetworkData(int);
        void WriteTotalNetworkData();

        // general simulation
        void TimeStep(); // time step of the simulation
        void EvolveNetwork(int, int); // evolves network for the required amount of timesteps
        double SumNodesStates();

        // spreading simulation
        int SpreadNodeToNode(int, int); // seedNode, targetNode // spreads activity from seedNode to targetNode
        Misc::AvgVar_t SpreadNodeToNodeManyTimes(int, int); // seedNode, targetNode // spreads activity from seedNode to targetNode nSim times
        SpreadAllToAll_ptr SpreadAllToAll_int; // runs the simulation
        void SpreadAllToAllInMemory(); // runs the simulation
        void SpreadAllToAllInFile(); // runs the simulation
        void WriteSpreadingTimeMatrix(); // writes spreadingTime matrix

        // phase transition
        Misc::AvgVar_t RunSingleTheta(double);
        void WritePhaseTransData();

    public:
        OneSeedSpreading(Input::InputParam_t, std::mt19937_64&);
        ~OneSeedSpreading();
        void EvolveAndRecordNodes(int, int);
        void SpreadAllToAll(); // runs the simulation
        void EvolveAndRecordNetwork(int, int); // runs the simulation
        void PhaseTransition();
    };
}

#endif
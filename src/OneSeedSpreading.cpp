#include <iostream>
#include <fstream>
#include <time.h>
#include <vector>
#include "Eigen/Core"
#include "Eigen/SparseCore" // Eigen; added by right click on project -> Properties -> C/C++ -> Additional Include Directories
#include "Elements.h"
#include "Couplings.h"
#include "Networks.h"
#include "Misc.h"
#include "OneSeedSpreading.h"
#include "Input.h"
#include "Output.h"

namespace Simulations
{
    /*OneSeedInputParam_t param = {
        -1, -1,
        {
            Elements::ElementType::SI, -1, 0.0, 0.0, -1
        },
        {
            Couplings::CouplingType::PulseCoupling, (*(new Elements::SIElement(-1, 0.0, 0.0))), (*(new Elements::SIElement(-1, 0.0, 0.0))), 0.0
        },
        {
            Networks::AdjacencyMatrixType::RandomGraph, -1, std::vector<int>(), -1, -1, 0.0, false, false, ""
        }
    };//*/

    OneSeedSpreading::OneSeedSpreading(Input::InputParam_t input, std::mt19937_64& rand)
        : rand(rand)
    {
        this->InitializeSimulation(input);
        std::cout << "OneSeedSimulation is prepared..." << std::endl;
    }
    OneSeedSpreading::~OneSeedSpreading()
    {
        int i, n = this->edges.size();
        for (i = 0; i < n; i++)
            delete this->edges[i];
        for (i = 0; i < this->param.N; i++)
            delete this->nodes[i];
    }

    // parses input params
    void OneSeedSpreading::InitializeSimulation(Input::InputParam_t input)
    {
        OneSeedSpreading::param = input;
        if (this->param.writeOnRun)
        {
            this->TimeStepAndRecordNet = &OneSeedSpreading::TimeStepAndRecordNetInFile;
            this->TimeStepAndRecordNodes = &OneSeedSpreading::TimeStepAndRecordNodesInFile;
            this->SpreadAllToAll_int = &OneSeedSpreading::SpreadAllToAllInFile;
        }
        else
        {
            this->TimeStepAndRecordNet = &OneSeedSpreading::TimeStepAndRecordNetInMemory;
            this->TimeStepAndRecordNodes = &OneSeedSpreading::TimeStepAndRecordNodesInMemory;
            this->SpreadAllToAll_int = &OneSeedSpreading::SpreadAllToAllInMemory;
        }
        this->isNetFullyActive = false;
        this->dist = Distribution::GetDistribution(Distribution::GetDefaultDistributionParam(Distribution::DistributionType::Uniform), rand);
        SetupNetwork();
    }

    // creates nodes
    void OneSeedSpreading::CreateNodes()
    {
        Elements::ElementParam_t p = {};
        p.ind = -1;
        p.nRefStates = param.nRefStates;
        p.threshold = param.theta;
        p.type = param.nodeType;
        p.s0 = 0.0;
        nodes = std::vector<Elements::IElement*>(param.N);
        int i = 0;
        while (i < param.N)
        {
            p.ind = i;
            nodes[i] = Elements::GetElement(p);
            i++;
        }
    }

    // creates nodes
    void OneSeedSpreading::SetupNetwork()
    {
        Eigen::SparseMatrix<double> A = CreateAdjacencyMatrix(); // matrix of col vectors
        CreateNodes();
        edges = std::vector<Couplings::ICoupling*>(A.nonZeros());
        int i = 0;
        for (int k = 0; k < A.outerSize(); ++k) // k-th column - copied from https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) // for each nonzero element of the k-th column
            {
                Couplings::CouplingParam_t p = {
                    param.coupType, *nodes[it.row()], *nodes[k], it.value()
                };
                edges[i] = Couplings::GetCoupling(p);
                nodes[k]->AddNeighbor(*edges[i]); // adds the last edge (just created) to the in-neighbors of the k-th node
                i++;
                /*it.value();
                it.row();   // row index
                it.col();   // col index (here it is equal to k)
                it.index(); // inner index, here it is equal to it.row()//*/
            }
            (nodes[k]->*(nodes[k]->SetTotalNeighWeight))();
        }
    }

    // creates AdjacencyMatrix
    Eigen::SparseMatrix<double> OneSeedSpreading::CreateAdjacencyMatrix()
    {
        Networks::AdjacencyMatrixParam_t pm;
        pm.type = param.netType;
        pm.nElems = param.N;
        pm.isDirected = param.isDirected;
        pm.isPeriodic = param.isPeriodic;
        pm.rewiringProb = param.rewiringProb;
        pm.numOfEdgesForNewElem = param.numOfEdgesBA;
        pm.numOfNeighbours = param.numOfNeighWS;
        pm.netFileName = param.netFile;
        pm.L = std::vector<int>(3); pm.L.push_back(param.Lx); pm.L.push_back(param.Ly); pm.L.push_back(param.Lz);
        pm.weightAmp = param.weightAmp;
        pm.weightType = param.weightType;
        Eigen::SparseMatrix<double> mat = Networks::GetAdjacencyMatrix(pm, this->rand)->BuildAndGetMatrix();
        if (this->param.netType == Networks::AdjacencyMatrixType::FromFile)
        {
            this->param.N = mat.rows();
        }
        return mat;
    }

    void OneSeedSpreading::PrepareToRecordNodes(int tTotal)
    {
        int i;
        this->param.tTotal = tTotal;
        if (this->param.writeOnRun)
        {
            this->outputFileNodes << this->timeStep;
            i = 0;
            while (i < param.N)
            {
                this->outputFileNodes << " " << this->nodes[i]->GetState();
                i++;
            }
            this->outputFileNodes << std::endl;
        }
        else
        {
            this->nodesData = Eigen::ArrayXXd(tTotal + 1, this->param.N);
            i = 0;
            while (i < param.N)
            {
                this->nodesData(0, i) = nodes[i]->GetState();
                i++;
            }
        }
    }

    void OneSeedSpreading::TimeStepAndRecordNodesInFile()
    {
        int i = 0, n = this->edges.size();
        while (i < n)
        {
            this->edges[i]->Step();
            i++;
        }
        this->outputFileNodes << this->timeStep;
        i = 0;
        while (i < param.N)
        {
            this->nodes[i]->Step();
            this->outputFileNodes << " " << this->nodes[i]->GetState();
            i++;
        }
        this->outputFileNodes << std::endl;
    }

    // time step of the simulation
    void OneSeedSpreading::TimeStepAndRecordNodesInMemory()
    {
        int i = 0, n = this->edges.size();
        while (i < n)
        {
            this->edges[i]->Step();
            i++;
        }
        i = 0;
        while (i < param.N)
        {
            this->nodes[i]->Step();
            this->nodesData(timeStep, i) = this->nodes[i]->GetState();
            i++;
        }
    }
    
    // evaluate the specified number of time steps of the network
    void OneSeedSpreading::EvolveAndRecordNodesOnce(int seed, int tTotal, int simInd)
    {
        if (this->param.writeOnRun)
        {
            //this->outputFileNodesName = "1s_dyn_nSim" + std::stringstream(simInd).str() + ".dat";
            this->outputFileNodesName = this->GetOutputFileName(std::stringstream(simInd).str());
            this->outputFileNodes = Output::OpenFile(this->outputFileNodesName, Input::GetInputParamString(param) + "# t    nodes");
        }

        this->SetInitialCondition(seed);

        int i = 0;
        this->timeStep = 0;
        this->PrepareToRecordNodes(tTotal);
        //std::cout << nodesData.row(timeStep) << std::endl;
        //std::cin.get();
        this->timeStep = 1;
        while (this->timeStep <= tTotal)
        {
            (this->*TimeStepAndRecordNodes)();
            //std::cout << nodesData.row(timeStep) << std::endl;
            //std::cin.get();
            this->timeStep++;
        }
        if (this->param.writeOnRun)
        {
            
            this->outputFileNodes.close();
            //Output::WriteFirstLine(this->outputFileNodesName, "# " + Input::inputCmdLine + "\n" + simTime);
        }
    }

    void OneSeedSpreading::EvolveAndRecordNodes(int seed, int tTotal)
    {
        this->initSimTime = std::clock();
        int i = 0;
        while (i < this->param.nSim)
        {
            this->EvolveAndRecordNodesOnce(seed, tTotal, i);
            this->WriteNodesData(i);
            this->ResetNodesAndEdges();
            this->ClearNetworkData();
            i++;
        }
        this->endSimTime = std::clock();
        std::string simTime = "# " + Misc::GetSimulationTime(this->initSimTime, this->endSimTime);
        std::cout << simTime << std::endl;
    }

    void OneSeedSpreading::PrepareToRecordNet(int tTotal)
    {
        int i = 0;
        this->param.tTotal = tTotal;
        double sum = 0.0;
        while (i < param.N)
        {
            sum += this->nodes[i]->GetState();
            i++;
        }
        if (this->param.writeOnRun)
        {
            this->outputFileNetDyn << this->timeStep << " " << sum << std::endl;
        }
        else
        {
            this->netDynData = Eigen::ArrayXd(tTotal + 1);
            this->netDynData(0) = sum;
        }
    }

    void OneSeedSpreading::TimeStepAndRecordNetInFile()
    {
        int i = 0, n = this->edges.size();
        while (i < n)
        {
            this->edges[i]->Step();
            i++;
        }
        i = 0;
        double  sum = 0.0;
        while (i < param.N)
        {
            this->nodes[i]->Step();
            sum += this->nodes[i]->GetState();
            i++;
        }
        this->outputFileNetDyn << this->timeStep << " " << sum << std::endl;
    }

    // time step of the simulation
    void OneSeedSpreading::TimeStepAndRecordNetInMemory()
    {
        int i = 0, n = this->edges.size();
        while (i < n)
        {
            this->edges[i]->Step();
            i++;
        }
        i = 0;
        double sum = 0.0;
        while (i < param.N)
        {
            this->nodes[i]->Step();
            sum += this->nodes[i]->GetState();
            i++;
        }
        this->netDynData(timeStep) = sum;
    }

    void OneSeedSpreading::EvolveAndRecordNetworkOnce(int seed, int tTotal, int simInd)
    {
        if (this->param.writeOnRun)
        {
            //this->outputFileNetDynName = "1s_net_nSim" + std::stringstream(simInd).str() + ".dat";
            this->outputFileNetDynName = this->GetOutputFileName(std::stringstream(simInd).str());
            this->outputFileNetDyn = Output::OpenFile(this->outputFileNetDynName, Input::GetInputParamString(param) + "# t    sum_of_nodes_signal");
        }

        this->SetInitialCondition(seed);

        int i = 0;
        this->timeStep = 0;
        this->PrepareToRecordNet(tTotal);
        //std::cout << nodesData.row(timeStep) << std::endl;
        //std::cin.get();
        this->timeStep = 1;
        while (this->timeStep <= tTotal)
        {
            (this->*TimeStepAndRecordNet)();
            //std::cout << nodesData.row(timeStep) << std::endl;
            //std::cin.get();
            this->timeStep++;
        }
        if (this->param.writeOnRun)
        {
            this->outputFileNetDyn.close();
            //Output::WriteFirstLine(this->outputFileNodesName, "# " + Input::inputCmdLine + "\n" + simTime);
        }
    }

    void OneSeedSpreading::EvolveAndRecordNetwork(int seed, int tTotal)
    {
        if (this->param.outputType == Output::OutputType::NetworkDynamics)
        {
            this->initSimTime = std::clock();
            int i = 0;
            while (i < this->param.nSim)
            {
                this->EvolveAndRecordNetworkOnce(seed, tTotal, i);
                this->WriteNetworkData(i);
                this->ResetNodesAndEdges();
                this->ClearNetworkData();
                i++;
            }
            this->endSimTime = std::clock();
            std::string simTime = "# " + Misc::GetSimulationTime(this->initSimTime, this->endSimTime);
            std::cout << simTime << std::endl;
        }
        else if (this->param.outputType == Output::OutputType::TotalNetworkSum)
        {
            if (this->param.writeOnRun)
            {
                this->outputFileTotNetSimName = this->GetOutputFileName("");
                this->outputFileTotNetSim = Output::OpenFile(this->outputFileTotNetSimName, Input::GetInputParamString(param) + "# nSim    sum_of_network_signal");
            }
            else
            {
                this->totalNetSimData = Eigen::ArrayXd(this->param.nSim);
            }

            this->initSimTime = std::clock();
            double s;
            int i = 0;
            while (i < this->param.nSim)
            {
                this->EvolveNetwork(seed, tTotal);
                s = this->SumNodesStates();
                if (this->param.writeOnRun)
                    this->outputFileTotNetSim << i << " " << s << std::endl;
                else
                    this->totalNetSimData(i) = s;
                i++;
            }
            this->endSimTime = std::clock();
            if (this->param.writeOnRun)
            {
                std::string simTime = "# " + Misc::GetSimulationTime(this->initSimTime, this->endSimTime);
                std::cout << simTime;
                this->outputFileTotNetSim.close();
                //Output::WriteFirstLine(this->outputFileNodesName, "# " + Input::inputCmdLine + "\n" + simTime);
            }
            else
            {
                this->WriteTotalNetworkData();
            }
        }
        else
        {
            throw std::invalid_argument("input.outputType not supported for OneSeedSpreading::EvolveAndRecordNetwork function");
        }
    }

    double OneSeedSpreading::SumNodesStates()
    {
        double sum = 0.0;
        int i = 0;
        while (i < this->param.N)
        {
            sum += this->nodes[i]->GetState();
            i++;
        }
        return sum;
    }

    void OneSeedSpreading::SetNodesThreshold(double theta)
    {
        int i = 0;
        while (i < this->param.N)
        {
            this->nodes[i]->SetThreshold(theta);
            i++;
        }
    }

    // time step
    void OneSeedSpreading::TimeStep()
    {
        int i = 0, n = this->edges.size();
        while (i < n)
        {
            this->edges[i]->Step();
            i++;
        }
        i = 0;
        this->isNetFullyActive = this->nodes[0]->GetState() > 0.0;
        while (i < this->param.N)
        {
            this->nodes[i]->Step();
            this->isNetFullyActive = this->isNetFullyActive && (this->nodes[i]->GetState() > 0.0);
            i++;
        }
    }

    void OneSeedSpreading::EvolveNetwork(int seed, int tTotal)
    {
        this->ResetNodesAndEdges();
        this->SetInitialCondition(seed);
        this->param.tTotal = tTotal;
        this->timeStep = 0;
        while (this->timeStep < tTotal)
        {
            this->TimeStep();
            this->timeStep++;
            if (this->isNetFullyActive)
                break;
        }
    }

    Misc::AvgVar_t OneSeedSpreading::RunSingleTheta(double theta)
    {
        this->SetNodesThreshold(theta);
        int i = 0;
        while (i < this->param.nSim)
        {
            this->EvolveNetwork(this->param.seedIndex, this->param.tTotal);
            this->totalNetSimDataPT[i] = this->SumNodesStates() / (double)this->param.N;
            i++;
        }
        return Misc::CalcAverageVariance(this->totalNetSimDataPT);
    }

    void OneSeedSpreading::PhaseTransition()
    {
        if (this->param.writeOnRun)
        {
            this->outputFilePhaseTransName = this->GetOutputFileName("");
            this->outputFilePhaseTrans = Output::OpenFile(this->outputFilePhaseTransName, Input::GetInputParamString(param) + "# theta    rho    Var(rho)");
        }
        else
        {
            this->thetaRhoChiData = Eigen::ArrayXXd(this->param.ntheta, 3);
        }

        this->initSimTime = std::clock();
        Eigen::ArrayXd theta = Eigen::ArrayXd::LinSpaced(this->param.ntheta, this->param.thetaMin, this->param.thetaMax);
        this->totalNetSimDataPT = std::vector<double>(this->param.nSim);
        int i = 0;
        while (i < this->param.ntheta)
        {
            Misc::AvgVar_t v = this->RunSingleTheta(theta(i));
            std::cout << "(" << i << "/"<< this->param.ntheta << ") theta = " << theta(i) << std::endl;
            if (this->param.writeOnRun)
                this->outputFilePhaseTrans << theta(i) << " " << v.Average << " " << v.Variance << std::endl;
            else
            {
                this->thetaRhoChiData(i, 0) = theta(i);
                this->thetaRhoChiData(i, 1) = v.Average;
                this->thetaRhoChiData(i, 2) = v.Variance;
            }
            i++;
        }
        this->endSimTime = std::clock();

        if (this->param.writeOnRun)
        {
            std::string simTime = "# " + Misc::GetSimulationTime(this->initSimTime, this->endSimTime);
            std::cout << simTime;
            this->outputFilePhaseTrans.close();
            Output::WriteFirstLine(this->outputFilePhaseTransName, "# " + Input::inputCmdLine + "\n" + simTime);
        }
        else
        {
            this->WritePhaseTransData();
        }
    }

    void OneSeedSpreading::SetInitialCondition(int seed)
    {
        if (this->param.simType == Input::SimulationType::FullyActive)
        {
            int i = 0;
            while (i < this->param.N)
            {
                this->nodes[i]->SetIC(1.0);
                i++;
            }
        }
        else
        {
            if (seed == -1) // random seed
            {
                seed = (int)std::floor((double)this->param.N * this->dist->NextNumber());
            }
            this->nodes[seed]->SetIC(1.0);
        }
    }

    // runs the simulation
    void OneSeedSpreading::SpreadAllToAllInFile()
    {
        this->outputFileSpAvgName = this->GetOutputFileName("avg");
        this->outputFileSpVarName = this->GetOutputFileName("var");
        this->outputFileSpAvg = Output::OpenFile(this->outputFileSpAvgName, Input::GetInputParamString(param) + "# seed(lines) Vs.target(columns) : time to spread from one to the other");
        this->outputFileSpVar = Output::OpenFile(this->outputFileSpVarName, Input::GetInputParamString(param) + "# seed(lines) Vs.target(columns) : variance of time to spread from one to the other");
        initSimTime = std::clock();
        Misc::AvgVar_t v;
        int i = 0, j;
        while (i < param.N)
        {
            j = 0;
            while (j < param.N)
            {
                v = SpreadNodeToNodeManyTimes(i, j);
                this->outputFileSpAvg << " " << v.Average;
                this->outputFileSpVar << " " << v.Variance;
                j++;
            }
            this->outputFileSpAvg << std::endl;
            this->outputFileSpVar << std::endl;
            i++;
        }
        endSimTime = std::clock();
        std::string simTime = "# " + Misc::GetSimulationTime(this->initSimTime, this->endSimTime);
        std::cout << simTime;
        this->outputFileSpAvg.close();
        this->outputFileSpVar.close();
        Output::WriteFirstLine(this->outputFileSpAvgName, "# " + Input::inputCmdLine + "\n" + simTime);
        Output::WriteFirstLine(this->outputFileSpVarName, "# " + Input::inputCmdLine + "\n" + simTime);
    }

    // runs the simulation
    void OneSeedSpreading::SpreadAllToAllInMemory()
    {
        initSimTime = std::clock();
        spreadingTimeAvg = Eigen::ArrayXXd(param.N, param.N);
        spreadingTimeVar = Eigen::ArrayXXd(param.N, param.N);
        Misc::AvgVar_t v;
        int i = 0, j;
        while (i < param.N)
        {
            j = 0;
            while (j < param.N)
            {
                v = SpreadNodeToNodeManyTimes(i, j);
                spreadingTimeAvg(i, j) = v.Average;
                spreadingTimeVar(i, j) = v.Variance;
                j++;
            }
            i++;
        }
        endSimTime = std::clock();
    }

    void OneSeedSpreading::SpreadAllToAll()
    {
        (this->*SpreadAllToAll_int)();
        this->WriteSpreadingTimeMatrix();
    }

    void OneSeedSpreading::ResetNodesAndEdges()
    {
        int i = 0, n = edges.size();
        while (i < n)
        {
            edges[i]->Reset();
            i++;
        }
        i = 0;
        while (i < param.N)
        {
            nodes[i]->SetIC(0.0);
            i++;
        }
        this->timeStep = 0;
    }

    void OneSeedSpreading::ClearNetworkData()
    {
        //nodes = std::vector<std::unique_ptr<Elements::IElement>>();
        //edges = std::vector<std::unique_ptr<Couplings::ICoupling>>();
        spreadingTimeAvg = Eigen::ArrayXXd();
        spreadingTimeVar = Eigen::ArrayXXd();
        nodesData = Eigen::ArrayXXd();
        netDynData = Eigen::ArrayXd();

        /*int i = 0, n = nodes.size();
        while (i < n)
        {
            delete nodes[i];
            i++;
        }
        i = 0; n = edges.size();
        while (i < n)
        {
            delete edges[i];
            i++;
        }//*/
    }

    // seedNode, targetNode // spreads activity from seedNode to targetNode
    int OneSeedSpreading::SpreadNodeToNode(int seedNode, int targetNode)
    {
        ResetNodesAndEdges();
        nodes[seedNode]->SetIC(1.0);
        timeStep = 0;
        while (timeStep < param.tTotal)
        {
            if (nodes[targetNode]->GetState() > 0.0)
            {
                return this->timeStep;
            }

            this->TimeStep();

            timeStep++;
        }
        return param.tTotal;
    }

    // seedNode, targetNode // spreads activity from seedNode to targetNode nSim times
    Misc::AvgVar_t OneSeedSpreading::SpreadNodeToNodeManyTimes(int seedNode, int targetNode)
    {
        if (seedNode == targetNode)
        {
            Misc::AvgVar_t v;
            v.Average = 0;
            v.Variance = 0;
            return v;
        }
        std::vector<double> times(param.nSim);
        int k = 0;
        while (k < param.nSim)
        {
            times[k] = (double)SpreadNodeToNode(seedNode, targetNode);
            k++;
        }
        return Misc::CalcAverageVariance(times);
    }

    std::string OneSeedSpreading::GetOutputFileName(std::string middle = "")
    {
        std::string s = "";
        if (this->param.simType == Input::SimulationType::OneActive)
            s += "1s";
        else
            s += "s";
        if (this->param.outputType == Output::OutputType::NodesDynamics)
            s += "_dyn";
        else if (this->param.outputType == Output::OutputType::TotalNetworkSum)
            s += "_ava";
        else if (this->param.outputType == Output::OutputType::NetworkDynamics)
            s += "_net";
        else if (this->param.outputType == Output::OutputType::PhaseTransition)
            s += "_bif";
        else if (this->param.outputType == Output::OutputType::SpreadingMatrix)
            s += "_spr";
        if (!middle.empty())
            s += "_" + middle;
        if (!this->param.outFileSuf.empty())
            s += "_" + this->param.outFileSuf;
        s += ".dat";
        return s;
    }

    void OneSeedSpreading::WritePhaseTransData()
    {
        if (!this->param.writeOnRun)
        {
            std::string headerTxt = "# " + Input::inputCmdLine + "\n" + Input::GetInputParamString(param) + "# " + Misc::GetSimulationTime(initSimTime, endSimTime) + "\n# theta    rho    Var(rho)";
            Output::WriteData(this->thetaRhoChiData, headerTxt, this->GetOutputFileName(""));//*/
        }
    }

    // writes spreadingTime matrix
    void OneSeedSpreading::WriteSpreadingTimeMatrix()
    {
        if (!this->param.writeOnRun)
        {
            std::string headerTxt = "# " + Input::inputCmdLine + "\n" + Input::GetInputParamString(param) + "# " + Misc::GetSimulationTime(initSimTime, endSimTime) + "\n# seed (lines) Vs. target (columns): time to spread from one to the other";
            Output::WriteData(spreadingTimeAvg, headerTxt, this->GetOutputFileName("avg"));
            headerTxt = "# " + Input::inputCmdLine + "\n" + Input::GetInputParamString(param) + "# " + Misc::GetSimulationTime(initSimTime, endSimTime) + "\n# seed (lines) Vs. target (columns): variance of time to spread from one to the other";
            Output::WriteData(spreadingTimeVar, headerTxt, this->GetOutputFileName("var"));
        }
    }

    void OneSeedSpreading::WriteTotalNetworkData()
    {
        if (!this->param.writeOnRun)
        {
            std::string headerTxt = "# " + Input::inputCmdLine + "\n" + Input::GetInputParamString(param) + "# " + Misc::GetSimulationTime(initSimTime, endSimTime) + "\n# nSim    sum_of_network_activity";
            Output::WriteData(this->totalNetSimData, headerTxt, this->GetOutputFileName(""), true);//*/
        }
    }

    void OneSeedSpreading::WriteNetworkData(int i = 0)
    {
        if (!this->param.writeOnRun)
        {
            std::string headerTxt = "# " + Input::inputCmdLine + "\n" + Input::GetInputParamString(param) + "# " + Misc::GetSimulationTime(initSimTime, endSimTime) + "\n# t    sum_of_nodes_activity";
            Output::WriteData(this->netDynData, headerTxt, this->GetOutputFileName(std::stringstream(i).str()), true);//*/
        }
    }

    void OneSeedSpreading::WriteNodesData(int i = 0)
    {
        if (!this->param.writeOnRun)
        {
            std::string headerTxt = "# " + Input::inputCmdLine + "\n" + Input::GetInputParamString(param) + "# " + Misc::GetSimulationTime(initSimTime, endSimTime) + "\n# t    nodes";
            Output::WriteData(this->nodesData, headerTxt, this->GetOutputFileName(std::stringstream(i).str()), true);//*/
        }
    }
}
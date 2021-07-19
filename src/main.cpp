#define _CRT_SECURE_NO_WARNINGS
#include <memory>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector>
#include "Eigen/Core"
#include "Eigen/SparseCore" // Eigen; added by right click on project -> Properties -> C/C++ -> Additional Include Directories
#include "Elements.h"
#include "Couplings.h"
#include "Networks.h"
#include "Misc.h"
#include "Input.h"
#include "Output.h"
#include "OneSeedSpreading.h"

void TestSIElements();
void TestNeighborSIElements();
void TestSIRDetElement();
void TestNeighborSIRDetElement();
void TestSparseMatrix();
void TestAdjMatFromFile();
void TestAdjMatrices(std::mt19937_64&);
void TestEnumFunc();//*/
void TestOutputFile();
void TestRunProcess();
void TestSimpleNetworkNetData(std::mt19937_64&, bool, Output::OutputType);
void TestSimpleNetwork(std::mt19937_64&, bool);
void TestSpreadingOneSeed(std::mt19937_64&, bool);
void TestPhaseTransition(std::mt19937_64&, bool);

int main(int argCount, char* args[])
{
    std::mt19937_64 rand;
    unsigned rSeed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
    rand.seed(rSeed);

    /*TestOutputFile();//*/
    /*TestSpreadingOneSeed(rand, true);
    TestSimpleNetwork(rand, true);
    TestSimpleNetworkNetData(rand, true, Output::OutputType::TotalNetworkSum);
    TestSimpleNetworkNetData(rand, true, Output::OutputType::NetworkDynamics);
    TestPhaseTransition(rand, true);//*/
    /*TestEnumFunc();
    TestAdjMatrices(rand);
    TestRunProcess();
    TestAdjMatFromFile();
    TestSparseMatrix();//*/
    /*TestSIElements();
    TestNeighborSIElements();
    TestSIRDetElement();
    TestNeighborSIRDetElement();//*/

    Input::InputParam_t input;

    #ifndef _DEBUG
    try
    {
    #endif
        input = Input::GetInputArgs(argCount, args);
    #ifndef _DEBUG
    }
    catch (const std::exception& e)
    {
        std::cout << "Input ERROR: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    #endif

    #ifndef _DEBUG
    try
    {
    #endif
        if ((input.simType == Input::SimulationType::OneActive) || (input.simType == Input::SimulationType::FullyActive))
        {
            Simulations::OneSeedSpreading sim(input, rand);
            if (input.outputType == Output::OutputType::SpreadingMatrix)
            {
                sim.SpreadAllToAll();
            }
            else if (input.outputType == Output::OutputType::NodesDynamics)
            {
                sim.EvolveAndRecordNodes(input.seedIndex, input.tTotal);
            }
            else if ((input.outputType == Output::OutputType::NetworkDynamics) || (input.outputType == Output::OutputType::TotalNetworkSum))
            {
                sim.EvolveAndRecordNetwork(input.seedIndex, input.tTotal);
            }
            else if (input.outputType == Output::OutputType::PhaseTransition)
            {
                sim.PhaseTransition();
            }
            exit(EXIT_SUCCESS);
        }
        else
        {
            throw std::invalid_argument("SimulationType not implemented");
        }
    #ifndef _DEBUG
    }
    catch (const std::exception& e)
    {
        std::cout << "Simulation ERROR: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    #endif//*/

#ifdef _DEBUG
    std::cin.get();
#endif
}

void TestPhaseTransition(std::mt19937_64& rand, bool wOnRun)
{
    Input::InputParam_t input = Input::GetDefaultValues();
    input.Lx = 3;
    input.Ly = 3;
    input.Lz = 1;
    input.N = 9;
    input.netType = Networks::AdjacencyMatrixType::SquareLattice;
    input.weightType = Networks::AdjacencyMatrixWeightType::Homogeneous;
    input.outputType = Output::OutputType::PhaseTransition;
    input.weightAmp = 0.5;
    input.isPeriodic = false;
    input.theta = 0.1;
    input.writeOnRun = wOnRun;
    input.nSim = 5;
    input.ntheta = 5;
    input.thetaMin = 0.01;
    input.thetaMax = 1.01;
    input.tTotal = 1000;
    input.seedIndex = -1;
    input.simType = Input::SimulationType::OneActive;
    Simulations::OneSeedSpreading sim(input, rand);
    sim.PhaseTransition();
}

void TestSpreadingOneSeed(std::mt19937_64& rand, bool wOnRun)
{
    Input::InputParam_t input = Input::GetDefaultValues();
    input.Lx = 3;
    input.Ly = 3;
    input.Lz = 1;
    input.N = 9;
    input.netType = Networks::AdjacencyMatrixType::SquareLattice;
    input.weightType = Networks::AdjacencyMatrixWeightType::Homogeneous;
    input.outputType = Output::OutputType::SpreadingMatrix;
    input.weightAmp = 0.5;
    input.isPeriodic = false;
    input.theta = 0.1;
    input.writeOnRun = wOnRun;
    Simulations::OneSeedSpreading sim(input, rand);
    sim.SpreadAllToAll();
}

void TestSimpleNetwork(std::mt19937_64& rand, bool wOnRun)
{
    Input::InputParam_t input = Input::GetDefaultValues();
    input.Lx = 3;
    input.Ly = 3;
    input.Lz = 1;
    input.N = 9;
    input.netType = Networks::AdjacencyMatrixType::SquareLattice;
    input.weightType = Networks::AdjacencyMatrixWeightType::Homogeneous;
    input.outputType = Output::OutputType::NodesDynamics;
    input.weightAmp = 0.5;
    input.isPeriodic = false;
    input.theta = 0.1;
    input.writeOnRun = wOnRun;
    input.nSim = 2;
    Simulations::OneSeedSpreading sim(input, rand);
    sim.EvolveAndRecordNodes(4, 10);
}

void TestSimpleNetworkNetData(std::mt19937_64& rand, bool wOnRun, Output::OutputType ot)
{
    Input::InputParam_t input = Input::GetDefaultValues();
    input.Lx = 3;
    input.Ly = 3;
    input.Lz = 1;
    input.N = 9;
    input.netType = Networks::AdjacencyMatrixType::SquareLattice;
    input.weightType = Networks::AdjacencyMatrixWeightType::Homogeneous;
    input.outputType = ot;
    input.weightAmp = 0.5;
    input.isPeriodic = false;
    input.theta = 0.1;
    input.writeOnRun = wOnRun;
    input.nSim = 2;
    Simulations::OneSeedSpreading sim(input, rand);
    sim.EvolveAndRecordNetwork(4, 10);
}

void TestEnumFunc()
{
    Input::SimulationType s = Input::SimulationType::nActive;
    std::cout << "Possible values: " << std::endl;
    std::cout << Misc::EnumStrList<Input::SimulationType>(Input::SimulationTypeChar) << std::endl;
    std::cout << "Standard assigned value = " << Misc::EnumToStr(s, Input::SimulationTypeChar) << std::endl;
    s = Misc::StrToEnum<Input::SimulationType>("FullyActive", Input::SimulationTypeChar);
    std::cout << "String assigned value = " << Misc::EnumToStr(s, Input::SimulationTypeChar) << std::endl;
    std::cin.get();
}

void TestRunProcess()
{
    std::string op;
    op = Misc::RunProcess("cmd /c dir");
    std::cout << op << std::endl;
    std::cin.get();
}

void TestAdjMatFromFile()
{
    Networks::AdjacencyMatrixFromFile* A = new Networks::AdjacencyMatrixFromFile("matrix.txt", false);
    Eigen::SparseMatrix<double> M;
    M = A->BuildAndGetMatrix();
    std::cout << M.toDense() << std::endl;
    std::cin.get();
    delete A;
}

void TestSparseMatrix()
{
    Eigen::SparseMatrix<double> M(5, 5);
    std::vector<Eigen::Triplet<double>> el;
    el.reserve(8);
    el.push_back(Eigen::Triplet<double>(0, 0, 1.0));
    el.push_back(Eigen::Triplet<double>(1, 3, 2.0));
    el.push_back(Eigen::Triplet<double>(2, 4, 2.0));
    el.push_back(Eigen::Triplet<double>(2, 2, 3.0));
    el.push_back(Eigen::Triplet<double>(3, 1, 4.0));
    el.push_back(Eigen::Triplet<double>(4, 0, 5.0));
    el.push_back(Eigen::Triplet<double>(3, 2, 4.0));
    el.push_back(Eigen::Triplet<double>(4, 4, 5.0));
    M.setFromTriplets(el.begin(), el.end());
    std::cout << M.toDense() << std::endl;
    std::cin.get();
    std::cout << std::endl;
    std::cout << M.triangularView<Eigen::Upper>().toDense() << std::endl;
    std::cin.get();
    std::cout << std::endl;
    std::cout << M.triangularView<Eigen::Lower>().toDense() << std::endl;
    std::cin.get();
    std::cout << std::endl;
    std::cout << M.triangularView<Eigen::StrictlyUpper>().toDense() << std::endl;
    std::cin.get();
    std::cout << std::endl;
    std::cout << M.triangularView<Eigen::StrictlyLower>().toDense() << std::endl;
    std::cin.get();
}

void TestSIRDetElement()
{
    std::vector<double> tt, x1, x2;
    std::vector<std::unique_ptr<Elements::IElement>> nodes;
    nodes.push_back(std::make_unique<Elements::SIRSDetElement>(0, 0.5, 0.0, 3, false));
    x1.push_back(dynamic_cast<Elements::SIRSDetElement&>(*nodes[0]).GetStateId());
    x2.push_back(nodes[0]->GetState());
    int t = 0;
    tt.push_back(t);
    while (t < 100)
    {
        nodes[0]->Step((t == 20 ? 1.0 : 0.0));
        t++;
        x1.push_back(dynamic_cast<Elements::SIRSDetElement&>(*nodes[0]).GetStateId());
        x2.push_back(nodes[0]->GetState());
        tt.push_back(t);
    }
    std::string headerTxt = "# t   xSID   x";
    std::string fileName = "testNodes03.dat";
    Output::WriteData(tt, x1, x2, headerTxt, fileName);
    //Output::WriteData(tt, x1, x2, static_cast<std::string>("# t   xSID   x"), static_cast<std::string>("testNodes03.dat"));
    //delete nodes[0];
}

void TestNeighborSIRDetElement()
{
    std::vector<double> tt, x1, x2;
    std::vector<std::unique_ptr<Elements::IElement>> nodes;
    std::vector<std::unique_ptr<Couplings::ICoupling>> edges;
    nodes.push_back(std::make_unique<Elements::SIRSDetElement>(0, 0.5, 0.0, 3, false));
    nodes.push_back(std::make_unique<Elements::SIRSDetElement>(1, 0.5, 0.0, 3, false));
    edges.push_back(std::make_unique<Couplings::PulseCoupling>(*nodes[0], *nodes[1], 1.0));
    nodes[1]->AddNeighbor(*edges[0]);
    x1.push_back(nodes[0]->GetState());
    x2.push_back(nodes[1]->GetState());
    int t = 0;
    tt.push_back(t);
    while (t < 100)
    {
        edges[0]->Step();
        nodes[0]->Step((t == 20 ? 1.0 : 0.0));
        nodes[1]->Step();
        t++;
        x1.push_back(nodes[0]->GetState());
        x2.push_back(nodes[1]->GetState());
        tt.push_back(t);
    }
    std::string headerTxt = "# t   x1   x2";
    std::string fileName = "testNodes04.dat";
    Output::WriteData(tt, x1, x2, headerTxt, fileName);
    //Output::WriteData(tt, x1, x2, static_cast<std::string>("# t   x1   x2"), static_cast<std::string>("testNodes04.dat"));
    //delete nodes[0], nodes[1], edges[0];
}

void TestSIElements()
{
    std::vector<double> tt, x1, x2;
    std::vector<std::unique_ptr<Elements::IElement>> nodes;
    nodes.push_back(std::make_unique<Elements::SIElement>(0, 0.5, 0.0, false));
    nodes.push_back(std::make_unique<Elements::SIElement>(1, 0.5, 1.0, false));
    x1.push_back(nodes[0]->GetState());
    x2.push_back(nodes[1]->GetState());
    int t = 0;
    tt.push_back(t);
    while (t < 100)
    {
        nodes[0]->Step((t == 20 ? 1.0 : 0.0));
        nodes[1]->Step((t == 20 ? 1.0 : 0.0));
        t++;
        x1.push_back(nodes[0]->GetState());
        x2.push_back(nodes[1]->GetState());
        tt.push_back(t);
    }
    std::string headerTxt = "# t   x1   x2";
    std::string fileName = "testNodes01.dat";
    Output::WriteData(tt, x1, x2, headerTxt, fileName);
    //Output::WriteData(tt, x1, x2, static_cast<std::string>("# t   x1   x2"), static_cast<std::string>("testNodes01.dat"));
    //delete nodes[0], nodes[1];
}

void TestNeighborSIElements()
{
    std::vector<double> tt, x1, x2;
    std::vector<std::unique_ptr<Elements::IElement>> nodes;
    std::vector<std::unique_ptr<Couplings::ICoupling>> edges;
    nodes.push_back(std::make_unique<Elements::SIElement>(0, 0.5, 0.0, false));
    nodes.push_back(std::make_unique<Elements::SIElement>(1, 0.5, 0.0, false));
    edges.push_back(std::make_unique<Couplings::PulseCoupling>(*nodes[0], *nodes[1], 1.0));
    nodes[1]->AddNeighbor(*edges[0]);
    x1.push_back(nodes[0]->GetState());
    x2.push_back(nodes[1]->GetState());
    int t = 0;
    tt.push_back(t);
    while (t < 100)
    {
        edges[0]->Step();
        nodes[0]->Step((t == 20 ? 1.0 : 0.0));
        nodes[1]->Step();
        t++;
        x1.push_back(nodes[0]->GetState());
        x2.push_back(nodes[1]->GetState());
        tt.push_back(t);
    }
    std::string headerTxt = "# t   x1   x2";
    std::string fileName = "testNodes02.dat";
    Output::WriteData(tt, x1, x2, headerTxt, fileName);
    //Output::WriteData(tt, x1, x2, static_cast<std::string>("# t   x1   x2"), static_cast<std::string>("testNodes02.dat"));
    //delete nodes[0], nodes[1], edges[0];
}

void TestAdjMatrices(std::mt19937_64& rand)
{
    std::ofstream fs;

    Networks::AdjacencyMatrixParam_t par;
    par.type = Networks::AdjacencyMatrixType::LinearLattice;
    par.nElems = 10;
    par.L.push_back(10);
    par.L.push_back(1);
    par.L.push_back(1);
    par.numOfNeighbours = 2;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("linearFree.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::LinearLattice;
    par.nElems = 10;
    par.L.push_back(10);
    par.L.push_back(1);
    par.L.push_back(1);
    par.numOfNeighbours = 2;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = true;
    par.netFileName = "";
    fs.open("linearPeriodic.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::SquareLattice;
    par.nElems = 16;
    par.L.push_back(4);
    par.L.push_back(4);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("squareFree.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::SquareLattice;
    par.nElems = 16;
    par.L.push_back(4);
    par.L.push_back(4);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = true;
    par.netFileName = "";
    fs.open("squarePeriodic.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::SquareLattice;
    par.nElems = 15;
    par.L.push_back(3);
    par.L.push_back(5);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("squareIFree.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::SquareLattice;
    par.nElems = 15;
    par.L.push_back(3);
    par.L.push_back(5);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = true;
    par.netFileName = "";
    fs.open("squareIPeriodic.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::BarabasiAlbertGraph;
    par.nElems = 15;
    par.L.push_back(15);
    par.L.push_back(1);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("ba.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::WattsStrogatzGraph;
    par.nElems = 15;
    par.L.push_back(15);
    par.L.push_back(1);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("ws.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::ConnectedWattsStrogatzGraph;
    par.nElems = 15;
    par.L.push_back(15);
    par.L.push_back(1);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("wsc.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::CubicLattice;
    par.nElems = 27;
    par.L.push_back(3);
    par.L.push_back(3);
    par.L.push_back(3);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("cubicFree.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::CubicLattice;
    par.nElems = 27;
    par.L.push_back(3);
    par.L.push_back(3);
    par.L.push_back(3);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("cubicPeriodic.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::RandomGraph;
    par.nElems = 15;
    par.L.push_back(15);
    par.L.push_back(1);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("random.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();

    par.L.clear();
    par.type = Networks::AdjacencyMatrixType::CompleteGraph;
    par.nElems = 15;
    par.L.push_back(15);
    par.L.push_back(1);
    par.L.push_back(1);
    par.numOfNeighbours = 4;
    par.numOfEdgesForNewElem = 3;
    par.rewiringProb = 0.02;
    par.isDirected = false;
    par.isPeriodic = false;
    par.netFileName = "";
    fs.open("complete.txt", std::ios::trunc);
    fs << Networks::GetAdjacencyMatrix(par, rand)->BuildAndGetMatrix().toDense();
    fs.close();
}

void TestOutputFile()
{
    //std::fstream fs = Output::OpenFile(static_cast<std::string>("test.txt"), "teste linha 1");
    std::string fileName = "test.txt";
    std::fstream fs = Output::OpenFile(fileName, "teste linha 1");
    fs.close();
    fs = std::fstream("test.txt", std::fstream::out);
    fs << "teste linha 2";
    fs.close();
}
//*/

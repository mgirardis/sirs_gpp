#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdexcept>
#include <math.h>
#include "Input.h"
#include <iostream>
#include <sstream> 
#include <fstream>
#include <sys/stat.h>
#include "Misc.h"
#include "Elements.h"
#include "Couplings.h"
#include "Networks.h"
#include "Output.h"
//#include "SIRSSim.h"
//#include "Distribution_c.h"

namespace Input
{
    /***
     ************************************************
     ***
     *** HOW TO ADD INPUT PARAMETERS TO THE PROGRAM
     ***
     ************************************************
     ***
     *** 1) create a string variable par_NAME = "NAME" (NAME is the name of the input command line parameter) in the Input.h file
     ***
     *** 2) add NAME as a field of the InputParam_t struct
     ***
     *** 3) add its default value to the function Input::GetDefaultValues
     ***
     *** 4) add an if statement to the function Input::GetInputArgs in order to obtain the parameter value from the input command line
     ***
     *** 5) add the parameter to the function Input::GetInputParamString
     ***
     *** 6) add parameter description to the Input::PrintHelp function
     ***
     *** 7) if boolean, add an OR condition to the function Input::IsBooleanParameter
     ***
     *** 8) add validation to the function Input::ValidateInputArgs
     ***
     ***/
    std::string inputCmdLine = "";

    void PrintMsg(std::string message)
    {
        std::cout << message << std::endl;
    }

    void PrintError(std::string message)
    {
        PrintMsg("#... ERROR: " + message);
    }

    void PrintHelp(char* args[])
    {
        InputParam_t input = GetDefaultValues();
        std::cout << "#..." << std::endl;
        std::cout << "#... Simulates a Susceptible-Infected(-Refratory(-Susceptible)) model connected with the desired coupling in the desired network" << std::endl;
        std::cout << "#... according to Misic et al (2015), Spreading process. Neuron." << std::endl;
        std::cout << "#..." << std::endl;
        std::cout << "#... Usage: " << args[0] << " netType NET_TYPE theta DOUBLE nSim INT ..." << std::endl;
        std::cout << "#... all args are optional" << std::endl;
        std::cout << "#..." << std::endl;
        std::cout << "#... where:" << std::endl;
        std::cout << "#..." << std::endl;
        std::cout << "#... " << par_simType.c_str() << " (" << Misc::EnumToStr(input.simType, SimulationTypeChar).c_str() << ")\t\t\t propagates activity starting from one, two, n or no seed node" << std::endl;
        std::cout << "#... \t\t Other values: " << Misc::EnumStrList<SimulationType>(SimulationTypeChar) << std::endl;
        std::cout << "#... " << par_seedType.c_str() << " (" << Misc::EnumToStr(input.seedType, SeedTypeChar).c_str() << ")\t\t\t defines the behavior of the superposition of many seeds activity" << std::endl;
        std::cout << "#... \t\t Other values: " << Misc::EnumStrList<SeedType>(SeedTypeChar) << std::endl;
        std::cout << "#... " << par_seedIndex.c_str() << " (" << input.seedIndex << ")\t\t\t index of seed (for OutputType == NetworkDynamics or NodeDynamics)" << std::endl;
        std::cout << "#... " << par_nSim.c_str() << " (" << input.nSim << ")\t\t\t number of realizations for each seed" << std::endl;
        std::cout << "#... " << par_nodeType.c_str() << " (" << Misc::EnumToStr(input.nodeType, Elements::ElementTypeChar).c_str() << ")\t\t\t selects node type" << std::endl;
        std::cout << "#... \t\t Other values: " << Misc::EnumStrList<Elements::ElementType>(Elements::ElementTypeChar) << std::endl;
        std::cout << "#... " << par_theta.c_str() << " (" << input.theta << ")\t\t\t threshold for activation of a node (homogeneous)" << std::endl;
        std::cout << "#... " << par_nRefStates.c_str() << " (" << input.nRefStates << ")\t\t\t number of refractory states in case of SIR(S) element" << std::endl;
        std::cout << "#... " << par_coupType.c_str() << " (" << Misc::EnumToStr(input.coupType, Couplings::CouplingTypeChar).c_str() << ")\t\t\t selects coupling type" << std::endl;
        std::cout << "#... \t\t Other values: " << Misc::EnumStrList<Couplings::CouplingType>(Couplings::CouplingTypeChar) << std::endl;
        std::cout << "#... " << par_netType.c_str() << " (" << Misc::EnumToStr(input.netType, Networks::AdjacencyMatrixTypeChar).c_str() << ")\t\t\t selects network type" << std::endl;
        std::cout << "#... \t\t Other values: " << Misc::EnumStrList<Networks::AdjacencyMatrixType>(Networks::AdjacencyMatrixTypeChar) << std::endl;
        std::cout << "#... " << par_N.c_str() << " (" << input.N << ")\t\t\t number of nodes in the network" << std::endl;
        std::cout << "#... " << par_Lx.c_str() << " (" << input.Lx << ")\t\t\t number of nodes in x-direction for lattices" << std::endl;
        std::cout << "#... " << par_Ly.c_str() << " (" << input.Ly << ")\t\t\t number of nodes in y-direction for lattices" << std::endl;
        std::cout << "#... " << par_Lz.c_str() << " (" << input.Lz << ")\t\t\t number of nodes in z-direction for lattices" << std::endl;
        std::cout << "#... " << par_numOfNeighWS.c_str() << " (" << input.numOfNeighWS << ")\t\t\t number of neighbors for Watts-Strogatz network" << std::endl;
        std::cout << "#... " << par_numOfEdgesBA.c_str() << " (" << input.numOfEdgesBA << ")\t\t\t number of new edges for Barabasi-Albert network" << std::endl;
        std::cout << "#... " << par_rewiringProb.c_str() << " (" << input.rewiringProb << ")\t\t\t rewiring probability for WS net and wiring prob for Random net" << std::endl;
        std::cout << "#... " << par_isDirected.c_str() << " (" << input.isDirected << ")\t\t\t determines if the network is directed" << std::endl;
        std::cout << "#... " << par_isPeriodic.c_str() << " (" << input.isPeriodic << ")\t\t\t determines if the network is periodic" << std::endl;
        std::cout << "#... " << par_netFile.c_str() << " (" << input.netFile << ")\t\t\t input file containing adjacency matrix (weighted or not) if netType == FromFile" << std::endl;
        std::cout << "#... " << par_weightType.c_str() << " (" << Misc::EnumToStr(input.weightType, Networks::AdjacencyMatrixWeightTypeChar).c_str() << ")\t\t\t selects whether adjacency matrix is weighted" << std::endl;
        std::cout << "#... \t\t Other values: " << Misc::EnumStrList<Networks::AdjacencyMatrixWeightType>(Networks::AdjacencyMatrixWeightTypeChar) << std::endl;
        std::cout << "#... " << par_weightAmp.c_str() << " (" << input.weightAmp << ")\t\t\t if UniformRandom weightType, then weights are distributed in [0,weightAmp]; if Homogeneous, all weights are weightAmp" << std::endl;
        std::cout << "#... " << par_tTotal.c_str() << " (" << input.tTotal << ")\t\t\t total time of simulation (if required)" << std::endl;
        std::cout << "#... " << par_writeOnRun.c_str() << " (" << input.writeOnRun << ")\t\t\t if true, writes output files on the fly; if false, writes only after all the calculations" << std::endl;
        std::cout << "#... " << par_outputType.c_str() << " (" << Misc::EnumToStr(input.outputType, Output::OutputTypeChar).c_str() << ")\t\t\t selects whether adjacency matrix is weighted" << std::endl;
        std::cout << "#... \t\t Other values: " << Misc::EnumStrList<Output::OutputType>(Output::OutputTypeChar) << std::endl;
        std::cout << "#... " << par_thetaMin.c_str() << " (" << input.thetaMin << ")\t\t\t minimum element threshold for phase transition simulation" << std::endl;
        std::cout << "#... " << par_thetaMax.c_str() << " (" << input.thetaMax << ")\t\t\t maximum element threshold for phase transition simulation" << std::endl;
        std::cout << "#... " << par_ntheta.c_str() << " (" << input.ntheta << ")\t\t\t number of thresholds to simulate, inclusive in [thetaMin;thetaMax]" << std::endl;
        std::cout << "#... " << par_outFileSuf.c_str() << " (" << input.outFileSuf << ")\t\t\t output filename suffix" << std::endl;
        std::cout << "#..." << std::endl;
        std::cout << "#... h, help, -h = prints help" << std::endl;
        std::cout << "#..." << std::endl;
        return;
    }//*/

    bool IsBooleanParameter(std::string& par)
    {
        return (par == par_isDirected) || (par == par_isPeriodic) || (par == par_writeOnRun);
    }

    InputParam_t GetDefaultValues()
    {
        InputParam_t input;// (std::mt19937_64());
        input.simType = SimulationType::OneActive;
        input.seedType = SeedType::Cooperative;
        input.nSim = 100;
        input.nodeType = Elements::ElementType::SI;
        input.theta = 0.01;
        input.nRefStates = 5;
        input.coupType = Couplings::CouplingType::PulseCoupling;
        input.netType = Networks::AdjacencyMatrixType::SquareLattice;
        input.N = 400;
        input.Lx = 20;
        input.Ly = 20;
        input.Lz = 1;
        input.numOfNeighWS = 4;
        input.numOfEdgesBA = 3;
        input.rewiringProb = 0.02;
        input.isDirected = false;
        input.isPeriodic = false;
        input.netFile = "";
        input.weightType = Networks::AdjacencyMatrixWeightType::Binary;
        input.weightAmp = 1.0;
        input.tTotal = std::numeric_limits<int>::max();
        input.writeOnRun = false;
        input.outputType = Output::OutputType::SpreadingMatrix;
        input.seedIndex = -1;
        input.thetaMin = 0.0;
        input.thetaMax = 1.0;
        input.ntheta = 100;
        input.outFileSuf = "";
        return input;
    }

    std::string GetInputParamString(InputParam_t input)
    {
        std::stringstream ss;
        ss << "# simType = " << Misc::EnumToStr(input.simType, SimulationTypeChar) << std::endl;
        ss << "# seedType = " << Misc::EnumToStr(input.seedType, SeedTypeChar) << std::endl;
        ss << "# seedIndex = " << input.seedIndex << std::endl;
        ss << "# nSim = " << input.nSim << std::endl;
        ss << "# nodeType = " << Misc::EnumToStr(input.nodeType, Elements::ElementTypeChar) << std::endl;
        ss << "# theta = " << input.theta << std::endl;
        ss << "# nRefStates = " << input.nRefStates << std::endl;
        ss << "# coupType = " << Misc::EnumToStr(input.coupType, Couplings::CouplingTypeChar) << std::endl;
        ss << "# netType = " << Misc::EnumToStr(input.netType, Networks::AdjacencyMatrixTypeChar) << std::endl;
        ss << "# N = " << input.N << std::endl;
        ss << "# Lx = " << input.Lx << std::endl;
        ss << "# Ly = " << input.Ly << std::endl;
        ss << "# Lz = " << input.Lz << std::endl;
        ss << "# numOfNeighWS = " << input.numOfNeighWS << std::endl;
        ss << "# numOfEdgesBA = " << input.numOfEdgesBA << std::endl;
        ss << "# rewiringProb = " << input.rewiringProb << std::endl;
        ss << "# isDirected = " << input.isDirected << std::endl;
        ss << "# isPeriodic = " << input.isPeriodic << std::endl;
        ss << "# netFile = " << input.netFile << std::endl;
        ss << "# weightType = " << Misc::EnumToStr(input.weightType, Networks::AdjacencyMatrixWeightTypeChar) << std::endl;
        ss << "# weightAmp = " << input.weightAmp << std::endl;
        ss << "# tTotal = " << input.tTotal << std::endl;
        ss << "# writeOnRun = " << input.writeOnRun << std::endl;
        ss << "# outputType = " << Misc::EnumToStr(input.outputType, Output::OutputTypeChar) << std::endl;
        ss << "# thetaMin = " << input.thetaMin << std::endl;
        ss << "# thetaMax = " << input.thetaMax << std::endl;
        ss << "# ntheta = " << input.ntheta << std::endl;
        ss << "# outFileSuf = " << input.outFileSuf << std::endl;
        return ss.str();
    }

    bool ValidateInputArgs(InputParam_t input)
    {
        if (input.nSim <= 0)
        {
            throw std::invalid_argument("nSim > 0");
        }
        if (input.nRefStates < 0)
        {
            throw std::invalid_argument("nRefStates >= 0");
        }
        if (input.N <= 0)
        {
            throw std::invalid_argument("N > 0");
        }
        if ((input.Lx*input.Ly*input.Lz) != input.N)
        {
            throw std::invalid_argument("Lx*Ly*Lz == N");
        }
        if (input.numOfNeighWS <= 0)
        {
            throw std::invalid_argument("numOfNeighWS > 0");
        }
        if (input.numOfEdgesBA <= 0)
        {
            throw std::invalid_argument("numOfEdgesBA > 0");
        }
        if (input.rewiringProb < 0.0)
        {
            throw std::invalid_argument("rewiringProb >= 0.0");
        }
        if (input.tTotal <= 0)
        {
            throw std::invalid_argument("tTotal > 0");
        }
        if (((input.seedIndex < 0) || (input.seedIndex >= input.N)) && (input.seedIndex != -1))
        {
            throw std::invalid_argument("0 <= seedIndex < N || seedIndex == -1 (random)");
        }
        if (input.ntheta <= 0)
        {
            throw std::invalid_argument("ntheta > 0");
        }
        if (input.thetaMin > input.thetaMax)
        {
            throw std::invalid_argument("thetaMin <= thetaMax");
        }
        return true;
    }

    InputParam_t GetInputArgs(int argCount, char* args[])
    {
        InputParam_t input = GetDefaultValues();
        std::string firstArg = std::string(args[1]);
        if ((firstArg == "h") || (firstArg == "-h") || (firstArg == "help") || (firstArg == "-help") || (firstArg == "--help")) // the user is asking for help
        {
            PrintHelp(args);
            exit(EXIT_FAILURE);
        }
        else
        {
            int i = 1;
            inputCmdLine = std::string(args[0]);
            while (i < argCount)
            {
                //std::string temp, tempNext;
                //std::stringstream ss1, ss2, ss3;
                //ss1 << args[i];
                //ss1 >> temp;
                std::string temp, tempNext;
                temp = std::string(args[i]);
                if ((firstArg == "h") || (firstArg == "-h") || (firstArg == "help") || (firstArg == "-help") || (firstArg == "--help"))
                {
                    PrintHelp(args);
                    exit(EXIT_FAILURE);
                }
                //if ((i + 1) >= argCount)
                if ((i + 1) < argCount)
                {
                    //PrintError("Missing value for " + temp);
                    //exit(EXIT_FAILURE);
                    tempNext = std::string(args[i + 1]);
                }
                else
                {
                    tempNext = "";
                }

                inputCmdLine = inputCmdLine + " " + temp;
                //ss3 << args[i + 1];
                //ss3 >> tempNext;
                //PrintMsg( ":: Ajustando -> " + temp + " = " + tempNext);
                if (!Misc::IsValidArgValue(tempNext))
                {
                    if (temp == par_simType)
                    {
                        input.simType = Misc::StrToEnum<SimulationType>(std::string(args[++i]), SimulationTypeChar);
                    }
                    else if (temp == par_seedType)
                    {
                        input.seedType = Misc::StrToEnum<SeedType>(std::string(args[++i]), SeedTypeChar);
                    }
                    else if (temp == par_nodeType)
                    {
                        input.nodeType = Misc::StrToEnum<Elements::ElementType>(std::string(args[++i]), Elements::ElementTypeChar);
                    }
                    else if (temp == par_coupType)
                    {
                        input.coupType = Misc::StrToEnum<Couplings::CouplingType>(std::string(args[++i]), Couplings::CouplingTypeChar);
                    }
                    else if (temp == par_netType)
                    {
                        input.netType = Misc::StrToEnum<Networks::AdjacencyMatrixType>(std::string(args[++i]), Networks::AdjacencyMatrixTypeChar);
                    }
                    else if (temp == par_weightType)
                    {
                        input.weightType = Misc::StrToEnum<Networks::AdjacencyMatrixWeightType>(std::string(args[++i]), Networks::AdjacencyMatrixWeightTypeChar);
                    }
                    else if (temp == par_outputType)
                    {
                        input.outputType = Misc::StrToEnum<Output::OutputType>(std::string(args[++i]), Output::OutputTypeChar);
                    }
                    else if (temp == par_netFile)
                    {
                        input.netFile = tempNext;
                        i++;
                    }
                    else if (temp == par_outFileSuf)
                    {
                        input.outFileSuf = tempNext;
                        i++;
                    }
                    else
                    {
                        if (!IsBooleanParameter(temp))
                        {
                            PrintError("Invalid value for argument " + temp + ". Impossible to set: " + temp + " = " + tempNext);
                            exit(EXIT_FAILURE);
                        }
                    }

                    if (!IsBooleanParameter(temp))
                    {
                        PrintMsg("# Set: " + temp + " = " + tempNext);
                        inputCmdLine = inputCmdLine + " " + tempNext;
                        i++;
                        continue;
                    }
                }
                else
                {
                    inputCmdLine = inputCmdLine + " " + tempNext;
                }
                if (temp == par_nSim)
                {
                    input.nSim = atoi(args[++i]);
                }
                else if (temp == par_theta)
                {
                    input.theta = atof(args[++i]);
                }
                else if (temp == par_nRefStates)
                {
                    input.nRefStates = atoi(args[++i]);
                }
                else if (temp == par_N)
                {
                    input.N = atoi(args[++i]);
                }
                else if (temp == par_Lx)
                {
                    input.Lx = atoi(args[++i]);
                }
                else if (temp == par_Ly)
                {
                    input.Ly = atoi(args[++i]);
                }
                else if (temp == par_Lz)
                {
                    input.Lz = atoi(args[++i]);
                }
                else if (temp == par_seedIndex)
                {
                    input.seedIndex = atoi(args[++i]);
                }
                else if (temp == par_numOfNeighWS)
                {
                    input.numOfNeighWS = atoi(args[++i]);
                }
                else if (temp == par_numOfEdgesBA)
                {
                    input.numOfEdgesBA = atoi(args[++i]);
                }
                else if (temp == par_rewiringProb)
                {
                    input.rewiringProb = atof(args[++i]);
                }
                else if (temp == par_weightAmp)
                {
                    input.weightAmp = atof(args[++i]);
                }
                else if (temp == par_tTotal)
                {
                    input.tTotal = atoi(args[++i]);
                }
                else if (temp == par_ntheta)
                {
                    input.ntheta = atoi(args[++i]);
                }
                else if (temp == par_thetaMin)
                {
                    input.thetaMin = atof(args[++i]);
                }
                else if (temp == par_thetaMax)
                {
                    input.thetaMax = atof(args[++i]);
                }
                else if (temp == par_isDirected)
                {
                    std::stringstream ss(tempNext);
                    int a;
                    if (!(ss >> a))
                    {
                        input.isDirected = true;
                        tempNext = "true";
                    }
                    else
                    {
                        ++i;
                        input.isDirected = a == 1;
                    }
                }
                else if (temp == par_isPeriodic)
                {
                    std::stringstream ss(tempNext);
                    int a;
                    if (!(ss >> a))
                    {
                        input.isPeriodic = true;
                        tempNext = "true";
                    }
                    else
                    {
                        ++i;
                        input.isPeriodic = a == 1;
                    }
                }
                else if (temp == par_writeOnRun)
                {
                    std::stringstream ss(tempNext);
                    int a;
                    if (!(ss >> a))
                    {
                        input.writeOnRun = true;
                        tempNext = "true";
                    }
                    else
                    {
                        ++i;
                        input.writeOnRun = a == 1;
                    }
                }
                else
                {
                    PrintError("Unrecognized parameter");
                    std::cout << "position  = " << i << "; parameter = " << temp << ";" << std::endl;
                    exit(EXIT_FAILURE);
                }
                PrintMsg("# Set: " + temp + " = " + tempNext);
                i++;
            }
        }
        if (!ValidateInputArgs(input))
        {
            PrintError("Error validating input parameters");
            exit(EXIT_FAILURE);
        }
        return input;//*/
    }
}

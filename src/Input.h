#ifndef _SIMINPUT_H
#define _SIMINPUT_H

#include <string>
#include <random>
#include "Elements.h"
#include "Couplings.h"
#include "Networks.h"
#include "Output.h"
//#include "Distribution_c.h"

namespace Output
{
    enum class OutputType;
}

namespace Elements
{
    enum class ElementType;
}

namespace Couplings
{
    enum class CouplingType;
}

namespace Networks
{
    enum class AdjacencyMatrixType;
}

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
    extern std::string inputCmdLine;

    enum class SimulationType
    {
        begin,
        Inactive,
        OneActive,
        TwoActive,
        nActive,
        FullyActive,
        end
    };
    static const char* SimulationTypeChar[] =
    {
        "begin",
        "Inactive",
        "OneActive",
        "TwoActive",
        "nActive",
        "FullyActive",
        "end"
    };

    enum class SeedType
    {
        begin,
        Cooperative,
        Competitive,
        end
    };
    static const char* SeedTypeChar[] =
    {
        "begin",
        "Cooperative",
        "Competitive",
        "end"
    };

    struct InputParam_t
    {
        SimulationType simType;
        SeedType seedType;
        int nSim;
        Elements::ElementType nodeType;
        double theta;
        int nRefStates;
        Couplings::CouplingType coupType;
        Networks::AdjacencyMatrixType netType;
        int N;
        int Lx;
        int Ly;
        int Lz;
        int numOfNeighWS;
        int numOfEdgesBA;
        double rewiringProb;
        bool isDirected;
        bool isPeriodic;
        std::string netFile;
        Networks::AdjacencyMatrixWeightType weightType;
        double weightAmp;
        int tTotal;
        bool writeOnRun;
        Output::OutputType outputType;
        int seedIndex;
        double thetaMin;
        double thetaMax;
        int ntheta;
        std::string outFileSuf;
    };

    static const std::string par_simType = "simType";
    static const std::string par_seedType = "seedType";
    static const std::string par_nSim = "nSim";
    static const std::string par_nodeType = "nodeType";
    static const std::string par_theta = "theta";
    static const std::string par_nRefStates = "nRefStates";
    static const std::string par_coupType = "coupType";
    static const std::string par_netType = "netType";
    static const std::string par_N = "N";
    static const std::string par_Lx = "Lx";
    static const std::string par_Ly = "Ly";
    static const std::string par_Lz = "Lz";
    static const std::string par_numOfNeighWS = "numOfNeigh";
    static const std::string par_numOfEdgesBA = "numOfEdgesBA";
    static const std::string par_rewiringProb = "rewiringProb";
    static const std::string par_isDirected = "directed";
    static const std::string par_isPeriodic = "periodic";
    static const std::string par_netFile = "netFile";
    static const std::string par_weightType = "weightType";
    static const std::string par_weightAmp = "weightAmp";
    static const std::string par_tTotal = "tTotal";
    static const std::string par_writeOnRun = "writeOnRun";
    static const std::string par_outputType = "outputType";
    static const std::string par_seedIndex = "seedIndex";
    static const std::string par_thetaMin = "thetaMin";
    static const std::string par_thetaMax = "thetaMax";
    static const std::string par_ntheta = "ntheta";
    static const std::string par_outFileSuf = "outFileSuf";

    // input parsing functions
    bool IsBooleanParameter(std::string&);
    void PrintHelp(char*[]);
    void PrintError(std::string);
    void PrintMsg(std::string);
    InputParam_t GetDefaultValues();
    InputParam_t GetInputArgs(int, char*[]);
    std::string GetInputParamString(InputParam_t);
    bool ValidateInputArgs(InputParam_t);
}

#endif
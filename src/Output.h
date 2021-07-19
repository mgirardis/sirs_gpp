#ifndef _SIMOUTPUT_H
#define _SIMOUTPUT_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include "Eigen/Core"
#include "Eigen/SparseCore" // Eigen; added by right click on project -> Properties -> C/C++ -> Additional Include Directories
#include "Misc.h"

namespace Output
{
    enum class OutputType
    {
        begin,
        NodesDynamics,
        NetworkDynamics,
        TotalNetworkSum,
        SpreadingMatrix,
        PhaseTransition,
        end
    };
    static const char* OutputTypeChar[] =
    {
        "begin",
        "NodesDynamics",
        "NetworkDynamics",
        "TotalNetworkSum",
        "SpreadingMatrix",
        "PhaseTransition",
        "end"
    };

    static inline std::fstream OpenFile(std::string& fileName, std::string headerTxt = "", bool verbose = true)
    {
        fileName = Misc::CheckAndGetFileName(fileName);
        if (verbose)
            std::cout << "Opening file to write: " << fileName << std::endl;
        std::fstream fs(fileName.c_str(), std::fstream::out | std::fstream::trunc);
        if (!headerTxt.empty())
            fs << headerTxt << std::endl;
        fs.precision(8);
        fs << std::scientific;
        return fs;
    }

    static inline void WriteFirstLine(std::string fileName, std::string line)
    {
        std::ifstream inputFile(fileName);
        std::string rfn = Misc::GetFileNameWithoutExtension(fileName) + ".tmp";
        std::fstream fs = OpenFile(rfn, "", false);
        fs << line << std::endl;
        fs << inputFile.rdbuf();
        fs.close();
        inputFile.close();
        std::remove(fileName.c_str());
        std::rename(rfn.c_str(), fileName.c_str());
    }

    static inline void WriteData(Eigen::SparseMatrix<double>& x, std::string& headerTxt, std::string& fileName)
    {
        fileName = Misc::CheckAndGetFileName(fileName);
        std::cout << "Writing " << fileName << std::endl;
        std::fstream fs(fileName.c_str(), std::fstream::out | std::fstream::trunc);
        if (!headerTxt.empty())
            fs << headerTxt << std::endl;
        fs << x.toDense();
        fs.close();
    }

    static inline void WriteData(Eigen::ArrayXd& x, std::string& headerTxt, std::string& fileName, bool firstColIsTime = false)
    {
        fileName = Misc::CheckAndGetFileName(fileName);
        std::cout << "Writing " << fileName << std::endl;
        std::fstream fs(fileName.c_str(), std::fstream::out | std::fstream::trunc);
        if (!headerTxt.empty())
            fs << headerTxt << std::endl;
        fs.precision(8);
        fs << std::scientific;
        if (firstColIsTime)
        {
            Eigen::ArrayXXd y(x.rows(), 2);
            y.col(0) = Eigen::ArrayXd::LinSpaced(x.rows(), 0.0, x.rows() - 1.0);
            y.col(1) = x;
            fs << y;
        }
        else
        {
            fs << x;
        }
        fs.close();
    }

    static inline void WriteData(Eigen::ArrayXXd& x, std::string& headerTxt, std::string& fileName, bool firstColIsTime = false)
    {
        fileName = Misc::CheckAndGetFileName(fileName);
        std::cout << "Writing " << fileName << std::endl;
        std::fstream fs(fileName.c_str(), std::fstream::out | std::fstream::trunc);
        if (!headerTxt.empty())
            fs << headerTxt << std::endl;
        fs.precision(8);
        fs << std::scientific;
        if (firstColIsTime)
        {
            int nCols = x.cols();
            x.conservativeResize(Eigen::NoChange, nCols + 1);
            Eigen::ArrayXXd y = x.leftCols(nCols);
            x.rightCols(nCols) = y;
            x.col(0) = Eigen::ArrayXd::LinSpaced(x.rows(), 0.0, x.rows() - 1.0);
        }
        fs << x;
        fs.close();
    }

    static inline void WriteData(std::vector<double>& x, std::string& headerTxt, std::string& fileName)
    {
        fileName = Misc::CheckAndGetFileName(fileName);
        std::cout << "Writing " << fileName << std::endl;
        std::fstream fs(fileName.c_str(), std::fstream::out | std::fstream::trunc);
        if (!headerTxt.empty())
            fs << headerTxt << std::endl;
        fs.precision(8);
        fs << std::scientific;
        int i = 0, n = (int)x.size();
        while (i < n)
        {
            fs << x[i] << "\t" << std::endl;
            i++;
        }
        fs.close();
    }

    static inline void WriteData(std::vector<double>& x, std::vector<double>& y, std::string& headerTxt, std::string& fileName)
    {
        fileName = Misc::CheckAndGetFileName(fileName);
        std::cout << "Writing " << fileName << std::endl;
        std::fstream fs(fileName.c_str(), std::fstream::out | std::fstream::trunc);
        if (!headerTxt.empty())
            fs << headerTxt << std::endl;
        fs.precision(8);
        fs << std::scientific;
        int i = 0, n = (int)x.size();
        while (i < n)
        {
            fs << x[i] << "\t" << y[i] << std::endl;
            i++;
        }
        fs.close();
    }

    static inline void WriteData(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::string& headerTxt, std::string& fileName)
    {
        fileName = Misc::CheckAndGetFileName(fileName);
        std::cout << "Writing " << fileName << std::endl;
        std::fstream fs(fileName.c_str(), std::fstream::out | std::fstream::trunc);
        if (!headerTxt.empty())
            fs << headerTxt << std::endl;
        fs.precision(8);
        fs << std::scientific;
        int i = 0, n = (int)x.size();
        while (i < n)
        {
            fs << x[i] << "\t" << y[i] << "\t" << z[i] << std::endl;
            i++;
        }
        fs.close();
    }
}

#endif
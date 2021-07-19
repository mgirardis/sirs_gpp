#include <memory>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <numeric>
#include <string>
#include <stdexcept>
#include <random>
#include <chrono>
#include "Eigen/Core"
#include "Eigen/SparseCore" // Eigen; added by right click on project -> Properties -> C/C++ -> Additional Include Directories
#include "Networks.h"
#include "Misc.h"
#include "Distribution_c.h"

namespace Networks
{
    //IAdjacencyMatrix* GetAdjacencyMatrix(AdjacencyMatrixParam_t par)
    std::unique_ptr<IAdjacencyMatrix> GetAdjacencyMatrix(AdjacencyMatrixParam_t par, std::mt19937_64& rand)
    {
        if (par.type == AdjacencyMatrixType::LinearLattice)
        {
            //return new LinearLatticeAdjacencyMatrix(par.nElems, par.numOfNeighbours, par.isDirected, par.isPeriodic);
            //return std::make_unique<LinearLatticeAdjacencyMatrix>(par.nElems, par.numOfNeighbours, par.isDirected, par.isPeriodic, par.weightType, par.weightAmp, rand);
            return std::unique_ptr<LinearLatticeAdjacencyMatrix>(new LinearLatticeAdjacencyMatrix(par.nElems, par.numOfNeighbours, par.isDirected, par.isPeriodic, par.weightType, par.weightAmp, rand));
        }
        else if (par.type == AdjacencyMatrixType::SquareLattice)
        {
            if (par.L[0] == par.L[1])
            {
                //return new SquareLatticeAdjacencyMatrix(par.nElems, par.numOfNeighbours, par.isDirected, par.isPeriodic);
                //return std::make_unique<SquareLatticeAdjacencyMatrix>(par.nElems, par.numOfNeighbours, par.isDirected, par.isPeriodic, par.weightType, par.weightAmp, rand);
                return std::unique_ptr<SquareLatticeAdjacencyMatrix>(new SquareLatticeAdjacencyMatrix(par.nElems, par.numOfNeighbours, par.isDirected, par.isPeriodic, par.weightType, par.weightAmp, rand));
            }
            else
            {
                std::vector<int> LL;
                LL.push_back(par.L[0]);
                LL.push_back(par.L[1]);
                //return new AdjacencyMatrixGridGraph(LL, par.isPeriodic, par.isDirected);
                //return std::make_unique<AdjacencyMatrixGridGraph>(LL, par.isPeriodic, par.isDirected, par.weightType, par.weightAmp, rand);
                return std::unique_ptr<AdjacencyMatrixGridGraph>(new AdjacencyMatrixGridGraph(LL, par.isPeriodic, par.isDirected, par.weightType, par.weightAmp, rand));
            }
        }
        else if (par.type == AdjacencyMatrixType::CompleteGraph)
        {
            //return new AdjacencyMatrixCompleteGraph(par.nElems, par.isDirected);
            //return std::make_unique<AdjacencyMatrixCompleteGraph>(par.nElems, par.isDirected, par.weightType, par.weightAmp, rand);
            return std::unique_ptr<AdjacencyMatrixCompleteGraph>(new AdjacencyMatrixCompleteGraph(par.nElems, par.isDirected, par.weightType, par.weightAmp, rand));
        }
        else if (par.type == AdjacencyMatrixType::RandomGraph)
        {
            //return new AdjacencyMatrixRandomGraph(par.nElems, par.rewiringProb, par.isDirected);
            //return std::make_unique<AdjacencyMatrixRandomGraph>(par.nElems, par.rewiringProb, par.isDirected, par.weightType, par.weightAmp, rand);
            return std::unique_ptr<AdjacencyMatrixRandomGraph>(new AdjacencyMatrixRandomGraph(par.nElems, par.rewiringProb, par.isDirected, par.weightType, par.weightAmp, rand));
        }
        else if (par.type == AdjacencyMatrixType::BarabasiAlbertGraph)
        {
            //return new AdjacencyMatrixBarabasiAlbert(par.nElems, par.numOfEdgesForNewElem, par.isDirected);
            //return std::make_unique<AdjacencyMatrixBarabasiAlbert>(par.nElems, par.numOfEdgesForNewElem, par.isDirected, par.weightType, par.weightAmp, rand);
            return std::unique_ptr<AdjacencyMatrixBarabasiAlbert>(new AdjacencyMatrixBarabasiAlbert(par.nElems, par.numOfEdgesForNewElem, par.isDirected, par.weightType, par.weightAmp, rand));
        }
        else if (par.type == AdjacencyMatrixType::WattsStrogatzGraph)
        {
            //return new AdjacencyMatrixWattsStrogatz(par.nElems, par.numOfNeighbours, par.rewiringProb, par.isDirected);
            //return std::make_unique<AdjacencyMatrixWattsStrogatz>(par.nElems, par.numOfNeighbours, par.rewiringProb, par.isDirected, par.weightType, par.weightAmp, rand);
            return std::unique_ptr<AdjacencyMatrixWattsStrogatz>(new AdjacencyMatrixWattsStrogatz(par.nElems, par.numOfNeighbours, par.rewiringProb, par.isDirected, par.weightType, par.weightAmp, rand));
        }
        else if (par.type == AdjacencyMatrixType::ConnectedWattsStrogatzGraph)
        {
            //return new AdjacencyMatrixWattsStrogatzConn(par.nElems, par.numOfNeighbours, par.rewiringProb, par.isDirected);
            //return std::make_unique<AdjacencyMatrixWattsStrogatzConn>(par.nElems, par.numOfNeighbours, par.rewiringProb, par.isDirected, par.weightType, par.weightAmp, rand);
            return std::unique_ptr<AdjacencyMatrixWattsStrogatzConn>(new AdjacencyMatrixWattsStrogatzConn(par.nElems, par.numOfNeighbours, par.rewiringProb, par.isDirected, par.weightType, par.weightAmp, rand));
        }
        else if (par.type == AdjacencyMatrixType::CubicLattice)
        {
            //return new AdjacencyMatrixGridGraph(par.L, par.isPeriodic, par.isDirected);
            //return std::make_unique<AdjacencyMatrixGridGraph>(par.L, par.isPeriodic, par.isDirected, par.weightType, par.weightAmp, rand);
            return std::unique_ptr<AdjacencyMatrixGridGraph>(new AdjacencyMatrixGridGraph(par.L, par.isPeriodic, par.isDirected, par.weightType, par.weightAmp, rand));
        }
        else if (par.type == AdjacencyMatrixType::FromFile)
        {
            if (par.netFileName == "") throw std::invalid_argument("No file specified for adjacency matrix");
            //return new AdjacencyMatrixFromFile(par.netFileName, par.isDirected);
            //return std::make_unique<AdjacencyMatrixFromFile>(par.netFileName, par.isDirected);
            return std::unique_ptr<AdjacencyMatrixFromFile>(new AdjacencyMatrixFromFile(par.netFileName, par.isDirected));
        }
        throw std::invalid_argument("Unrecognized AdjacencyMatrixType");//*/
    }

    IAdjacencyMatrix::IAdjacencyMatrix(AdjacencyMatrixWeightType wt, double wAmp, std::mt19937_64& rand)
        : rand(rand)
    {
        this->IsWeighted = wt != AdjacencyMatrixWeightType::Binary;
        this->weightType = wt;
        this->weightAmp = wAmp;
    }
    IAdjacencyMatrix::~IAdjacencyMatrix(){}
    Eigen::SparseMatrix<double> IAdjacencyMatrix::DeleteLowerTriangularElem(Eigen::SparseMatrix<double>& M)
    {
        return M.triangularView<Eigen::StrictlyUpper>();
    }
    Eigen::SparseMatrix<double> IAdjacencyMatrix::AssignWeights(Eigen::SparseMatrix<double>& M)
    {
        std::vector<double> s(M.nonZeros());
        if (this->weightType == AdjacencyMatrixWeightType::Homogeneous)
        {
            for (int j = 0; j < M.nonZeros(); j++) s[j] = this->weightAmp;
        }
        else if (this->weightType == AdjacencyMatrixWeightType::UniformRandom)
        {
            Distribution::DistributionParam_t dPar = Distribution::GetDefaultDistributionParam(Distribution::DistributionType::Uniform);
            dPar.x0 = 0.0;
            dPar.x1 = this->weightAmp;
            s = Distribution::GetDistribution(dPar, this->rand)->GetSample(M.nonZeros());
        }
        int i = 0;
        for (int k = 0; k < M.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
            {
                it.valueRef() = s[i];
                i++;
            }
        }
        return M;
    }

    AdjacencyMatrixGridGraph::AdjacencyMatrixGridGraph(std::vector<int> L, bool isPeriodic, bool isDirected, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : AdjacencyMatrixFromNetworkX("grid_graph", isDirected, wType, wAmp, rand)
    {
        this->NElems = std::accumulate(L.begin(), L.end(), 1, std::multiplies<int>());
        this->L = L;
        this->isPeriodic = isPeriodic;
        this->CallNetworkX();
        AdjacencyMatrixFromFile::Initialize(this->FileName);
    }
    AdjacencyMatrixGridGraph::~AdjacencyMatrixGridGraph(){}
    std::string AdjacencyMatrixGridGraph::GetNetworkXGeneratorMethodCall()
    {
        std::stringstream ss;
        ss.precision(16);
        int i = 0, nn = this->L.size();
        while (i < nn)
        {
            ss << this->L[i] << ",";
            i++;
        }
        std::string Lstr = ss.str();
        Lstr = Lstr.substr(0, Lstr.size() - 1);
        ss.str("");
        ss.clear();
        ss << "G=nx." << this->GraphMethod << "(dim=[" << Lstr << "],periodic=" << (this->isPeriodic ? "True" : "False") << ")";
        return ss.str();
    }
    std::string AdjacencyMatrixGridGraph::GetNetworkXAdjMatrix()
    {
        return "nx.adj_matrix(G,nodelist=sorted(G.nodes()))";
    }

    AdjacencyMatrixCompleteGraph::AdjacencyMatrixCompleteGraph(int numOfElems, bool isDirected, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : AdjacencyMatrixFromNetworkX("complete_graph", isDirected, wType, wAmp, rand)
    {
        this->NElems = numOfElems;
        this->CallNetworkX();
        AdjacencyMatrixFromFile::Initialize(this->FileName);
    }
    AdjacencyMatrixCompleteGraph::~AdjacencyMatrixCompleteGraph(){}
    std::string AdjacencyMatrixCompleteGraph::GetNetworkXGeneratorMethodCall()
    {
        return AdjacencyMatrixFromNetworkX::GetNetworkXGeneratorMethodCall();
    }
    std::string AdjacencyMatrixCompleteGraph::GetNetworkXAdjMatrix()
    {
        return AdjacencyMatrixFromNetworkX::GetNetworkXAdjMatrix();
    }

    AdjacencyMatrixWattsStrogatzConn::AdjacencyMatrixWattsStrogatzConn(int numOfElems, int numOfNeighbours, double rewiringProb, bool isDirected, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : AdjacencyMatrixFromNetworkX("connected_watts_strogatz_graph", isDirected, wType, wAmp, rand)
    {
        this->NElems = numOfElems;
        this->k = numOfNeighbours;
        this->p = rewiringProb;
        this->CallNetworkX();
        AdjacencyMatrixFromFile::Initialize(this->FileName);
    }
    AdjacencyMatrixWattsStrogatzConn::~AdjacencyMatrixWattsStrogatzConn(){}
    std::string AdjacencyMatrixWattsStrogatzConn::GetNetworkXGeneratorMethodCall()
    {
        std::stringstream ss;
        ss.precision(16);
        ss << "G=nx." << this->GraphMethod << "(" << this->NElems << "," << this->k << "," << this->p << ")";
        return ss.str();
    }
    std::string AdjacencyMatrixWattsStrogatzConn::GetNetworkXAdjMatrix()
    {
        return AdjacencyMatrixFromNetworkX::GetNetworkXAdjMatrix();
    }

    AdjacencyMatrixWattsStrogatz::AdjacencyMatrixWattsStrogatz(int numOfElems, int numOfNeighbours, double rewiringProb, bool isDirected, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : AdjacencyMatrixFromNetworkX("watts_strogatz_graph", isDirected, wType, wAmp, rand)
    {
        this->NElems = numOfElems;
        this->k = numOfNeighbours;
        this->p = rewiringProb;
        this->CallNetworkX();
        AdjacencyMatrixFromFile::Initialize(this->FileName);
    }
    AdjacencyMatrixWattsStrogatz::~AdjacencyMatrixWattsStrogatz(){}
    std::string AdjacencyMatrixWattsStrogatz::GetNetworkXGeneratorMethodCall()
    {
        std::stringstream ss;
        ss.precision(16);
        ss << "G=nx." << this->GraphMethod << "(" << this->NElems << "," << this->k << "," << this->p << ")";
        return ss.str();
    }
    std::string AdjacencyMatrixWattsStrogatz::GetNetworkXAdjMatrix()
    {
        return AdjacencyMatrixFromNetworkX::GetNetworkXAdjMatrix();
    }

    AdjacencyMatrixRandomGraph::AdjacencyMatrixRandomGraph(int numOfElems, double pOfEdgeCreation, bool isDirected, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : AdjacencyMatrixFromNetworkX("fast_gnp_random_graph", isDirected, wType, wAmp, rand)
    {
        this->NElems = numOfElems;
        this->p = pOfEdgeCreation;
        this->CallNetworkX();
        AdjacencyMatrixFromFile::Initialize(this->FileName);
    }
    AdjacencyMatrixRandomGraph::~AdjacencyMatrixRandomGraph(){}
    std::string AdjacencyMatrixRandomGraph::GetNetworkXGeneratorMethodCall()
    {
        std::stringstream ss;
        ss.precision(16);
        ss << "G=nx." << this->GraphMethod << "(" << this->NElems << "," << this->p << ", directed=" << (this->IsDirected ? "True" : "False") << ")";
        return ss.str();
    }
    std::string AdjacencyMatrixRandomGraph::GetNetworkXAdjMatrix()
    {
        return AdjacencyMatrixFromNetworkX::GetNetworkXAdjMatrix();
    }

    AdjacencyMatrixBarabasiAlbert::AdjacencyMatrixBarabasiAlbert(int numOfElems, int numOfEdgesForNewElem, bool isDirected, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : AdjacencyMatrixFromNetworkX("barabasi_albert_graph", isDirected, wType, wAmp, rand)
    {
        this->NElems = numOfElems;
        this->m = numOfEdgesForNewElem;
        this->CallNetworkX();
        AdjacencyMatrixFromFile::Initialize(this->FileName);
    }
    AdjacencyMatrixBarabasiAlbert::~AdjacencyMatrixBarabasiAlbert(){}
    std::string AdjacencyMatrixBarabasiAlbert::GetNetworkXGeneratorMethodCall()
    {
        std::stringstream ss;
        ss.precision(16);
        ss << "G=nx." << this->GraphMethod << "(" << this->NElems << "," << this->m << ")";
        return ss.str();
    }
    std::string AdjacencyMatrixBarabasiAlbert::GetNetworkXAdjMatrix()
    {
        return AdjacencyMatrixFromNetworkX::GetNetworkXAdjMatrix();
    }

    AdjacencyMatrixFromNetworkX::AdjacencyMatrixFromNetworkX(std::string graphMethod, bool isDirected, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : AdjacencyMatrixFromFile(isDirected, wType, wAmp, rand)
    {
        this->GraphMethod = graphMethod;
        this->FileName = this->GraphMethod + ".tmp";
        this->FileName = Misc::CheckAndGetFileName(this->FileName);
    }
    AdjacencyMatrixFromNetworkX::~AdjacencyMatrixFromNetworkX(){}
    std::string AdjacencyMatrixFromNetworkX::GetNetworkXGeneratorMethodCall()
    {
        std::stringstream ss;
        ss.precision(16);
        ss << "G=nx." << this->GraphMethod << "(" << this->NElems << ")";
        return ss.str();
    }
    std::string AdjacencyMatrixFromNetworkX::GetNetworkXAdjMatrix()
    {
        return "nx.adj_matrix(G)";
    }
    void AdjacencyMatrixFromNetworkX::CallNetworkX()
    {
        std::string prefix = ""; // prefix used for linux or Mac OS X
#ifdef WIN32
        prefix = "cmd /c "; // in windows, one must first call the cmd terminal to execute python commands
#endif
        std::stringstream ss;
        ss << prefix << "python -c \"import networkx as nx; import numpy as np; " << this->GetNetworkXGeneratorMethodCall() << "; np.savetxt('" << this->FileName << "', " << this->GetNetworkXAdjMatrix() <<",fmt='%1d')\"";
        std::string op;
        // creating process
#ifndef _DEBUG
        try
        {
#endif
            op = Misc::RunProcess(ss.str());
#ifndef _DEBUG
        }
        catch (const std::exception& e)
        {
            std::cout << "EXCEPTION in AdjacencyMatrixFromNetworkX::CallNetworkX" << std::endl << e.what();
        }
#endif
        if (Misc::StringContains(op, static_cast<std::string>("Traceback")))
        {
            std::cout << "ERROR Runing networkX..." << std::endl << op;
            throw std::runtime_error("An error occurred during the graph generation process...");
        }//*/
    }
    Eigen::SparseMatrix<double> AdjacencyMatrixFromNetworkX::BuildAndGetMatrix()
    {
        if (!Misc::FileExists(this->FileName))
            this->CallNetworkX();
        Eigen::SparseMatrix<double> adj = AdjacencyMatrixFromFile::BuildAndGetMatrix();
        if (Misc::FileExists(this->FileName))
        {
            if (remove(this->FileName.c_str()) != 0)
            {
                std::cout << "WARNING: file " << this->FileName << " could not be removed..." << std::endl;
            }
        }
        return adj;
    }

    AdjacencyMatrixFromFile::AdjacencyMatrixFromFile(void)
        : IAdjacencyMatrix(AdjacencyMatrixWeightType::Homogeneous, 1.0, std::mt19937_64())
    {
        this->IsFromFile = true;
        this->NElems = -1;
    }
    AdjacencyMatrixFromFile::AdjacencyMatrixFromFile(std::string fileName, bool isDirected)
        : IAdjacencyMatrix(AdjacencyMatrixWeightType::Homogeneous, 1.0, std::mt19937_64())
    {
        this->IsDirected = isDirected;
        this->IsFromFile = true;
        this->Initialize(fileName);
    }
    AdjacencyMatrixFromFile::AdjacencyMatrixFromFile(bool isDirected, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : IAdjacencyMatrix(wType, wAmp, rand)
    {
        this->IsFromFile = false;
        this->IsDirected = isDirected;
    }
    AdjacencyMatrixFromFile::~AdjacencyMatrixFromFile(){}
    void AdjacencyMatrixFromFile::SetNElems(void)
    {
        std::ifstream fs;
        std::string line;
        fs.open(this->FileName);
        while (!fs.eof())
        {
            std::getline(fs, line);
            line = Misc::StringTrimAndReduce(line);
            if ((line[0] == '#') || (line.empty()))
                continue;
            else
                break;
        }
        this->NElems = Misc::StringSplit(line, static_cast<std::string>(" ")).size();
        fs.close();
    }
    void AdjacencyMatrixFromFile::Initialize(std::string fileName)
    {
        this->FileName = fileName;
#ifndef _DEBUG
        try
        {
#endif
            this->SetNElems();
#ifndef _DEBUG
        }
        catch (const std::exception& e)
        {
            std::cout << "EXCEPTION in AdjacencyMatrixFromFile::Initialize" << std::endl << e.what();
        }
#endif
    }
    Eigen::SparseMatrix<double> AdjacencyMatrixFromFile::BuildAndGetMatrix()
    {
        if (this->NElems == -1)
            throw std::invalid_argument("AdjacencyMatrixFromFile has not been initialized properly!");

        Eigen::SparseMatrix<double> adj(this->NElems, this->NElems);
        std::vector<Eigen::Triplet<double>> el;
#ifndef _DEBUG
        try
        {
#endif
            // openning temporary file for reading
            std::ifstream fs;
            fs.open(this->FileName);

            std::vector<std::string> col;
            std::string line;
            double val;
            int j, n, i = 0;
            while (!fs.eof())
            {
                std::getline(fs, line);
                line = Misc::StringTrimAndReduce(line);
                if ((line[0] == '#') || (line.empty()))
                    continue;
                col = Misc::StringSplit(line, " ");
                n = col.size();
                j = 0;
                while (j < n)
                {
                    val = atof(col[j].c_str());
                    if (val != 0.0)
                    {
                        el.push_back(Eigen::Triplet<double>(i, j, val));
                    }
                    j++;
                }
                i++;
            }
            fs.close();
            adj.setFromTriplets(el.begin(), el.end());
#ifndef _DEBUG
        }
        catch (const std::exception& e)
        {
            std::cout << "EXCEPTION in AdjacencyMatrixFromFile::BuildAndGetMatrix" << std::endl << e.what();
        }
#endif
        if (this->IsDirected) // transforms the matrix into a triangular matrix by simply deleting the lower diagonal
        {
            adj = this->DeleteLowerTriangularElem(adj);
        }
        if ((this->IsWeighted) && (!this->IsFromFile))
        {
            adj = this->AssignWeights(adj);
        }
        return adj;
    }

    LatticeAdjacencyMatrix::LatticeAdjacencyMatrix(int nElems, int nElemsOnARow, int nNeighbours, bool isDirected, bool isPeriodic, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : IAdjacencyMatrix(wType, wAmp, rand)
    {
        this->NNeighbours = nNeighbours;
        this->NElemsOnARow = nElemsOnARow;
        this->NElems = nElems;
        this->IsDirected = isDirected;
        this->IsPeriodic = isPeriodic;
    }
    LatticeAdjacencyMatrix::~LatticeAdjacencyMatrix(){}

    SquareLatticeAdjacencyMatrix::SquareLatticeAdjacencyMatrix(int nElems, int nNeighbours, bool isDirected, bool isPeriodic, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : LatticeAdjacencyMatrix(nElems, nElems, 4, isDirected, isPeriodic, wType, wAmp, rand)
    {
        double s = std::sqrt(nElems);
        if (std::floor(s) == s)
        {
            this->NElemsOnARow = static_cast<int>(s);
        }
        else
        {
            throw std::invalid_argument("The specified number of elements does not have an exact sqrt, so it cannot form a sqr net");
        }
    }
    SquareLatticeAdjacencyMatrix::~SquareLatticeAdjacencyMatrix(){}
    Eigen::SparseMatrix<double> SquareLatticeAdjacencyMatrix::BuildAndGetMatrix()
    {
        if (this->IsPeriodic)
            return this->BuildAndGetPeriodicMatrix();
        else
            return this->BuildAndGetFreeMatrix();
    }
    Eigen::SparseMatrix<double> SquareLatticeAdjacencyMatrix::BuildAndGetPeriodicMatrix()
    {
        Eigen::SparseMatrix<double> matrix(this->NElems, this->NElems);
        matrix.reserve(Eigen::VectorXi::Constant(this->NElems, 2 * this->NNeighbours));
        int i, j, k, m;
        i = 0;
        int* n = new int[this->NNeighbours]; // the index of the neurons connected to each k neuron (i is the row # and j is the column # of the network matrix, which is LxL, here L = nNeuronsOnARow)
        // the first row doesn't connect to n0, the last row doesn't connect to n1, the first column doesn't connect to n2, the last column doesn't connect to n3
        // n0 = (i-1) * L + j
        // n1 = (i+1) * L + j
        // n2 = i * L + (j-1)
        // n3 = i * L + (j+1)
        while (i < this->NElemsOnARow)
        {
            j = 0;
            while (j < this->NElemsOnARow)
            {
                // the index of the neuron at site i,j on the network
                k = i * this->NElemsOnARow + j;

                // the indeces of the neurons around it
                n[0] = (i - 1) * this->NElemsOnARow + j;
                n[1] = (i + 1) * this->NElemsOnARow + j;
                n[2] = i * this->NElemsOnARow + (j - 1);
                n[3] = i * this->NElemsOnARow + (j + 1);

                // the matrix should be symmetric (Aij = Aji, because every connection is a two-way connection)
                m = 0;
                while (m < this->NNeighbours)
                {
                    if ((i == 0) && (m == 0)) // first row connects to last row
                    {
                        n[0] = (this->NElemsOnARow - 1) * this->NElemsOnARow + j; // i = nNeuronsOnARow - 1
                    }
                    if ((i == (this->NElemsOnARow - 1)) && (m == 1)) // last row connects to the first row
                    {
                        n[1] = j; // i = 0
                    }
                    if ((j == 0) && (m == 2)) // first column connects to the last column
                    {
                        n[2] = i * this->NElemsOnARow + (this->NElemsOnARow - 1); // j = nNeuronsOnARow - 1
                    }
                    if ((j == (this->NElemsOnARow - 1)) && (m == 3)) // last column connects to the first column
                    {
                        n[3] = i * this->NElemsOnARow; // j = 0
                    }
                    matrix.insert(k, n[m]) = 1.0;
                    //AMatrix.addElement(n[m], k, true);
                    m++;
                }
                j++;
            }
            i++;
        }
        if (this->IsDirected) // transforms the matrix into a triangular matrix by simply deleting the lower diagonal
        {
            matrix = this->DeleteLowerTriangularElem(matrix);
        }
        if (this->IsWeighted)
        {
            matrix = this->AssignWeights(matrix);
        }
        delete n;
        return matrix;
    }
    Eigen::SparseMatrix<double> SquareLatticeAdjacencyMatrix::BuildAndGetFreeMatrix()
    {
        Eigen::SparseMatrix<double> matrix(this->NElems, this->NElems);
        matrix.reserve(Eigen::VectorXi::Constant(this->NElems, 2 * this->NNeighbours));
        int i, j, k, m;
        i = 0;
        int* n = new int[this->NNeighbours]; // the index of the neurons connected to each k neuron (i is the row # and j is the column # of the network matrix, which is LxL, here L = nNeuronsOnARow)
        // the first row doesn't connect to n0, the last row doesn't connect to n1, the first column doesn't connect to n2, the last column doesn't connect to n3
        // n0 = (i-1) * L + j
        // n1 = (i+1) * L + j
        // n2 = i * L + (j-1)
        // n3 = i * L + (j+1)
        while (i < this->NElemsOnARow)
        {
            j = 0;
            while (j < this->NElemsOnARow)
            {
                // the index of the neuron at site i,j on the network
                k = i * this->NElemsOnARow + j;

                // the indeces of the neurons around it
                n[0] = (i - 1) * this->NElemsOnARow + j;
                n[1] = (i + 1) * this->NElemsOnARow + j;
                n[2] = i * this->NElemsOnARow + (j - 1);
                n[3] = i * this->NElemsOnARow + (j + 1);

                // the matrix should be symmetric (Aij = Aji, because every connection is a two-way connection)
                m = 0;
                while (m < this->NNeighbours)
                {
                    if ((i == 0) && (m == 0)) // first row doesn't connect to n0
                    {
                        m++;
                        continue;
                    }
                    if ((i == (this->NElemsOnARow - 1)) && (m == 1)) // last row doesn't connect to n1
                    {
                        m++;
                        continue;
                    }
                    if ((j == 0) && (m == 2)) // first column doesn't connect to n2
                    {
                        m++;
                        continue;
                    }
                    if ((j == (this->NElemsOnARow - 1)) && (m == 3)) // last column doesn't connect to n3
                    {
                        m++;
                        continue;
                    }
                    matrix.insert(k, n[m]) = 1.0;
                    //AMatrix.addElement(n[m], k, true);
                    m++;
                }
                j++;
            }
            i++;
        }
        if (this->IsDirected) // transforms the matrix into a diagonal matrix by simply deleting the lower diagonal
        {
            matrix = this->DeleteLowerTriangularElem(matrix);
        }
        if (this->IsWeighted)
        {
            matrix = this->AssignWeights(matrix);
        }
        delete n;
        return matrix;
    }

    LinearLatticeAdjacencyMatrix::LinearLatticeAdjacencyMatrix(int nElems, int nNeighbours, bool isDirected, bool isPeriodic, AdjacencyMatrixWeightType wType, double wAmp, std::mt19937_64& rand)
        : LatticeAdjacencyMatrix(nElems, nElems, 2, isDirected, isPeriodic, wType, wAmp, rand)
    {}
    LinearLatticeAdjacencyMatrix::~LinearLatticeAdjacencyMatrix(){}
    Eigen::SparseMatrix<double> LinearLatticeAdjacencyMatrix::BuildAndGetMatrix()
    {
        if (this->IsPeriodic)
            return this->BuildAndGetPeriodicMatrix();
        else
            return this->BuildAndGetFreeMatrix();
    }
    Eigen::SparseMatrix<double> LinearLatticeAdjacencyMatrix::BuildAndGetPeriodicMatrix()
    {
        Eigen::SparseMatrix<double> matrix(this->NElems, this->NElems);
        matrix.reserve(Eigen::VectorXi::Constant(this->NElems, 2 * this->NNeighbours));
        int i = 0;
        int* n = new int[this->NNeighbours]; // the index of the neurons connected to each k neuron (i is the row # and j is the column # of the network matrix, which is LxL, here L = nNeuronsOnARow)
        // the first row doesn't connect to n0, the last row doesn't connect to n1, the first column doesn't connect to n2, the last column doesn't connect to n3
        // n0 = (i-1) * L + j
        // n1 = (i+1) * L + j
        // n2 = i * L + (j-1)
        // n3 = i * L + (j+1)
        while (i < this->NElems)
        {
            if (i != 0)
            {
                n[0] = i - 1;
            }
            else
            {
                n[0] = NElems - 1;
            }
            if (i != NElems - 1)
            {
                n[1] = i + 1;
            }
            else
            {
                n[1] = 0;
            }
            matrix.insert(i, n[0]) = 1.0;
            matrix.insert(i, n[1]) = 1.0;
            i++;
        }
        if (this->IsDirected) // transforms the matrix into a diagonal matrix by simply deleting the lower diagonal
        {
            matrix = this->DeleteLowerTriangularElem(matrix);
        }
        if (this->IsWeighted)
        {
            matrix = this->AssignWeights(matrix);
        }
        delete n;
        return matrix;
    }
    Eigen::SparseMatrix<double> LinearLatticeAdjacencyMatrix::BuildAndGetFreeMatrix()
    {
        Eigen::SparseMatrix<double> matrix(this->NElems, this->NElems);
        matrix.reserve(Eigen::VectorXi::Constant(this->NElems, 2 * this->NNeighbours));
        int i = 0;
        int* n = new int[this->NNeighbours]; // the index of the neurons connected to each k neuron (i is the row # and j is the column # of the network matrix, which is LxL, here L = nNeuronsOnARow)
        // the first row doesn't connect to n0, the last row doesn't connect to n1, the first column doesn't connect to n2, the last column doesn't connect to n3
        // n0 = (i-1) * L + j
        // n1 = (i+1) * L + j
        // n2 = i * L + (j-1)
        // n3 = i * L + (j+1)
        while (i < this->NElems)
        {
            if (i != 0)
            {
                n[0] = i - 1; // (i - 1) * nNeurons + j;
                matrix.insert(i, n[0]) = 1.0;
            }
            if (i != NElems - 1)
            {
                n[1] = i + 1; // (i + 1) * nNeurons + j;
                matrix.insert(i, n[1]) = 1.0;
            }
            i++;
        }
        if (this->IsDirected) // transforms the matrix into a diagonal matrix by simply deleting the lower diagonal
        {
            matrix = this->DeleteLowerTriangularElem(matrix);
        }
        if (this->IsWeighted)
        {
            matrix = this->AssignWeights(matrix);
        }
        delete n;
        return matrix;
    }
}
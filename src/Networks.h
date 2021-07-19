#ifndef _NETWORKS_H
#define _NETWORKS_H

#include <memory>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include "Eigen/Core"
#include "Eigen/SparseCore" // Eigen; added by right click on project -> Properties -> C/C++ -> Additional Include Directories

namespace Networks
{
    enum class AdjacencyMatrixType
    {
        begin,
        /// <summary>
        /// generates a regular lattice with linear geometry, 2 neighbours and periodic boundary condition
        /// </summary>
        LinearLattice,

        /// <summary>
        /// generates a regular lattice with square geometry, 4 neighbours and free boundary condition
        /// </summary>
        SquareLattice,

        /// <summary>
        /// generates a cubic lattic with periodic boundary condition - networkX
        /// </summary>
        CubicLattice,

        /// <summary>
        /// generates a mean-field graph (every element is connected to all the others -- networkX)
        /// </summary>
        CompleteGraph,

        /// <summary>
        /// generates a random graph (networkX)
        /// </summary>
        RandomGraph,

        /// <summary>
        /// generates a Barabasi-Albert Scale-free graph (networkX)
        /// </summary>
        BarabasiAlbertGraph,

        /// <summary>
        /// generates a Watts-Strogatz Small-World graph (it's not guaranteed that all elements are connected in one tree... there may be lonely clusters - poor of them -- networkX)
        /// </summary>
        WattsStrogatzGraph,

        /// <summary>
        /// generates a Watts-Strogatz Small-World graph with every element connected in a single tree... The process can throw an exception if the amount of tries is exceeded (which is 100 -- networkX)
        /// </summary>
        ConnectedWattsStrogatzGraph,

        /// <summary>
        /// specifies that the adjacency matrix should be read from file
        /// </summary>
        FromFile,
        end
    };
    static const char* AdjacencyMatrixTypeChar[] =
    {
        "begin",
        "LinearLattice",
        "SquareLattice",
        "CubicLattice",
        "CompleteGraph",
        "RandomGraph",
        "BarabasiAlbertGraph",
        "WattsStrogatzGraph",
        "ConnectedWattsStrogatzGraph",
        "FromFile",
        "end"
    };

    enum class AdjacencyMatrixWeightType
    {
        begin,
        Binary,
        Homogeneous,
        UniformRandom,
        end
    };
    static const char* AdjacencyMatrixWeightTypeChar[] =
    {
        "begin",
        "Binary",
        "Homogeneous",
        "UniformRandom",
        "end"
    };

    struct AdjacencyMatrixParam_t
    {
        AdjacencyMatrixType type;
        int nElems;
        std::vector<int> L;
        int numOfNeighbours;
        int numOfEdgesForNewElem;
        double rewiringProb;
        bool isDirected;
        bool isPeriodic;
        AdjacencyMatrixWeightType weightType;
        double weightAmp;
        std::string netFileName;
    };

    class IAdjacencyMatrix
    {
    protected:
        /// <summary>
        /// number of elements
        /// </summary>
        int NElems;
        bool IsWeighted; // determines if is weighted adjacency matrix
        double weightAmp;
        AdjacencyMatrixWeightType weightType;
        //std::function<Eigen::SparseMatrix<double>(void)> BuildAndGetMatrix_m; // ponteiro para uma funcao
        std::mt19937_64& rand; //random number generator to create randomly distributed weights

    public:
        IAdjacencyMatrix(AdjacencyMatrixWeightType, double, std::mt19937_64&);
        ~IAdjacencyMatrix(void);
        /// <summary>
        /// builds and gets the matrix (sparse matrix)
        /// </summary>
        /// <returns>a sparse matrix of type ElemType containing the adjacency matrix</returns>
        virtual Eigen::SparseMatrix<double> BuildAndGetMatrix() = 0;
        
        
        /// <summary>
        /// utilized to create directed graphs
        /// </summary>
        Eigen::SparseMatrix<double> DeleteLowerTriangularElem(Eigen::SparseMatrix<double>&);

        /// <summary>
        /// assign weigths to the matrix entries according to the selected distribution type
        /// </summary>
        Eigen::SparseMatrix<double> AssignWeights(Eigen::SparseMatrix<double>&);
    };

    /// <summary>
    /// Get the specified adjacency matrix...
    /// </summary>
    /// <param name="type">the type of the desired matrix</param>
    /// <param name="nElems">total number of elements in the network</param>
    /// <param name="L">number of elements in each direction of the network (useful for lattice networks)</param>
    /// <param name="numOfNeighbours">number of neighbours of each element (for a Watts-Strogatz, it's the initial number of neighbours of each element)</param>
    /// <param name="numOfEdgesForNewElem">the "m" parameter for a Barabasi-Albert Graph (num of edges for each new attached node) - only important for Barabasi-Albert Graph</param>
    /// <param name="rewiringProb">the probability to replace a connection with another one (only important to create Watts-Strogatz graph)</param>
    /// <param name="netFileName">name of the file with the adjacency matrix, in case of FromFile selected</param>
    /// <returns>the desired adjacency matrix constructor</returns>
    ///                                    type, nElems, L, numOfNeighbours, numOfEdgesForNewElem, rewiringProb, isDirected, netFileName
    //IAdjacencyMatrix* GetAdjacencyMatrix(AdjacencyMatrixParam_t);
    std::unique_ptr<IAdjacencyMatrix> GetAdjacencyMatrix(AdjacencyMatrixParam_t, std::mt19937_64&);

    class LatticeAdjacencyMatrix : public IAdjacencyMatrix
    {
    protected:
        int NElemsOnARow;
        int NNeighbours;
        bool IsDirected;
        bool IsPeriodic;
    public:
        LatticeAdjacencyMatrix(int, int, int, bool, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&); //(int nElems, int nElemsOnARow, int nNeighbours, bool isDirected, bool isPeriodic)
        ~LatticeAdjacencyMatrix();
        virtual Eigen::SparseMatrix<double> BuildAndGetMatrix() = 0;
    };

    class LinearLatticeAdjacencyMatrix : public LatticeAdjacencyMatrix
    {
    private:
        Eigen::SparseMatrix<double> BuildAndGetPeriodicMatrix();
        Eigen::SparseMatrix<double> BuildAndGetFreeMatrix();
    public:
        LinearLatticeAdjacencyMatrix(int, int, bool, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&);// nElems, nNeighbours, isDirected, isPeriodic
        ~LinearLatticeAdjacencyMatrix();
        virtual Eigen::SparseMatrix<double> BuildAndGetMatrix();
    };

    class SquareLatticeAdjacencyMatrix : public LatticeAdjacencyMatrix
    {
    private:
        Eigen::SparseMatrix<double> BuildAndGetPeriodicMatrix();
        Eigen::SparseMatrix<double> BuildAndGetFreeMatrix();
    public:
        SquareLatticeAdjacencyMatrix(int, int, bool, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&); // (int nElems, int nNeighbours, bool isDirected, bool isPeriodic);
        ~SquareLatticeAdjacencyMatrix();
        virtual Eigen::SparseMatrix<double> BuildAndGetMatrix();
    };

    class AdjacencyMatrixFromFile : public IAdjacencyMatrix
    {
    private:
        /// <summary>
        /// adjusts the amount of elements within this adjacency matrix by counting columns of the first line of the specified file
        /// </summary>
        void SetNElems(void);
        bool IsFromFile; // determines whether this adjacency matrix is derectly from an input file

    protected:
        bool IsDirected;
        /// <summary>
        /// the temporary filename where will be stored the adjacency matrix
        /// </summary>
        std::string FileName;
        /// <summary>
        /// empty constructor
        /// </summary>
        AdjacencyMatrixFromFile(void);
        /// <summary>
        /// initializes this adjacency matrix if it has been inherited
        /// </summary>
        /// <param name="fileName">the name of the file with the adjacency matrix</param>
        void Initialize(std::string); // fileName
    public:
        /// <summary>
        /// constructor of the base class
        /// </summary>
        /// <param name="fileName">the name of the file with the adjacency matrix</param>
        AdjacencyMatrixFromFile(std::string, bool); //(String fileName, Boolean isDirected);
        AdjacencyMatrixFromFile(bool, AdjacencyMatrixWeightType, double, std::mt19937_64&); // isDirected
        ~AdjacencyMatrixFromFile();
        /// <summary>
        /// builds and gets the matrix (sparse matrix)
        /// </summary>
        /// <returns>a sparse matrix of type T containing the adjacency matrix</returns>
        virtual Eigen::SparseMatrix<double> BuildAndGetMatrix();
    };

    class AdjacencyMatrixFromNetworkX : public AdjacencyMatrixFromFile
    {
    protected:
        /// <summary>
        /// networkX method used to build the graph
        /// </summary>
        std::string GraphMethod;
        /// <summary>
        /// generates the method call according to the chosen GraphMethod and to the NetworkX Docs
        /// </summary>
        /// <returns>a string containing the NetworkX method call</returns>
        virtual std::string GetNetworkXGeneratorMethodCall();
        /// <summary>
        /// the python command line to get the adjacency matrix of the selected method from networkx
        /// </summary>
        /// <returns>the python command line to get the adjacency matrix of the selected method from networkx</returns>
        virtual std::string GetNetworkXAdjMatrix();
        /// <summary>
        /// runs the process with the networkx
        /// </summary>
        void CallNetworkX();

    public:
        /// <summary>
        /// constructor of the base class
        /// </summary>
        /// <param name="graphMethod">the graph generator method as specified by the NetworkX Documentation</param>
        AdjacencyMatrixFromNetworkX(std::string, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&); //(std::string graphMethod, bool isDirected);
        ~AdjacencyMatrixFromNetworkX();
        virtual Eigen::SparseMatrix<double> BuildAndGetMatrix();
    };

    class AdjacencyMatrixBarabasiAlbert : public AdjacencyMatrixFromNetworkX
    {
    private:
        /// <summary>
        /// number of edges that a new attached element will have initially
        /// </summary>
        int m;
    protected:
        virtual std::string GetNetworkXGeneratorMethodCall();
        virtual std::string GetNetworkXAdjMatrix();
    public:
        /// <summary>
        /// constructor - creates the temporary file with the graph's adj matrix
        /// </summary>
        /// <param name="numOfElems">number of elements in the graph</param>
        /// <param name="numOfEdgesForNewElem">initial number of connections that a newly attached element will have</param>
        AdjacencyMatrixBarabasiAlbert(int, int, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&); //(int numOfElems, int numOfEdgesForNewElem, bool isDirected)
        ~AdjacencyMatrixBarabasiAlbert();
    };

    class AdjacencyMatrixRandomGraph : public AdjacencyMatrixFromNetworkX
    {
    private:
        /// <summary>
        /// the rewiring probability
        /// </summary>
        double p;
    protected:
        virtual std::string GetNetworkXGeneratorMethodCall();
        virtual std::string GetNetworkXAdjMatrix();
    public:
        /// <summary>
        /// constructor - creates the temporary file with the graph's adj matrix
        /// </summary>
        /// <param name="numOfElems">number of elements in the graph</param>
        /// <param name="numOfEdgesForNewElem">initial number of connections that a newly attached element will have</param>
        AdjacencyMatrixRandomGraph(int, double, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&);//(int numOfElems, double pOfEdgeCreation, bool isDirected)
        ~AdjacencyMatrixRandomGraph();
    };

    class AdjacencyMatrixWattsStrogatz : public AdjacencyMatrixFromNetworkX
    {
    private:
        /// <summary>
        /// the initial number of neighbours of each element
        /// </summary>
        int k;
        /// <summary>
        /// the rewiring probability
        /// </summary>
        double p;
    protected:
        virtual std::string GetNetworkXGeneratorMethodCall();
        virtual std::string GetNetworkXAdjMatrix();
    public:
        /// <summary>
        /// constructor - creates the temporary file with the graph's adj matrix
        /// </summary>
        /// <param name="numOfElems">number of elements in the graph</param>
        /// <param name="numOfEdgesForNewElem">initial number of connections that a newly attached element will have</param>
        AdjacencyMatrixWattsStrogatz(int, int, double, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&);//(int numOfElems, int numOfNeighbours, double rewiringProb, bool isDirected)
        ~AdjacencyMatrixWattsStrogatz();
    };

    class AdjacencyMatrixWattsStrogatzConn : public AdjacencyMatrixFromNetworkX
    {
    private:
        /// <summary>
        /// the initial number of neighbours of each element
        /// </summary>
        int k;
        /// <summary>
        /// the rewiring probability
        /// </summary>
        double p;
    protected:
        virtual std::string GetNetworkXGeneratorMethodCall();
        virtual std::string GetNetworkXAdjMatrix();
    public:
        /// <summary>
        /// constructor - creates the temporary file with the graph's adj matrix
        /// </summary>
        /// <param name="numOfElems">number of elements in the graph</param>
        /// <param name="numOfEdgesForNewElem">initial number of connections that a newly attached element will have</param>
        AdjacencyMatrixWattsStrogatzConn(int, int, double, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&);//(int numOfElems, int numOfNeighbours, double rewiringProb, bool isDirected)
        ~AdjacencyMatrixWattsStrogatzConn();
    };

    class AdjacencyMatrixCompleteGraph : public AdjacencyMatrixFromNetworkX
    {
    protected:
        virtual std::string GetNetworkXGeneratorMethodCall();
        virtual std::string GetNetworkXAdjMatrix();
    public:
        /// <summary>
        /// constructor - creates the temporary file with the graph's adj matrix
        /// </summary>
        /// <param name="numOfElems">number of elements in the graph</param>
        /// <param name="numOfEdgesForNewElem">initial number of connections that a newly attached element will have</param>
        AdjacencyMatrixCompleteGraph(int, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&);//(int numOfElems, bool isDirected)
        ~AdjacencyMatrixCompleteGraph();
    };

    class AdjacencyMatrixGridGraph : public AdjacencyMatrixFromNetworkX
    {
    private:
        /// <summary>
        /// Array containing the amount of elements on each dimension; dim of the graph == dim of this array
        /// </summary>
        std::vector<int> L;
        /// <summary>
        /// is it a periodic graph
        /// </summary>
        bool isPeriodic;
    protected:
        virtual std::string GetNetworkXGeneratorMethodCall();
        virtual std::string GetNetworkXAdjMatrix();
    public:
        /// <summary>
        /// constructor - creates the temporary file with the graph's adj matrix
        /// </summary>
        /// <param name="numOfElems">number of elements in the graph</param>
        /// <param name="numOfEdgesForNewElem">initial number of connections that a newly attached element will have</param>
        AdjacencyMatrixGridGraph(std::vector<int>, bool, bool, AdjacencyMatrixWeightType, double, std::mt19937_64&);//(std::vector<int> L, bool isPeriodic, bool isDirected)
        ~AdjacencyMatrixGridGraph();
    };
}

#endif
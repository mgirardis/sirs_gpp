#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <memory>
#include <math.h>
#include <vector>
#include <random>
#include "Distribution_c.h"
#include <stdexcept>

namespace Distribution
{
    std::unique_ptr<IDistribution> GetDistribution(DistributionParam_t par, std::mt19937_64& rng)
    {
        if (par.type == DistributionType::Uniform)
        {
            return std::unique_ptr<UniformDistribution>(new UniformDistribution(par.x0, par.x1, rng));
        }
        else if (par.type == DistributionType::Parabolic)
        {
            return std::unique_ptr<ParabolicDistribution>(new ParabolicDistribution(par.d, par.x0, rng));
        }
        else if (par.type == DistributionType::Gaussian)
        {
            return std::unique_ptr<GaussianDistribution>(new GaussianDistribution(par.sigma, par.mu, rng));
        }
        else if (par.type == DistributionType::Lorentzian)
        {
            return std::unique_ptr<LorentzianDistribution>(new LorentzianDistribution(par.gamma, rng));
        }
        else if (par.type == DistributionType::DiscreteGaussian)
        {
            return std::unique_ptr<PiecewiseGaussianDistribution>(new PiecewiseGaussianDistribution(par.sigma, par.mu, par.x0GaussPW, par.m, rng));
        }
        else if (par.type == DistributionType::DiscreteLorentzian)
        {
            return std::unique_ptr<PiecewiseLorentzianDistribution>(new PiecewiseLorentzianDistribution(par.gamma, par.x0LorentzPW, par.x1LorentzPW, par.m, rng));
        }
        else
        {
            throw std::invalid_argument("Unrecognized DistributionType");
        }
    }
    
    DistributionParam_t GetDefaultDistributionParam(DistributionType type)
    {
        DistributionParam_t par = {};
        par.type = type;
        if (type == DistributionType::Uniform)
        {
            par.x0 = 0.0;
            par.x1 = 1.0;
        }
        else if (type == DistributionType::Parabolic)
        {
            par.d = 10.0;
            par.x0 = 2.0;
        }
        else if (type == DistributionType::Gaussian)
        {
            par.sigma = 10.0;
            par.mu = 50.0;
        }
        else if (type == DistributionType::Lorentzian)
        {
            par.gamma = 5.0;
        }
        else if (type == DistributionType::DiscreteGaussian)
        {
            par.sigma = 10.0;
            par.mu = 50.0;
            par.x0GaussPW = 20.0;
            par.m = 3;
        }
        else if (type == DistributionType::DiscreteLorentzian)
        {
            par.gamma = 5.0;
            par.x0LorentzPW = -160.0;
            par.x1LorentzPW = -15.0;
            par.m = 3;
        }
        else
        {
            throw std::invalid_argument("Unrecognized DistributionType");
        }
        return par;
    }

    /*void WriteVector(std::string fileName, std::vector<double>& v)
    {
        FILE* f = fopen(fileName.c_str(), "w");
        int i = 0;
        size_t n = v.size();
        while (i < n)
        {
            fprintf(f, "%17.8E\n", v[i]);
            i++;
        }
        fclose(f);
    }*/

    UniformDistribution::UniformDistribution(double x0, double x1, std::mt19937_64& rng)
        : Distribution_c(rng)
    {
        this->x0 = x0;
        this->x1 = x1;
        this->distribution = std::unique_ptr<std::uniform_real_distribution<double>>(new std::uniform_real_distribution<double>(x0, x1));
    }
    UniformDistribution::~UniformDistribution(){}

    ParabolicDistribution::ParabolicDistribution(double d, double x0, std::mt19937_64& rng)
        : Distribution_c(rng)
    {
        this->d = std::abs(d);
        this->x0 = std::abs(x0);
        this->distribution = std::unique_ptr<beta_distribution<double>>(new beta_distribution<double>(2.0, 2.0));
    }
    ParabolicDistribution::~ParabolicDistribution(){}
    double ParabolicDistribution::NextNumber()
    {
        double x;
        do
        {
            x = d * (2.0 * this->distribution->operator()(this->randomNumGen) - 1.0);
        } while (std::abs(x) > x0);
        return x;
    }

    GaussianDistribution::GaussianDistribution(double sigma, double mu, std::mt19937_64& rng)
        : Distribution_c(rng)
    {
        this->sigma = sigma;
        this->mu = mu;
        this->distribution = std::unique_ptr<std::normal_distribution<double>>(new std::normal_distribution<double>(mu, sigma));
    }
    GaussianDistribution::~GaussianDistribution(){}

    LorentzianDistribution::LorentzianDistribution(double gamma, std::mt19937_64& rng)
        : Distribution_c(rng)
    {
        this->gamma = gamma;
        this->distribution = std::unique_ptr<std::cauchy_distribution<double>>(new std::cauchy_distribution<double>(0.0, gamma / 2.0));
    }
    LorentzianDistribution::~LorentzianDistribution(){}

    PiecewiseGaussianDistribution::PiecewiseGaussianDistribution(double sigma, double mu, double x0, int m, std::mt19937_64& rng)
        : PiecewiseDistribution(m, rng)
    {
        this->sigma = sigma;
        this->mu = mu;
        this->x0 = x0;
        this->CalcIntervals();
        this->CalcWeights();
        /*WriteVector("gIntervals.dat", intervals);
        WriteVector("gWeights.dat", weights);*/
        this->distribution = std::unique_ptr<std::piecewise_constant_distribution<double>>(new std::piecewise_constant_distribution<double>(intervals.begin(), intervals.end(), weights.begin()));
    }
    PiecewiseGaussianDistribution::~PiecewiseGaussianDistribution(){}
    double PiecewiseGaussianDistribution::DistFunc(double x)
    {
        return exp(-(x - mu)*(x - mu) / (2.0 * sigma * sigma)) / (sqrt(2.0 * M_PI) * sigma);
    }
    void PiecewiseGaussianDistribution::CalcIntervals()
    {
        intervals.clear();
        int n = 2 * m + 1;
        double dO = 2 * (mu - x0) / n; // interval is from -3 to 3
        int i = 0;
        while (i <= n)
        {
            intervals.push_back(x0 + i * dO);
            i++;
        }
    }

    PiecewiseLorentzianDistribution::PiecewiseLorentzianDistribution(double gamma, double x0, double x1, int m, std::mt19937_64& rng)
        : PiecewiseDistribution(m, rng)
    {
        this->gamma = gamma;
        this->x0 = x0;
        this->x1 = x1;
        this->CalcIntervals();
        this->CalcWeights();
        /*WriteVector("lIntervals.dat", intervals);
        WriteVector("lWeights.dat", weights);*/
        this->distribution = std::unique_ptr<std::piecewise_constant_distribution<double>>(new std::piecewise_constant_distribution<double>(intervals.begin(), intervals.end(), weights.begin()));
    }
    PiecewiseLorentzianDistribution::~PiecewiseLorentzianDistribution(){}
    double PiecewiseLorentzianDistribution::DistFunc(double x)
    {
        return gamma / (2.0 * M_PI * (gamma * gamma / 4.0 + x*x));
    }
    void PiecewiseLorentzianDistribution::CalcIntervals()
    {
        double x2 = -x1;
        intervals.clear();
        int n = 2 * m + 1;
        double dO;
        dO = (-x0 + x1) / m; // interval is from -32 to -3
        int i = 0;
        while (i < m)
        {
            intervals.push_back(x0 + i * dO);
            i++;
        }
        dO = 2 * x2 / n;
        i = 0;
        while (i < n)
        {
            intervals.push_back(x1 + i * dO);
            i++;
        }
        dO = (-x0 + x1) / m;
        i = 0;
        while (i <= m)
        {
            intervals.push_back(x2 + i * dO);
            i++;
        }
    }

    PiecewiseDistribution::PiecewiseDistribution(int m, std::mt19937_64& rng)
        : Distribution_c(rng)
    {
        this->m = m;
        /*this->CalcIntervals();
        this->CalcWeights();*/
    }
    PiecewiseDistribution::~PiecewiseDistribution(){}
    void PiecewiseDistribution::CalcIntervals(){}
    void PiecewiseDistribution::CalcWeights()
    {
        weights.clear();
        int n = static_cast<int>(intervals.size()) - 1;
        double temp;
        double sum = 0.0;
        int i = 0;
        while (i < n)
        {
            temp = this->DistFunc((intervals[i + 1] + intervals[i]) / 2.0);
            weights.push_back(temp);
            sum += temp;
            i++;
        }
        i = 0;
        while (i < n)
        {
            weights[i] = weights[i] / sum;
            i++;
        }
    }
    double PiecewiseDistribution::DistFunc(double x){ return x;	}

    template<typename DistType>
    Distribution_c<DistType>::Distribution_c(std::mt19937_64& rng)
        : randomNumGen(rng)
    {
    }

    template<typename DistType>
    std::vector<double> Distribution_c<DistType>::GetSample(int N)
    {
        std::vector<double> v;
        int i = 0;
        while (i < N)
        {
            v.push_back(this->NextNumber());
            i++;
        }
        return v;
    }
    template<typename DistType>
    double Distribution_c<DistType>::NextNumber()
    {
        return this->distribution->operator()(this->randomNumGen);
    }

    template<typename DistType>
    Distribution_c<DistType>::~Distribution_c()
    {
        //delete this->distribution;
    }
}

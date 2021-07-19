#ifndef _DISTRIBUTION_C_H
#define _DISTRIBUTION_C_H
#include <memory>
#include <vector>
#include <random>
#include "beta_distribution.h"

namespace Distribution
{
    //void WriteVector(std::string, std::vector<double>&);

    enum class DistributionType
    {
        begin,
        Gaussian,
        DiscreteGaussian,
        Lorentzian,
        DiscreteLorentzian,
        Uniform,
        Parabolic,
        end
    };
    static const char* DistributionTypeChar[] =
    {
        "begin",
        "Gaussian",
        "DiscreteGaussian",
        "Lorentzian",
        "DiscreteLorentzian",
        "Uniform",
        "Parabolic",
        "end"
    };

    class IDistribution
    {
    public:
        virtual std::vector<double> GetSample(int) = 0;
        virtual double NextNumber(void) = 0;
    };

    struct DistributionParam_t
    {
        DistributionType type;
        double x0;
        double x1;
        double sigma;
        double mu;
        double gamma;
        double d;
        double x0GaussPW;
        double x0LorentzPW;
        double x1LorentzPW;
        int m;
    };

    std::unique_ptr<IDistribution> GetDistribution(DistributionParam_t, std::mt19937_64&);
    DistributionParam_t GetDefaultDistributionParam(DistributionType);

    template<typename DistType>
    class Distribution_c : public IDistribution
    {
    public:
        Distribution_c(std::mt19937_64&);
        ~Distribution_c();

        virtual std::vector<double> GetSample(int);
        virtual double NextNumber(void);

    protected:
        std::mt19937_64& randomNumGen;
        std::unique_ptr<DistType> distribution;
    };

    class UniformDistribution : public Distribution_c<std::uniform_real_distribution<double>>
    {
    public:
        UniformDistribution(double, double, std::mt19937_64&);
        ~UniformDistribution();
        //virtual double NextNumber(void);
        //virtual std::vector<double> GetSample(int);
    private:
        double x0;
        double x1;
        //std::uniform_real_distribution<double>* distribution;
    };

    class ParabolicDistribution : public Distribution_c<beta_distribution<double>>
    {
    public:
        ParabolicDistribution(double, double, std::mt19937_64&);
        ~ParabolicDistribution();
        //virtual std::vector<double> GetSample(int);
        virtual double NextNumber(void);
    private:
        double d;
        double x0;
        //beta_distribution<double>* distribution;
    };

    class GaussianDistribution : public Distribution_c<std::normal_distribution<double>>
    {
    public:
        GaussianDistribution(double, double, std::mt19937_64&);
        ~GaussianDistribution();
        //virtual double NextNumber(void);
        //virtual std::vector<double> GetSample(int);
    private:
        double sigma;
        double mu;
        //std::normal_distribution<double>* distribution;
    };

    class LorentzianDistribution : public Distribution_c<std::cauchy_distribution<double>>
    {
    public:
        LorentzianDistribution(double, std::mt19937_64&);
        ~LorentzianDistribution();
        //virtual double NextNumber(void);
        //virtual std::vector<double> GetSample(int);
    private:
        double gamma;
        //std::cauchy_distribution<double>* distribution;
    };

    class PiecewiseDistribution : public Distribution_c<std::piecewise_constant_distribution<double>>
    {
    public:
        PiecewiseDistribution(int, std::mt19937_64&);
        ~PiecewiseDistribution();
        //virtual double NextNumber(void);
        //virtual std::vector<double> GetSample(int);
        int m;
    protected:
        //std::piecewise_constant_distribution<double>* distribution;

        std::vector<double> intervals;
        std::vector<double> weights;

        void CalcWeights(void);
        virtual void CalcIntervals(void);
        virtual double DistFunc(double);
    };

    class PiecewiseGaussianDistribution : public PiecewiseDistribution
    {
    public:
        PiecewiseGaussianDistribution(double, double, double, int, std::mt19937_64&);
        ~PiecewiseGaussianDistribution();
        //virtual double NextNumber(void);
        //virtual std::vector<double> GetSample(int);
    private:
        double sigma;
        double mu;
        double x0;
    protected:
        void CalcIntervals(void) override;
        virtual double DistFunc(double) override;
    };

    class PiecewiseLorentzianDistribution : public PiecewiseDistribution
    {
    public:
        PiecewiseLorentzianDistribution(double, double, double, int, std::mt19937_64&);
        ~PiecewiseLorentzianDistribution();
        //virtual double NextNumber(void);
        //virtual std::vector<double> GetSample(int);
    private:
        double gamma;
        double x0;
        double x1;
    protected:
        void CalcIntervals(void) override;
        virtual double DistFunc(double) override;
    };
}

#endif
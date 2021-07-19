#ifndef _COUPLINGS_H
#define _COUPLINGS_H

#include "Elements.h"

namespace Elements
{
    class IElement;
}

namespace Couplings
{
    enum class CouplingType
    {
        begin,
        PulseCoupling,
        end
    };
    static const char* CouplingTypeChar[] =
    {
        "begin",
        "PulseCoupling",
        "end"
    };

    struct CouplingParam_t
    {
        CouplingType type;
        Elements::IElement& pre;
        Elements::IElement& post;
        double weight;
    };

    class ICoupling
    {
    protected:
        Elements::IElement& pre;
        Elements::IElement& post;
        double I;

    public:
        ICoupling(Elements::IElement&, Elements::IElement&);
        ~ICoupling(void);
        Elements::IElement& GetPreElement(void);
        Elements::IElement& GetPostElement(void);
        void Reset(void);
        virtual void Step(void) = 0;
        virtual double GetSignal(void) = 0;
        virtual double GetCouplingItensity(void) = 0;
    };

    ICoupling* GetCoupling(CouplingParam_t);
    //std::unique_ptr<ICoupling> GetCoupling(CouplingParam_t);

    class WeightedCoupling : public ICoupling
    {
    public:
        WeightedCoupling(Elements::IElement&, Elements::IElement&, double);
        ~WeightedCoupling(void);
        virtual void Step(void);
        virtual double GetSignal(void);
        virtual double GetCouplingItensity(void);
    protected:
        double w;
    };

    class PulseCoupling : public WeightedCoupling
    {
    public:
        PulseCoupling(Elements::IElement&, Elements::IElement&, double);
        ~PulseCoupling();
        virtual void Step(void);
        virtual double GetSignal(void);
    };
}

#endif
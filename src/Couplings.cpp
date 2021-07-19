#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <memory>
#include <stdexcept>
#include "Elements.h"
#include "Couplings.h"

namespace Couplings
{
    //std::unique_ptr<ICoupling> GetCoupling(CouplingParam_t par)
    ICoupling* GetCoupling(CouplingParam_t par)
    {
        if (par.type == CouplingType::PulseCoupling)
        {
            //return std::make_unique<PulseCoupling>(par.pre, par.post, par.weight);
            return new PulseCoupling(par.pre, par.post, par.weight);
        }
        else
        {
            throw std::invalid_argument("Unrecognized CouplingType");
        }
    }

    ICoupling::ICoupling(Elements::IElement& pre, Elements::IElement& post)
        : pre(pre), post(post)
    {
        this->Reset();
    }
    ICoupling::~ICoupling(){}
    Elements::IElement& ICoupling::GetPreElement()
    {
        return this->pre;
    }
    Elements::IElement& ICoupling::GetPostElement()
    {
        return this->post;
    }
    void ICoupling::Reset()
    {
        this->I = 0.0;
    }

    WeightedCoupling::WeightedCoupling(Elements::IElement& pre, Elements::IElement& post, double weight)
        : ICoupling(pre, post)
    {
        this->w = weight;
    }
    WeightedCoupling::~WeightedCoupling(){}
    void WeightedCoupling::Step()
    {
        return;
    }
    double WeightedCoupling::GetSignal()
    {
        return this->w;
    }
    double WeightedCoupling::GetCouplingItensity()
    {
        return this->w;
    }

    PulseCoupling::PulseCoupling(Elements::IElement& pre, Elements::IElement& post, double weight)
        : WeightedCoupling(pre, post, weight)
    {
    }
    PulseCoupling::~PulseCoupling(){}
    void PulseCoupling::Step()
    {
        this->I = this->w * this->pre.GetState();
    }
    double PulseCoupling::GetSignal()
    {
        return this->I;
    }
}
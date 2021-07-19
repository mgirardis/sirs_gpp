#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <vector>
#include "Elements.h"
#include <stdexcept>
#include <iostream>

namespace Elements
{
    //std::unique_ptr<IElement> GetElement(ElementParam_t par)
    IElement* GetElement(ElementParam_t par)
    {
        if (par.type == ElementType::SI)
        {
            //return std::make_unique<SIElement>(par.ind, par.threshold, par.s0);
            return new SIElement(par.ind, par.threshold, par.s0, par.normalizeInput);
        }
        else if (par.type == ElementType::SIRSDet)
        {
            //return std::make_unique<SIRSDetElement>(par.ind, par.threshold, par.s0, par.nRefStates);
            return new SIRSDetElement(par.ind, par.threshold, par.s0, par.nRefStates, par.normalizeInput);
        }
        throw std::invalid_argument("Unrecognized ElementType");
    }

    IElement::IElement(bool normalizeInput)
    {
        if (normalizeInput)
        {
            this->SetTotalNeighWeight = &IElement::SetTotalNeighWeight_Sum;
        }
        else
        {
            this->SetTotalNeighWeight = &IElement::SetTotalNeighWeight_One;
        }
    }
    IElement::~IElement()
    {
        /*for (size_t i = 0; i < this->inNeighbors.size(); i++)
        {
            delete this->inNeighbors[i];
        }//*/
    }
    int IElement::GetIndex()
    {
        return this->ind;
    }
    double IElement::GetState()
    {
        return this->s;
    }
    void IElement::AddNeighbor(Couplings::ICoupling& n)
    {
        this->inNeighbors.push_back(&n);
    }
    double IElement::GetNeighborsSignal()
    {
        size_t i = 0, n = this->inNeighbors.size();
        double sum = 0.0;
        while (i < n)
        {
            sum += this->inNeighbors[i]->GetSignal();
            i++;
        }
        return sum;
    }
    void IElement::SetTotalNeighWeight_Sum()
    {
        size_t i = 0, n = this->inNeighbors.size();
        double sum = 0.0;
        while (i < n)
        {
            sum += this->inNeighbors[i]->GetCouplingItensity();
            i++;
        }
        this->totalNeighWeight = sum;
    }
    void IElement::SetTotalNeighWeight_One()
    {
        this->totalNeighWeight = 1.0;
    }
    void IElement::SetIC(double s0)
    {
        this->s = s0;
    }

    ThresholdElement::ThresholdElement(int ind, double threshold, double s0, bool normalizeInput)
		: IElement(normalizeInput)
    {
        this->ind = ind;
        this->threshold = threshold;
        this->s = s0;
    }
    ThresholdElement::~ThresholdElement(){}
    void ThresholdElement::Step()
    {
        return;
    }
    void ThresholdElement::Step(double Iext)
    {
        return;
    }
    void ThresholdElement::SetThreshold(double theta)
    {
        this->threshold = theta;
        this->s = 0.0;
    }

    SIElement::SIElement(int ind, double threshold, double s0, bool normalizeInput)
        : ThresholdElement(ind, threshold, s0, normalizeInput)
    {}
    SIElement::~SIElement(){}
    void SIElement::Step()
    {
        if ((this->GetNeighborsSignal() / this->totalNeighWeight) >= this->threshold)
        {
            this->s = 1.0;
        }
    }
    void SIElement::Step(double Iext)
    {
        //std::cout << "node " << this->ind << " neighbors signal = " << this->GetNeighborsSignal() << std::endl;
        if ((Iext + (this->GetNeighborsSignal() / this->totalNeighWeight)) >= this->threshold)
        {
            this->s = 1.0;
        }
    }

	//                            (int    , double          , double   , int  , bool); // index, threshold, ic, number of ref states, normalizeInput
    SIRSDetElement::SIRSDetElement(int ind, double threshold, double s0, int m, bool normalizeInput)
        : ThresholdElement(ind, threshold, s0, normalizeInput)
    {
        this->m = m;
        this->sid = ( s0 > 0.0 ? 1 : 0 );
    }
    SIRSDetElement::~SIRSDetElement(){}
    int SIRSDetElement::GetStateId(void)
    {
        return this->sid;
    }
    void SIRSDetElement::Step()
    {
        this->Step(0.0);
    }
    void SIRSDetElement::Step(double Iext)
    {
        if (this->sid == 0)
        {
            if ((Iext + (this->GetNeighborsSignal() / this->totalNeighWeight)) >= this->threshold)
            {
                this->s = 1.0;
                this->sid = 1;
            }
        }
        else if (this->sid == 1)
        {
            this->s = 0.0;
            this->sid = 2;
        }
        else
        {
            if (this->sid == (m + 1))
            {
                this->sid = 0;
            }
            else
            {
                this->sid++;
            }
        }
    }

	// GLElement(int, double, double, double, double, double, double, bool)
	// index, V, Gamma, mu, VB, VR, threshold, normalizeInput
    GLElement::GLElement(int ind, double X0, double V0, double Gamma, double mu, double VB, double VR, double VT, bool normalizeInput, PhiFunctionType phiType, std::mt19937_64& rand)
        : ThresholdElement(ind, VT, X0, normalizeInput), rand(rand)
    {
		this->V = V0;
		this->Gamma = Gamma;
		this->mu = mu;
		this->VB = VB;
		this->VR = VR;
		if (phiType == PhiFunctionType::Linear)
		{
			this->Phi = &GLElement::Phi_Linear;
		}
		else if (phiType == PhiFunctionType::Logistic)
		{
			this->Phi = &GLElement::Phi_Logistic;
		}
	}
    GLElement::~GLElement(){}
	double GLElement::Phi_Linear(double V)
	{
		double tg = this->threshold + 1.0 / this->Gamma;
		return this->Gamma * (V - this->threshold) * (double)(V > this->threshold) * (double)(V < tg) + (double)(V >= tg);
	}
	double GLElement::Phi_Logistic(double V)
	{
		return this->Gamma * (V - this->threshold) * (double)(V > this->threshold) / (1.0 + this->Gamma * (V - this->threshold));
	}
    void GLElement::Step()
    {
		this->Step(0.0);
    }
    void GLElement::Step(double Iext)
    {
		// s is the X variable
		this->V = (this->mu * (this->V - this->VB) + this->VB + Iext + this->GetNeighborsSignal()/this->totalNeighWeight) * (1.0 - this->s) + this->s * this->VR;
		this->s = (double)(this->rand() < (this->*Phi)(V));
    }

    // GLElement(int, double, double, double, double, double, double, bool)
	// index, V, Gamma, mu, VB, VR, threshold, normalizeInput
    GLAdaptElement::GLAdaptElement(int ind, double X0, double V0, double Gamma, double mu, double VB, double VR, double VT, double tau, bool normalizeInput, PhiFunctionType phiType, std::mt19937_64& rand)
        : GLElement(ind, X0, V0, Gamma, mu, VB, VR, VT, normalizeInput, phiType, rand)
    {
		this->V = V0;
		this->Gamma = Gamma;
		this->mu = mu;
		this->VB = VB;
		this->VR = VR;
		this->invTau = 1.0 / tau;
	}
    GLAdaptElement::~GLAdaptElement(){}
    void GLAdaptElement::Step(double Iext)
    {
		// s is the X variable
		this->Gamma = (1.0 + this->invTau - this->s) * this->Gamma;
		this->V = (this->mu * (this->V - this->VB) + this->VB + Iext + this->GetNeighborsSignal()/this->totalNeighWeight) * (1.0 - this->s) + this->s * this->VR;
		this->s = (double)(this->rand() < (this->*Phi)(V));
    }
}
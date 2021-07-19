#ifndef _ELEMENTS_H
#define _ELEMENTS_H

#include <memory>
#include <random>
#include <vector>
#include "Couplings.h"

namespace Couplings
{
    class ICoupling;
}

namespace Elements
{
    /*template <class T> class VectorOfRefs : public std::vector<T *> {
    public:
        inline T & at(const uint64_t i) {
            T * x = std::vector<T *>::at(i);
            return *x;
        }
    };*/

	enum class ElementType
	{
		begin,
		SI,
		SIRSDet,
		GL,
		GLAdapt,
		end
	};
	static const char* ElementTypeChar[] =
	{
		"begin",
		"SI",
		"SIRSDet",
		"GL",
		"GLAdapt",
		"end"
	};

	enum class PhiFunctionType
	{
		begin,
		Logistic,
		Linear,
		end
	};
	static const char* PhiFunctionTypeChar[] =
	{
		"begin",
		"Logistic",
		"Linear",
		"end"
	};

    struct ElementParam_t
    {
        ElementType type;
        int ind;
        double threshold;
        double s0;
        int nRefStates;
		bool normalizeInput;
		PhiFunctionType phiType;
    };

    class IElement
    {
    protected:
		// type that points to void functions of IElement class
        int ind;
        double s;
        double totalNeighWeight;
        std::vector<Couplings::ICoupling*> inNeighbors;
        virtual double GetNeighborsSignal(void);
        void SetTotalNeighWeight_Sum();
        void SetTotalNeighWeight_One();

    public:
		typedef void(IElement::*SetTotalNeighWeightFn_ptr)();
        IElement(bool); // normalizeInput
        ~IElement();
        int GetIndex(void);
        double GetState(void);
        void AddNeighbor(Couplings::ICoupling&);
        virtual void SetIC(double);
        virtual void Step(void) = 0;
        virtual void Step(double) = 0;
        virtual void SetThreshold(double) = 0;
        SetTotalNeighWeightFn_ptr SetTotalNeighWeight;
    };

    //std::unique_ptr<IElement> GetElement(ElementParam_t);
    IElement* GetElement(ElementParam_t);

    class ThresholdElement : public IElement
    {
    protected:
        double threshold;

    public:
        ThresholdElement(int, double, double, bool); // index, threshold, ic, normalizeInput
        ~ThresholdElement();

        virtual void Step(void);
        virtual void Step(double);
        virtual void SetThreshold(double);
    };

    class SIElement : public ThresholdElement
    {
    public:
        SIElement(int, double, double, bool); // index, threshold, ic, normalizeInput
        ~SIElement();

        virtual void Step(void);
        virtual void Step(double);
    };

    class SIRSDetElement : public ThresholdElement
    {
    private:
        int m; // number of refractory states
        int sid; // state identifier

    public:
        SIRSDetElement(int, double, double, int, bool); // index, threshold, ic, number of ref states, normalizeInput
        ~SIRSDetElement();

        int GetStateId(void);
        virtual void Step(void);
        virtual void Step(double);
    };

    class GLElement : public ThresholdElement
    {
    protected:
        double V; // membrane potential used to calculate firing probability
        double Gamma; // Gain of the firing probability (parameter or variable)
        double mu; // Time constant for V integration (leakage)
        double VB; // Baseline potential (if no activity, V = VB)
        double VR; // Reset potential (V = VR if neuron spikes)

		double Phi_Logistic(double);
		double Phi_Linear(double);
		std::mt19937_64& rand;

    public:
        GLElement(int, double, double, double, double, double, double, double, bool, PhiFunctionType, std::mt19937_64&); // index, V, Gamma, mu, VB, VR, threshold, normalizeInput
        ~GLElement();

		typedef double(GLElement::* PhiFunc_ptr)(double);
        virtual void Step(void);
        virtual void Step(double);
        PhiFunc_ptr Phi;
    };

    class GLAdaptElement : public GLElement
    {
    private:
		double invTau;

    public:
        GLAdaptElement(int, double, double, double, double, double, double, double, double, bool, PhiFunctionType, std::mt19937_64&); // index, V, Gamma, mu, VB, VR, threshold, tau, normalizeInput
        ~GLAdaptElement();

        //virtual void Step(void);
        virtual void Step(double);
    };
}

#endif
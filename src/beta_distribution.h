#ifndef _BETA_DISTRIBUTION_H
#define _BETA_DISTRIBUTION_H

#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <functional>

namespace Distribution {

	//template <typename URNG, typename RealType = double>
	template <typename RealType = double>
	class beta_distribution
	{
	public:
		typedef RealType result_type;

		class param_type
		{
		public:
			typedef beta_distribution distribution_type;

			explicit param_type(RealType a = 2.0, RealType b = 2.0)
				: a_param(a), b_param(b) { }

			RealType a() const { return a_param; }
			RealType b() const { return b_param; }

			bool operator==(const param_type& other) const
			{
				return (a_param == other.a_param &&
					b_param == other.b_param);
			}

			bool operator!=(const param_type& other) const
			{
				return !(*this == other);
			}

		private:
			RealType a_param, b_param;
		};

		/*explicit beta_distribution(RealType d, RealType r, RealType a = 2.0, RealType b = 2.0)
			: a_gamma(a), b_gamma(b)
		{
			this->d = d;
			this->r = r;
			this->generate = std::bind(&beta_distribution::generateInterval, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		}*/

		explicit beta_distribution(RealType a = 2.0, RealType b = 2.0)
			: a_gamma(a), b_gamma(b)
		{
			//this->generate = std::bind(&beta_distribution::generateStd, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		}
		explicit beta_distribution(const param_type& param)
			: a_gamma(param.a()), b_gamma(param.b())
		{
			//this->generate = std::bind(&beta_distribution::generateStd, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		}

		void reset() { }

		param_type param() const
		{
			return param_type(a(), b());
		}

		void param(const param_type& param)
		{
			a_gamma = gamma_dist_type(param.a());
			b_gamma = gamma_dist_type(param.b());
		}

		template <typename URNG>
		result_type operator()(URNG& engine)
		{
			return generate(engine, a_gamma, b_gamma);
		}

		template <typename URNG>
		result_type operator()(URNG& engine, const param_type& param)
		{
			gamma_dist_type a_param_gamma(param.a()),
				b_param_gamma(param.b());
			return generate(engine, a_param_gamma, b_param_gamma);
		}

		result_type min() const { return this->r1; }
		result_type max() const { return this->r2; }

		RealType a() const { return a_gamma.alpha(); }
		RealType b() const { return b_gamma.alpha(); }

		bool operator==(const beta_distribution<result_type>& other) const
		{
			return (param() == other.param() &&
				a_gamma == other.a_gamma &&
				b_gamma == other.b_gamma);
		}

		bool operator!=(const beta_distribution<result_type>& other) const
		{
			return !(*this == other);
		}

	private:
		/*RealType r;
		RealType d;*/

		typedef std::gamma_distribution<result_type> gamma_dist_type;

		gamma_dist_type a_gamma, b_gamma;

		//std::function<result_type(URNG&, gamma_dist_type&, gamma_dist_type&)> generate;

		/*template <typename URNG>
		result_type generateInterval(URNG& engine,
			gamma_dist_type& x_gamma,
			gamma_dist_type& y_gamma)
		{
			result_type x;
			do
			{
				x = x_gamma(engine);
				x = this->d * (2 * x / (x + y_gamma(engine)) - 1);
			} while (abs(x) > this->r);
			return x;
		}*/

		template <typename URNG>
		result_type generate(URNG& engine,
			gamma_dist_type& x_gamma,
			gamma_dist_type& y_gamma)
		{
			result_type x = x_gamma(engine);
			return x / (x + y_gamma(engine));
		}
	};

	template <typename CharT, typename RealType>
	std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os,
		const beta_distribution<RealType>& beta)
	{
		os << "~Beta(" << beta.a() << "," << beta.b() << ")";
		return os;
	}

	template <typename CharT, typename RealType>
	std::basic_istream<CharT>& operator>>(std::basic_istream<CharT>& is,
		beta_distribution<RealType>& beta)
	{
		std::string str;
		RealType a, b;
		if (std::getline(is, str, '(') && str == "~Beta" &&
			is >> a && is.get() == ',' && is >> b && is.get() == ')') {
			beta = beta_distribution<RealType>(a, b);
		}
		else {
			is.setstate(std::ios::failbit);
		}
		return is;
	}

}

#endif
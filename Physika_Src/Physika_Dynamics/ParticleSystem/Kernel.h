#pragma once
#include "Physika_Core/Utilities/cuda_utilities.h"

namespace Physika {

	template<typename Real>
	class Kernel
	{
	public:
		__host__ __device__ Kernel() {};
		__host__ __device__ ~Kernel() {};

		__host__ __device__ inline virtual Real Weight(const Real r, const Real h)
		{
			return Real(0);
		}

		__host__ __device__ inline virtual Real Gradient(const Real r, const Real h)
		{
			return Real(0);
		}
	};

	//spiky kernel
	template<typename Real>
	class SpikyKernel : public Kernel<Real>
	{
	public:
		__host__ __device__ SpikyKernel() : Kernel<Real>() {};
		__host__ __device__ ~SpikyKernel() {};

		__host__ __device__ inline Real Weight(const Real r, const Real h) override
		{
			const Real q = r / h;
			if (q > 1.0f) return 0.0f;
			else {
				const Real d = 1.0 - q;
				const Real hh = h*h;
				return 15.0f / ((Real)M_PI * hh * h) * d * d * d;
			}
		}

		__host__ __device__ inline Real Gradient(const Real r, const Real h) override
		{
			const Real q = r / h;
			if (q > 1.0f) return 0.0;
			//else if (r==0.0f) return 0.0f;
			else {
				const Real d = 1.0 - q;
				const Real hh = h*h;
				return -45.0f / ((Real)M_PI * hh*h) *d*d;
			}
		}
	};

	//cubic kernel
	template<typename Real>
	class CubicKernel : public Kernel<Real>
	{
	public:
		__host__ __device__ CubicKernel() : Kernel<Real>() {};
		__host__ __device__ ~CubicKernel() {};

		__host__ __device__ inline Real Weight(const Real r, const Real h) override
		{
			const Real hh = h*h;
			const Real q = 2.0f*r / h;

			const Real alpha = 3.0f / (2.0f * (Real)M_PI * hh * h);

			if (q > 2.0f) return 0.0f;
			else if (q >= 1.0f)
			{
				//1/6*(2-q)*(2-q)*(2-q)
				const Real d = 2.0f - q;
				return alpha / 6.0f*d*d*d;
			}
			else
			{
				//(2/3)-q*q+0.5f*q*q*q
				const Real qq = q*q;
				const Real qqq = qq*q;
				return alpha*(2.0f / 3.0f - qq + 0.5f*qqq);
			}
		}

		__host__ __device__ inline Real Gradient(const Real r, const Real h) override
		{
			const Real hh = h*h;
			const Real q = 2.0f*r / h;

			const Real alpha = 3.0f / (2.0f * (Real)M_PI * hh * h);

			if (q > 2.0f) return 0.0f;
			else if (q >= 1.0f)
			{
				//-0.5*(2.0-q)*(2.0-q)
				const Real d = 2.0f - q;
				return -0.5f*alpha*d*d;
			}
			else
			{
				//-2q+1.5*q*q
				const Real qq = q*q;
				return alpha*(-2.0f*q + 1.5f*qq);
				//return alpha*(-0.5);
			}
		}
	};

	template<typename Real>
	class SmoothKernel : public Kernel<Real>
	{
	public:
		__host__ __device__ SmoothKernel() : Kernel<Real>() {};
		__host__ __device__ ~SmoothKernel() {};

		__host__ __device__ inline Real Weight(const Real r, const Real h) override
		{
			const Real q = r / h;
			if (q > 1.0f) return 0.0f;
			else {
				return 1.0f - q*q;
			}
		}

		__host__ __device__ inline Real Gradient(const Real r, const Real h) override
		{
			const Real q = r / h;
			if (q > 1.0f) return 0.0f;
			else {
				const Real hh = h*h;
				const Real dd = 1 - q*q;
				const Real alpha = 1.0f;// (Real) 945.0f / (32.0f * (Real)M_PI * hh *h);
				return -alpha * dd;
			}
		}
	};

	template<typename Real>
	class QuarticKernel : public Kernel<Real>
	{
	public:
		__host__ __device__ QuarticKernel() : Kernel<Real>() {};
		__host__ __device__ ~QuarticKernel() {};

		__host__ __device__ inline Real Weight(const Real r, const Real h) override
		{
			const Real hh = h*h;
			const Real q = 2.5f*r / h;
			if (q > 2.5) return 0.0f;
			else if (q > 1.5f)
			{
				const Real d = 2.5f - q;
				const Real dd = d*d;
				return 0.0255f*dd*dd / hh;
			}
			else if (q > 0.5f)
			{
				const Real d = 2.5f - q;
				const Real t = 1.5f - q;
				const Real dd = d*d;
				const Real tt = t*t;
				return 0.0255f*(dd*dd - 5.0f*tt*tt) / hh;
			}
			else
			{
				const Real d = 2.5f - q;
				const Real t = 1.5f - q;
				const Real w = 0.5f - q;
				const Real dd = d*d;
				const Real tt = t*t;
				const Real ww = w*w;
				return 0.0255f*(dd*dd - 5.0f*tt*tt + 10.0f*ww*ww) / hh;
			}
		}

		__host__ __device__ inline Real Gradient(const Real r, const Real h) override
		{
			const Real hh = h*h;
			const Real q = 2.5f*r / h;
			if (q > 2.5) return 0.0f;
			else if (q > 1.5f)
			{
				//0.102*(2.5-q)^3
				const Real d = 2.5f - q;
				return -0.102f*d*d*d / hh;
			}
			else if (q > 0.5f)
			{
				const Real d = 2.5f - q;
				const Real t = 1.5f - q;
				return -0.102f*(d*d*d - 5.0f*t*t*t) / hh;
			}
			else
			{
				const Real d = 2.5f - q;
				const Real t = 1.5f - q;
				const Real w = 0.5f - q;
				return -0.102f*(d*d*d - 5.0f*t*t*t + 10.0f*w*w*w) / hh;
			}
		}
	};
}

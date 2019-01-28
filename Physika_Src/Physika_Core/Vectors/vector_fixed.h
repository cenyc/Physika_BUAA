#pragma once

#include <cstddef>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "Physika_Core/Utilities/cuda_utilities.h"

namespace Physika
{
	template<class T, unsigned int N>
	class FixedVector
	{
	public:
		// type definitions
		typedef T              value_type;
		typedef T*             iterator;
		typedef const T*       const_iterator;
		typedef T&             reference;
		typedef const T&       const_reference;
		typedef unsigned int   size_type;

		CPU_GPU_FUNC_DECL FixedVector()
		{
		}

		/// Specific constructor for 1-element vectors.
		template<int NN = N, typename std::enable_if<NN == 1, int>::type = 0>
		CPU_GPU_FUNC_DECL explicit FixedVector(value_type r1)
		{
			static_assert(N == 1, "");
			this->data[0] = r1;
		}

		/// Specific constructor for 2-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 2, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2)
		{
			static_assert(N == 2, "");
			this->data[0] = r1;
			this->data[1] = r2;
		}

		/// Specific constructor for 3-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 3, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2, value_type r3)
		{
			static_assert(N == 3, "");
			this->data[0] = r1;
			this->data[1] = r2;
			this->data[2] = r3;
		}

		/// Specific constructor for 4-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 4, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2, value_type r3, value_type r4)
		{
			static_assert(N == 4, "");
			this->data[0] = r1;
			this->data[1] = r2;
			this->data[2] = r3;
			this->data[3] = r4;
		}

		/// Specific constructor for 5-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 5, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2, value_type r3, value_type r4, value_type r5)
		{
			static_assert(N == 5, "");
			this->data[0] = r1;
			this->data[1] = r2;
			this->data[2] = r3;
			this->data[3] = r4;
			this->data[4] = r5;
		}

		/// Specific constructor for 6-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 6, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2, value_type r3, value_type r4, value_type r5, value_type r6)
		{
			static_assert(N == 6, "");
			this->data[0] = r1;
			this->data[1] = r2;
			this->data[2] = r3;
			this->data[3] = r4;
			this->data[4] = r5;
			this->data[5] = r6;
		}

		/// Specific constructor for 7-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 7, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2, value_type r3, value_type r4, value_type r5, value_type r6, value_type r7)
		{
			static_assert(N == 7, "");
			this->data[0] = r1;
			this->data[1] = r2;
			this->data[2] = r3;
			this->data[3] = r4;
			this->data[4] = r5;
			this->data[5] = r6;
			this->data[6] = r7;
		}

		/// Specific constructor for 8-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 8, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2, value_type r3, value_type r4, value_type r5, value_type r6, value_type r7, value_type r8)
		{
			static_assert(N == 8, "");
			this->data[0] = r1;
			this->data[1] = r2;
			this->data[2] = r3;
			this->data[3] = r4;
			this->data[4] = r5;
			this->data[5] = r6;
			this->data[6] = r7;
			this->data[7] = r8;
		}

		/// Specific constructor for 9-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 9, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2, value_type r3, value_type r4, value_type r5, value_type r6, value_type r7, value_type r8, value_type r9)
		{
			static_assert(N == 9, "");
			this->data[0] = r1;
			this->data[1] = r2;
			this->data[2] = r3;
			this->data[3] = r4;
			this->data[4] = r5;
			this->data[5] = r6;
			this->data[6] = r7;
			this->data[7] = r8;
			this->data[8] = r9;
		}

		/// Specific constructor for 10-elements vectors.
		template<int NN = N, typename std::enable_if<NN == 10, int>::type = 0>
		CPU_GPU_FUNC_DECL FixedVector(value_type r1, value_type r2, value_type r3, value_type r4, value_type r5, value_type r6, value_type r7, value_type r8, value_type r9, value_type r10)
		{
			static_assert(N == 10, "");
			this->data[0] = r1;
			this->data[1] = r2;
			this->data[2] = r3;
			this->data[3] = r4;
			this->data[4] = r5;
			this->data[5] = r6;
			this->data[6] = r7;
			this->data[7] = r8;
			this->data[8] = r9;
			this->data[9] = r10;
		}


		// iterator support
		CPU_GPU_FUNC_DECL iterator begin()
		{
			return data;
		}
		CPU_GPU_FUNC_DECL const_iterator begin() const
		{
			return data;
		}
		CPU_GPU_FUNC_DECL iterator end()
		{
			return data + N;
		}
		CPU_GPU_FUNC_DECL const_iterator end() const
		{
			return data + N;
		}

		// operator[]
		CPU_GPU_FUNC_DECL reference operator[](size_type i)
		{
#ifndef NDEBUG
			assert(i < N && "index in fixed_array must be smaller than size");
#endif
			return data[i];
		}
		CPU_GPU_FUNC_DECL const_reference operator[](size_type i) const
		{
#ifndef NDEBUG
			assert(i < N && "index in fixed_array must be smaller than size");
#endif
			return data[i];
		}

		// front() and back()
		CPU_GPU_FUNC_DECL reference front()
		{
			return data[0];
		}
		CPU_GPU_FUNC_DECL const_reference front() const
		{
			return data[0];
		}
		CPU_GPU_FUNC_DECL reference back()
		{
			return data[N - 1];
		}
		CPU_GPU_FUNC_DECL const_reference back() const
		{
			return data[N - 1];
		}

		// size is constant
		CPU_GPU_FUNC_DECL static size_type size()
		{
			return N;
		}
		CPU_GPU_FUNC_DECL static bool empty()
		{
			return N == 0;
		}
		CPU_GPU_FUNC_DECL static size_type max_size()
		{
			return N;
		}

		// direct access to data
		CPU_GPU_FUNC_DECL const T* getDataPtr() const
		{
			return data;
		}

		// assignment with type conversion
		template <typename T2>
		CPU_GPU_FUNC_DECL FixedVector<T, N>& operator= (const FixedVector<T2, N>& rhs)
		{
			//std::copy(rhs.begin(),rhs.end(), begin());
			for (size_type i = 0; i < N; i++)
				data[i] = rhs[i];
			return *this;
		}

		// assign one value to all elements
		CPU_GPU_FUNC_DECL inline void assign(const T& value)
		{
			//std::fill_n(begin(),size(),value);
			for (size_type i = 0; i < N; i++)
				data[i] = value;
		}

		//template<int NN = N, typename std::enable_if<NN>0,int>::type = 0>
		CPU_FUNC_DECL inline friend std::ostream& operator << (std::ostream& out, const FixedVector<T, N>& a)
		{
			static_assert(N > 0, "Cannot create a zero size arrays");
			for (size_type i = 0; i < N - 1; i++)
				out << a.elems[i] << " ";
			out << a.elems[N - 1];
			return out;
		}

		CPU_FUNC_DECL inline friend std::istream& operator >> (std::istream& in, FixedVector<T, N>& a)
		{
			for (size_type i = 0; i < N; i++)
				in >> a.elems[i];
			return in;
		}

		CPU_GPU_FUNC_DECL inline bool operator < (const FixedVector& v) const
		{
			for (size_type i = 0; i < N; i++)
			{
				if (data[i] < v[i])
					return true;  // (*this)<v
				else if (data[i] > v[i])
					return false; // (*this)>v
			}
			return false; // (*this)==v
		}

		public:
			T data[N];

	};

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 2> makeFixedVector(const T& v0, const T& v1)
	{
		FixedVector<T, 2> v;
		v[0] = v0;
		v[1] = v1;
		return v;
	}

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 3> makeFixedVector(const T& v0, const T& v1, const T& v2)
	{
		FixedVector<T, 3> v;
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		return v;
	}

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 4> makeFixedVector(const T& v0, const T& v1, const T& v2, const T& v3)
	{
		FixedVector<T, 4> v;
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
		return v;
	}

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 5> makeFixedVector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4)
	{
		FixedVector<T, 5> v;
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
		v[4] = v4;
		return v;
	}

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 6> makeFixedVector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5)
	{
		FixedVector<T, 6> v;
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
		v[4] = v4;
		v[5] = v5;
		return v;
	}

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 7> makeFixedVector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6)
	{
		FixedVector<T, 7> v;
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
		v[4] = v4;
		v[5] = v5;
		v[6] = v6;
		return v;
	}

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 8> makeFixedVector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7)
	{
		FixedVector<T, 8> v;
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
		v[4] = v4;
		v[5] = v5;
		v[6] = v6;
		v[7] = v7;
		return v;
	}

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 9> makeFixedVector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7, const T& v8)
	{
		FixedVector<T, 9> v;
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
		v[4] = v4;
		v[5] = v5;
		v[6] = v6;
		v[7] = v7;
		v[8] = v8;
		return v;
	}

	template<class T>
	CPU_GPU_FUNC_DECL inline FixedVector<T, 10> makeFixedVector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7, const T& v8, const T& v9)
	{
		FixedVector<T, 10> v;
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
		v[4] = v4;
		v[5] = v5;
		v[6] = v6;
		v[7] = v7;
		v[8] = v8;
		v[9] = v9;
		return v;
	}
} // namespace helper

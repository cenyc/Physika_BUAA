#include "Function1Pt.h"
#include "Physika_Core/Utilities/cuda_helper_math.h"
#include "Physika_Core/Utilities/cuda_utilities.h"
namespace Physika
{
	namespace Function1Pt
	{
		template<typename T1, typename T2>
		__global__ void KerLength(T1* lhs, T2* rhs, int num)
		{
			int pId = threadIdx.x + (blockIdx.x * blockDim.x);
			if (pId >= num) return;

			lhs[pId] = length(rhs[pId]);
		}

		template<typename T1, typename T2>
		void Length(DeviceArray<T1>& lhs, DeviceArray<T2>& rhs)
		{
			assert(lhs.size() == rhs.size());
			unsigned pDim = cudaGridSize(rhs.size(), BLOCK_SIZE);
			KerLength << <pDim, BLOCK_SIZE >> > (lhs.getDataPtr(), rhs.getDataPtr(), lhs.size());
		}
	}
}
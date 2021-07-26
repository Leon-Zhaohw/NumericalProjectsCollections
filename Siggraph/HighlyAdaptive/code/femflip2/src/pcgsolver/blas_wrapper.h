#ifndef BLAS_WRAPPER_H
#define BLAS_WRAPPER_H

// Simple placeholder code for BLAS calls - replace with calls to a real BLAS library

#include <vector>

namespace BLAS{
	
	// dot products ==============================================================
	template<class T>
	inline T dot(const std::vector<T> &x, const std::vector<T> &y)
	{ 
		//return cblas_ddot((int)x.size(), &x[0], 1, &y[0], 1); 
		
		T sum = 0;
		for(int i = 0; i < x.size(); ++i)
			sum += x[i]*y[i];
		return sum;
	}
	
	// inf-norm (maximum absolute value: index of max returned) ==================
	
	template<class T>
	inline int index_abs_max(const std::vector<T> &x)
	{ 
		//return cblas_idamax((int)x.size(), &x[0], 1); 
		int maxind = 0;
		T maxvalue = 0;
		for(int i = 0; i < x.size(); ++i) {
			if(fabs(x[i]) > maxvalue) {
				maxvalue = fabs(x[i]);
				maxind = i;
			}
		}
		return maxind;
	}
	
	// inf-norm (maximum absolute value) =========================================
	// technically not part of BLAS, but useful
	
	template<class T>
	inline T abs_max(const std::vector<T> &x)
	{ return std::fabs(x[index_abs_max(x)]); }
	
	// saxpy (y=alpha*x+y) =======================================================
	
	template<class T>
	inline void add_scaled(T alpha, const std::vector<T> &x, std::vector<T> &y)
	{ 
		//cblas_daxpy((int)x.size(), alpha, &x[0], 1, &y[0], 1); 
		for(int i = 0; i < x.size(); ++i)
			y[i] += alpha*x[i];
	}
	
}
#endif

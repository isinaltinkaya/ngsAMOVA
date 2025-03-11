#ifndef __EIGEN_H__
#define __EIGEN_H__
#include <cstddef>  // size_t
#include <Eigen/Eigenvalues>

namespace Eigen {
	template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
	class Matrix;
}  // namespace Eigen


using EigenMatrixRowMajorD = Eigen::Matrix<
	double,
	Eigen::Dynamic,
	Eigen::Dynamic,
	Eigen::RowMajor
>;

using EigenVectorD = Eigen::Matrix<
	double,
	Eigen::Dynamic,
	1,
	Eigen::ColMajor
>;

#endif
#include "Triangulation.h"
#include <iostream>

Triangulation::Triangulation(
			const Eigen::MatrixXd& P0,
			const Eigen::MatrixXd& P1,
			const Eigen::VectorXd& x0,
			const Eigen::VectorXd& x1):
			P(P0, P1),
			x(x0, x1)
{
}



Triangulation::Triangulation(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& _P,
							const std::pair<Eigen::VectorXd, Eigen::VectorXd>& _x):
							P(_P),
							x(_x)
{
}




Eigen::Vector3d Triangulation::solve() const
{
	return Triangulation::solve(P, x);
}




Eigen::Vector3d Triangulation::solve(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& P,
									const std::pair<Eigen::VectorXd, Eigen::VectorXd>& x)
{
	Eigen::MatrixXd A = buildMatrixA(P, x);
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinV);
	Eigen::MatrixXd h = svd.matrixV().rightCols(1);
	h /= h(3, 0);
	Eigen::Vector3d X;
	X << h(0, 0), h(1, 0), h(2, 0);
	return X;
}




Eigen::MatrixXd Triangulation::buildMatrixA(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& P,
											const std::pair<Eigen::VectorXd, Eigen::VectorXd>& x)
{
	Eigen::MatrixXd A(4, 4);

	A.row(0) = x.first.x()  * P.first.row(2)  - P.first.row(0);
	A.row(1) = x.first.y()  * P.first.row(2)  - P.first.row(1);
	A.row(2) = x.second.x() * P.second.row(2) - P.second.row(0);
	A.row(3) = x.second.y() * P.second.row(2) - P.second.row(1);

	return A;
}


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




Eigen::Vector3d Triangulation::solve(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& _P,
									const std::pair<Eigen::VectorXd, Eigen::VectorXd>& _x)
{
	Eigen::MatrixXd A = buildMatrixA(_P, _x);
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinV);
	Eigen::MatrixXd h = svd.matrixV().rightCols(1);
	h /= h(3, 0);
	Eigen::Vector3d X;
	X << h(0, 0), h(1, 0), h(2, 0);
	return X;
}




Eigen::MatrixXd Triangulation::buildMatrixA(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& _P,
											const std::pair<Eigen::VectorXd, Eigen::VectorXd>& _x)
{
	Eigen::MatrixXd A(4, 4);

	Eigen::VectorXd v0 = _x.first;
	Eigen::VectorXd v1 = _x.second;

	Eigen::VectorXd p0r0 = _P.first.row(0);
	Eigen::VectorXd p0r1 = _P.first.row(1);
	Eigen::VectorXd p0r2 = _P.first.row(2);

	Eigen::VectorXd p1r0 = _P.second.row(0);
	Eigen::VectorXd p1r1 = _P.second.row(1);
	Eigen::VectorXd p1r2 = _P.second.row(2);

	Eigen::VectorXd v0xTimesP0R2MinusP0R0 = v0[0] * p0r2 - p0r0;
	v0xTimesP0R2MinusP0R0.normalize();

	Eigen::VectorXd v0yTimesP0R2MinusP0R1 = v0[1] * p0r2 - p0r1;
	v0yTimesP0R2MinusP0R1.normalize();

	Eigen::VectorXd v1xTimesP0R2MinusP0R0 = v1[0] * p1r2 - p1r0;
	v1xTimesP0R2MinusP0R0.normalize();

	Eigen::VectorXd v1yTimesP0R2MinusP0R1 = v1[1] * p1r2 - p1r1;
	v1yTimesP0R2MinusP0R1.normalize();

	A << v0xTimesP0R2MinusP0R0[0], v0xTimesP0R2MinusP0R0[1], v0xTimesP0R2MinusP0R0[2], v0xTimesP0R2MinusP0R0[3],
		v0yTimesP0R2MinusP0R1[0], v0yTimesP0R2MinusP0R1[1], v0yTimesP0R2MinusP0R1[2], v0yTimesP0R2MinusP0R1[3],
		v1xTimesP0R2MinusP0R0[0], v1xTimesP0R2MinusP0R0[1], v1xTimesP0R2MinusP0R0[2], v1xTimesP0R2MinusP0R0[3],
		v1yTimesP0R2MinusP0R1[0], v1yTimesP0R2MinusP0R1[1], v1yTimesP0R2MinusP0R1[2], v1yTimesP0R2MinusP0R1[3];

	return A;
}


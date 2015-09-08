#include "Reconstruction3D.h"
#include "Triangulation.h"
#include <iostream>
#include <limits>
#include <iomanip>



Eigen::MatrixXd Reconstruction3D::buildMatrixA(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts)
{
	int i = 0;
	Eigen::MatrixXd A(pts.size(), 9);		// 8 x 9, 2n x 9

	for (const auto p : pts)
	{
		A(i, 0) = p.second.x() * p.first.x();
		A(i, 1) = p.second.x() * p.first.y();
		A(i, 2) = p.second.x();
		A(i, 3) = p.second.y() * p.first.x();
		A(i, 4) = p.second.y() * p.first.y();
		A(i, 5) = p.second.y();
		A(i, 6) = p.first.x();
		A(i, 7) = p.first.y();
		A(i, 8) = 1;
		++i;
	}

	return A;
}



Eigen::MatrixXd Reconstruction3D::applyConstraint(const Eigen::MatrixXd& inputMat)
{
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner> svd(inputMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::VectorXd singularValues = svd.singularValues();
	double singularValue = (singularValues[0] + singularValues[1]) * 0.5;
	Eigen::DiagonalMatrix< double, 3, 3 > diagonal(singularValues(0), singularValues(1), 0.0);
	Eigen::MatrixXd D = diagonal.toDenseMatrix();
	Eigen::MatrixXd constrainedMat = svd.matrixU() * D * svd.matrixV().transpose();

	//std::cout << std::fixed << "In  Mat:" << std::endl << inputMat << std::endl << std::endl;
	//std::cout << std::fixed << "Out Mat:" << std::endl << constrainedMat << std::endl << std::endl;

	std::cout << std::fixed
		<< "In  determinant : " << inputMat.determinant() << std::endl
		<< "Out determinant : " << constrainedMat.determinant() << std::endl << std::endl;

	return constrainedMat;
}


Eigen::MatrixXd Reconstruction3D::computeF(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts)
{
	Eigen::MatrixXd A = Reconstruction3D::buildMatrixA(pts);

	//std::cout << std::fixed << "A:" << std::endl << A << std::endl << std::endl;

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
	Eigen::MatrixXd kernel = svd.matrixV().col(svd.matrixV().cols() - 1);

	Eigen::MatrixXd F(3, 3);
	F(0, 0) = kernel(0);
	F(0, 1) = kernel(1);
	F(0, 2) = kernel(2);
	F(1, 0) = kernel(3);
	F(1, 1) = kernel(4);
	F(1, 2) = kernel(5);
	F(2, 0) = kernel(6);
	F(2, 1) = kernel(7);
	F(2, 2) = kernel(8);

	return F;
}





Eigen::MatrixXd Reconstruction3D::computeE(const Eigen::MatrixXd& K, const Eigen::MatrixXd& F)
{
	Eigen::MatrixXd E = K.transpose() * F * K;
	E = applyConstraint(E);
	E /= E(2, 2);
	return E;
}




Eigen::MatrixXd Reconstruction3D::computeK(const Eigen::MatrixXd& F)
{
	Eigen::Matrix3d R, Q;
	//Reconstruction3D::RQdecomposition(P, R, Q);
	Reconstruction3D::QRdecomposition(F, R, Q);
	return R;
}





Eigen::Matrix3d Reconstruction3D::computeEpipoleMat(const Eigen::Matrix3d& F)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(F.transpose(), Eigen::ComputeFullV);
	Eigen::MatrixXd kernel = svd.matrixV().col(svd.matrixV().cols() - 1);

	Eigen::Matrix3d e(3, 3);
	e(0, 0) = kernel(0);
	e(0, 1) = kernel(1);
	e(0, 2) = kernel(2);
	e(1, 0) = kernel(3);
	e(1, 1) = kernel(4);
	e(1, 2) = kernel(5);
	e(2, 0) = kernel(6);
	e(2, 1) = kernel(7);
	e(2, 2) = kernel(8);

	return e;
}


Eigen::Vector3d Reconstruction3D::computeEpipole(const Eigen::Matrix3d& F)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(F.transpose(), Eigen::ComputeFullV);
	Eigen::MatrixXd kernel = svd.matrixV().col(svd.matrixV().cols() - 1);
	return Eigen::Vector3d( -kernel(5), kernel(2), -kernel(1));
}


Eigen::MatrixXd Reconstruction3D::computeP(const Eigen::MatrixXd& F, const Eigen::MatrixXd& eMat)
{
	Eigen::MatrixXd Fe(3, 4);
	Fe.block(0, 0, 3, 3) = F;
	Fe(0, 3) = -eMat(1, 2);
	Fe(1, 3) =  eMat(0, 2);
	Fe(2, 3) = -eMat(0, 1);

	return eMat * Fe;
}


Eigen::MatrixXd Reconstruction3D::computeP(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
											const Eigen::MatrixXd& E)
{
	Eigen::JacobiSVD< Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner > svd(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::MatrixXd Vt = svd.matrixV().transpose();


	Eigen::MatrixXd W(3, 3);
	W << 0.0, -1.0, 0.0,
		1.0, 0.0, 0.0,
		0.0, 0.0, 1.0;

	Eigen::VectorXd u3 = U.col(2);

	Eigen::MatrixXd P0(Eigen::MatrixXd::Identity(3, 4)); // Just K0: no translation or rotation.

	Eigen::MatrixXd P1noT = U * W * Vt; // P without translation part.

	Eigen::MatrixXd P1(3, 4);
	P1.block(0, 0, 3, 3) = P1noT;
	P1.block(0, 3, 3, 1) = u3;


	Eigen::MatrixXd bestP1;
	int numCorrectSolutions = 0;
	if (checkP(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
	}

	P1.block(0, 3, 3, 1) = -u3;


	if (checkP(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
	}

	P1noT = U * W.transpose() * Vt;
	P1.block(0, 0, 3, 3) = P1noT;
	P1.block(0, 3, 3, 1) = u3;

	if (checkP(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
	}

	P1.block(0, 3, 3, 1) = -u3;

	if (checkP(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
	}

	if (numCorrectSolutions < 1)
		std::cerr << "[Error]   : None of the results for P1 was accepted!" << std::endl << std::endl;
	else if (numCorrectSolutions > 1)
		std::cerr << "[Error]   : More than one solution found" << std::endl << std::endl;
	else
		std::cerr << "[Info]    : Solution found for P1: \n" << bestP1 << std::endl << std::endl;

	return bestP1;
}


bool Reconstruction3D::checkP(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
	const Eigen::MatrixXd& P0,
	const Eigen::MatrixXd& P1)
{
	Eigen::VectorXd p0 = pts[0].first;
	Eigen::VectorXd p1 = pts[0].second;

	std::vector< std::pair< Eigen::VectorXd, Eigen::VectorXd > > pair(1);
	pair[0] = std::pair< Eigen::VectorXd, Eigen::VectorXd >(p0, p1);


	Triangulation triangulation(std::make_pair(P0, P1), pts[0]);
	Eigen::Vector3d X = triangulation.solve();

	Eigen::VectorXd x0 = P0 * X;
	Eigen::VectorXd x1 = P1 * X;

	return (x0.z() > 0.0) && (x1.z() > 0.0);
}



double Reconstruction3D::computeError(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const Eigen::MatrixXd& F)
{
	double error = 0;

	for (const auto x : pts)
	{
		Eigen::Vector3d x0 = x.first.homogeneous();
		Eigen::Vector3d x1 = x.second.homogeneous();

		error += x1.transpose() * F * x0;
	}

	return error;
}


void Reconstruction3D::solveCS(float &c, float &s, float a, float b)
{
	float den = (a * a) + (b * b);
	c = -b / std::sqrt(den);
	s = a / std::sqrt(den);
}

void Reconstruction3D::RQdecomposition(Eigen::MatrixXd A, Eigen::Matrix3d &R, Eigen::Matrix3d &Q)
{
	float c, s;
	Eigen::Matrix3d Qx, Qy, Qz;
	Eigen::Matrix3d M;
	M << A(0, 0), A(0, 1), A(0, 2),
		A(1, 0), A(1, 1), A(1, 2),
		A(2, 0), A(2, 1), A(2, 2);

	R = M;
	Reconstruction3D::solveCS(c, s, R(2, 1), R(2, 2));
	Qx << 1, 0, 0,
		0, c, -s,
		0, s, c;
	R *= Qx;

	Reconstruction3D::solveCS(c, s, R(2, 0), -R(2, 2));
	Qy << c, 0, s,
		0, 1, 0,
		-s, 0, c;
	R *= Qy;

	Reconstruction3D::solveCS(c, s, R(1, 0), R(1, 1));
	Qz << c, -s, 0,
		s, c, 0,
		0, 0, 1;
	R *= Qz;

	if (std::abs(R(1, 0)) > 0.00005 || std::abs(R(2, 0)) > 0.00005 || std::abs(R(2, 1)) > 0.00005)
		std::cerr << "PROBLEM WITH RQdecomposition" << std::endl;;

	//    if(R(1,0)!= 0) R(1,0) = 0;
	//    if(R(2,0)!= 0) R(2,0) = 0;
	//    if(R(2,1)!= 0) R(2,1) = 0;


	Q = Qz.transpose()  * Qy.transpose() * Qx.transpose();

	if (R(0, 0) < 0)
	{
		R.col(0) *= -1;
		Q.row(0) *= -1;
	}
	if (R(1, 1) < 0)
	{
		R.col(1) *= -1;
		Q.row(1) *= -1;
	}
	if (R(2, 2) < 0)
	{
		R.col(2) *= -1;
		Q.row(2) *= -1;
	}
}

void Reconstruction3D::QRdecomposition(Eigen::MatrixXd A, Eigen::Matrix3d &R, Eigen::Matrix3d &Q)
{
	Eigen::Matrix3d m;
	m = A.topLeftCorner<3, 3>();
	//    Eigen::FullPivHouseholderQR<Eigen::Matrix3d> qr(m.rows(), m.cols());
	//    Eigen::HouseholderQR<Eigen::Matrix3d> qr(m);
	Eigen::ColPivHouseholderQR<Eigen::Matrix3d> qr(m);
	qr.solve(m);
	R = qr.householderQ().inverse() * qr.matrixQR();
	Q = qr.householderQ();
	//    std::cout << qr.matrixQ().inverse() * qr.matrixQR() << std::endl;
	
	std::cout
		<< std::endl << std::fixed
		<< "R: " << std::endl
		<< R << std::endl << std::endl;

	std::cout
		<< std::endl << std::fixed
		<< "Q: " << std::endl
		<< Q << std::endl << std::endl;
}
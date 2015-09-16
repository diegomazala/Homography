#include "Reconstruction3D.h"
#include "Triangulation.h"
#include "RansacDLT.h"
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
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(inputMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner> svd(inputMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::VectorXd singularValues = svd.singularValues();
	double singularValue = (singularValues[0] + singularValues[1]) * 0.5;
	Eigen::DiagonalMatrix< double, 3, 3 > diagonal(singularValue, singularValue, 0.0);
	Eigen::MatrixXd D = diagonal.toDenseMatrix();
	Eigen::MatrixXd constrainedMat = svd.matrixU() * D * svd.matrixV().transpose();

	//std::cout << std::fixed << "In  Mat:" << std::endl << inputMat << std::endl << std::endl;
	//std::cout << std::fixed << "Out Mat:" << std::endl << constrainedMat << std::endl << std::endl;

	//std::cout << std::fixed
	//	<< "F In  determinant : " << inputMat.determinant() << std::endl
	//	<< "F Out determinant : " << constrainedMat.determinant() << std::endl << std::endl;


	//MatrixXd DMat = D.toDenseMatrix();
	//MatrixXd restricted = svd.matrixU() * DMat * svd.matrixV().transpose();

	//std::cout << "========== E Restricion Appliance =============: " << endl
	//	<< "Singular diagonal: " << endl << DMat << endl << endl
	//	<< "U:" << endl << svd.matrixU() << endl << endl
	//	<< "V:" << endl << svd.matrixV() << endl << endl
	//	<< "E: " << endl << restricted << endl << endl
	//	<< "Determinant: " << restricted.determinant() << endl << endl
	//	<< "========== Restriction Appliance End =============" << endl << endl;

	//*m_resultH = restricted;

	return constrainedMat;
}


Eigen::MatrixXd Reconstruction3D::computeF(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts)
{
	Eigen::MatrixXd A = Reconstruction3D::buildMatrixA(pts);

	//std::cout << std::fixed << "F A:" << std::endl << A << std::endl << std::endl;

	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
	Eigen::MatrixXd kernel = svd.matrixV().col(svd.matrixV().cols() - 1);

	//std::cout << std::fixed << "kernel:" << std::endl << kernel << std::endl << std::endl;

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





Eigen::MatrixXd Reconstruction3D::computeE(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& K, const Eigen::MatrixXd& F)
{
	Eigen::MatrixXd E = K.second.transpose() * F * K.first;
	//E /= E(2, 2);	// o erro aumenta
	//Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner> svd(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
	//Eigen::VectorXd singularValues = svd.singularValues();
	//Eigen::DiagonalMatrix< double, 3, 3 > diagonal(1, 1, 0);
	//Eigen::MatrixXd D = diagonal.toDenseMatrix();
	//Eigen::MatrixXd constrainedMat = svd.matrixU() * D * svd.matrixV().transpose();
	//E = constrainedMat;

	//std::cout << std::fixed
	//	<< "E In  determinant : " << F.determinant() << std::endl
	//	<< "E Out determinant : " << constrainedMat.determinant() << std::endl << std::endl;

	
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


Eigen::MatrixXd Reconstruction3D::computeP(const Eigen::MatrixXd& F)
{
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner> svd(F.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(F.transpose(), Eigen::ComputeFullV);
	Eigen::MatrixXd kernel = svd.matrixV().col(svd.matrixV().cols() - 1);

	Eigen::Matrix3d eMat(3, 3);
	eMat(0, 0) = kernel(0);
	eMat(0, 1) = kernel(1);
	eMat(0, 2) = kernel(2);
	eMat(1, 0) = kernel(3);
	eMat(1, 1) = kernel(4);
	eMat(1, 2) = kernel(5);
	eMat(2, 0) = kernel(6);
	eMat(2, 1) = kernel(7);
	eMat(2, 2) = kernel(8);

	

	Eigen::Vector3d e(-kernel(5), kernel(2), -kernel(1));

	Eigen::MatrixXd Fe(3, 4);
	Fe.block(0, 0, 3, 3) = F;
	Fe.col(3) = e;
	//Fe(0, 3) = -e(1, 2);
	//Fe(1, 3) =  e(0, 2);
	//Fe(2, 3) = -e(0, 1);

	Eigen::MatrixXd P = eMat * Fe;

	std::cout
		<< std::fixed << std::endl
		<< "[Info]  Epipole             : " << std::endl << e.transpose() << std::endl << std::endl
		<< "[Info]  Epipole Skew Matrix : " << std::endl << eMat << std::endl << std::endl
		<< "[Info]  F : " << std::endl << F << std::endl << std::endl
		<< "[Info]  Fe : " << std::endl << Fe << std::endl << std::endl
		<< "[Info]  P'                  : " << std::endl << P << std::endl << std::endl;
		

	return P;
}

Eigen::MatrixXd Reconstruction3D::computeP(const Eigen::MatrixXd& F, const Eigen::MatrixXd& eMat)
{
	Eigen::MatrixXd Fe(3, 4);
	Fe.block(0, 0, 3, 3) = F;
	Fe(0, 3) = -eMat(1, 2);
	Fe(1, 3) = eMat(0, 2);
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

	//std::cout << "--- Computing P Matrix: " << std::endl
	//	<< "E: " << std::endl << E << std::endl << std::endl
	//	<< "U: " << std::endl << U << std::endl << std::endl
	//	<< "D: " << std::endl << svd.singularValues() << std::endl << std::endl
	//	<< "Vt: " << std::endl << Vt << std::endl << std::endl
	//	<< "W: " << std::endl << W << std::endl << std::endl
	//	<< "u3 :" << u3 << std::endl << std::endl;

	Eigen::MatrixXd P0(Eigen::MatrixXd::Identity(3, 4)); // Just K0: no translation or rotation.

	Eigen::MatrixXd P1noT = U * W * Vt; // P without translation part.

	Eigen::MatrixXd P1(3, 4);
	P1.block(0, 0, 3, 3) = P1noT;
	P1.block(0, 3, 3, 1) = u3;

	Eigen::MatrixXd bestP1 = P0;
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

	std::cout << "Number of correct solutions: " << numCorrectSolutions << std::endl << std::endl;

	if (numCorrectSolutions < 1)
		std::cerr << "[Error]   : None of the results for P1 was accepted!" << std::endl << std::endl;
	else if (numCorrectSolutions > 1)
		std::cerr << "[Error]   : More than one solution found : " << numCorrectSolutions << std::endl << std::endl;
//	else
//		std::cerr << "[Info]    : Solution found for P: \n" << bestP1 << std::endl << std::endl;

	return bestP1;
}


bool Reconstruction3D::checkP(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
	const Eigen::MatrixXd& P0,
	const Eigen::MatrixXd& P1)
{
	Eigen::VectorXd X = Triangulation::solve(std::make_pair(P0, P1), pts[0]);

	Eigen::VectorXd x0 = P0 * X;
	Eigen::VectorXd x1 = P1 * X;

	//std::cout << std::fixed	<< "P' determinant : " << P1.block(0,0,3,3).determinant() << std::endl << std::endl;


	//std::cout << "========== P1 Check =============: " << std::endl
	//	<< "P0:" << std::endl << P0 << std::endl << std::endl
	//	<< "P1:" << std::endl << P1 << std::endl << std::endl
	//	<< "x: " << std::endl << p0.transpose() << std::endl << std::endl
	//	<< "x': " << std::endl << p1.transpose() << std::endl << std::endl
	//	<< "3d point: " << std::endl << X.transpose() << std::endl << std::endl
	//	<< "Reprojected x: " << std::endl << x0.transpose() << std::endl << std::endl
	//	<< "Reprojected x': " << std::endl << x1.transpose() << std::endl << std::endl
	//	<< "========== P1 Check End =============: " << std::endl << std::endl;


	return (x0.z() > 0.0) && (x1.z() > 0.0);
}

void Reconstruction3D::computeP(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
								const Eigen::MatrixXd& E, 
								std::vector<Eigen::MatrixXd>& P_solutions)
{
	P_solutions.clear();

	Eigen::JacobiSVD< Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner > svd(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::MatrixXd Vt = svd.matrixV().transpose();


	Eigen::MatrixXd W(3, 3);
	W << 0.0, -1.0, 0.0,
		1.0, 0.0, 0.0,
		0.0, 0.0, 1.0;

	Eigen::VectorXd u3 = U.col(2);

	//std::cout << "--- Computing P Matrix: " << std::endl
	//	<< "E: " << std::endl << E << std::endl << std::endl
	//	<< "U: " << std::endl << U << std::endl << std::endl
	//	<< "D: " << std::endl << svd.singularValues() << std::endl << std::endl
	//	<< "Vt: " << std::endl << Vt << std::endl << std::endl
	//	<< "W: " << std::endl << W << std::endl << std::endl
	//	<< "u3 :" << u3 << std::endl << std::endl;

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

		std::cout << std::endl << "[Info]  Best solution for P: 1st" << std::endl << std::endl;
	}
	P_solutions.push_back(P1);

	P1.block(0, 3, 3, 1) = -u3;

	if (checkP(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
		std::cout << std::endl << "[Info]  Best solution for P: 2nd" << std::endl << std::endl;
	}
	P_solutions.push_back(P1);

	P1noT = U * W.transpose() * Vt;
	P1.block(0, 0, 3, 3) = P1noT;
	P1.block(0, 3, 3, 1) = u3;


	if (checkP(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
		std::cout << std::endl << "[Info]  Best solution for P: 3rd" << std::endl << std::endl;
	}
	P_solutions.push_back(P1);

	P1.block(0, 3, 3, 1) = -u3;


	if (checkP(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
		std::cout << std::endl << "[Info]  Best solution for P: 4th" << std::endl << std::endl;
	}
	P_solutions.push_back(P1);

	std::cout << "[Info]  Number of correct solutions: " << numCorrectSolutions << std::endl << std::endl;
}


Eigen::MatrixXd Reconstruction3D::selectBestP(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D,
												const std::pair<Eigen::Matrix3d, Eigen::Matrix3d>& K,
												const std::vector<Eigen::MatrixXd>& P_solutions)
{
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd>	P(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
	
	P.first = Eigen::MatrixXd::Identity(3, 4);
	P.first.block(0, 0, 3, 3) = K.first;

	
	int bestChoice = 0;
	double minError = FLT_MAX;

	int solution = 0;
	for (auto m : P_solutions)
	{
		double error = 0;
		//P.second = m;
		P.second = K.second * m;

		P.second /= P.second(2, 2);

		for (int i = 0; i < points2D.size(); ++i)
		{
			const auto& x = points2D[i];
			Eigen::Vector3d x0 = x.first.homogeneous();
			Eigen::Vector3d x1 = x.second.homogeneous();

			const Eigen::VectorXd& X = Triangulation::solve(P, x);

			Eigen::Vector3d xx0 = P.first * X;
			Eigen::Vector3d xx1 = P.second * X;

			xx0 /= xx0[2];
			xx1 /= xx1[2];

			//std::cout
			//	<< std::fixed << i << " : " << x0.transpose() << '\t' << x1.transpose() << std::endl
			//	<< std::fixed << i << " : " << xx0.transpose() << '\t' << xx1.transpose() << std::endl << std::endl;

			double d0 = (x0 - xx0).norm();
			double d1 = (x1 - xx1).norm();

			//error += d0 * d0 + d1 * d1;
			error += d0 + d1;
		}

		if (error < minError)
			bestChoice = solution;

		++solution;
		//std::cout << std::endl << "Error: " << error << std::endl << std::endl;
	}

	return P_solutions[bestChoice];
}


double Reconstruction3D::pointLineDistance(Eigen::Vector2d point, Eigen::Vector3d line)
{
	const double x = point[0];
	const double y = point[1];
	const double a = line[0];
	const double b = line[1];
	const double c = line[2];

	const double num = std::abs(a * x + b * y + c);
	const double den = std::sqrt(a * a + b * b);

	return num / den;
}


double Reconstruction3D::computeError(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const Eigen::MatrixXd& F)
{
	double error = 0;

	for (const auto x : pts)
	{
		Eigen::Vector3d x0 = x.first.homogeneous();
		Eigen::Vector3d x1 = x.second.homogeneous();

		Eigen::Vector3d l0 = F.transpose() * x1;
		Eigen::Vector3d l1 = F * x0;

		double d0 = pointLineDistance(x.first, l0);
		double d1 = pointLineDistance(x.second, l1);
		double d = d0 * d0 + d1 * d1;

		//std::cout
		//	<< std::fixed
		//	<< x0.transpose() << '\t' << x1.transpose() << std::endl
		//	<< l0.transpose() << '\t' << l1.transpose() << std::endl
		//	<< "distance: " << d << '\t' << d0 << '\t' << d1 << std::endl << std::endl;

		error += d;
	}

	return error;
}

#if 0
double Reconstruction3D::computeGeometricError(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& P, double outlierThreshold, int& inliers)
{
	double error = 0;
	inliers = 0;
	
	for (const auto x : pts)
	{
		Eigen::Vector3d x0 = x.first.homogeneous();
		Eigen::Vector3d x1 = x.second.homogeneous();

		Eigen::VectorXd X = Triangulation::solve(P, x);

		Eigen::Vector3d xx0 = P.first * X;
		Eigen::Vector3d xx1 = P.second * X;

		xx0 /= xx0[2];
		xx1 /= xx1[2];

		double d0 = Eigen::Vector3d(x0 - xx0).norm();
		double d1 = Eigen::Vector3d(x1 - xx1).norm();

		//error += d0 * d0 + d1 * d1;
		double localError = d0 + d1;

		std::cout << std::fixed 
			<< "x0, x1      : " << x0.transpose() << '\t' << x1.transpose() << std::endl
			<< "xx0, xx1    : " << xx0.transpose() << '\t' << xx1.transpose() << std::endl
			<< "Local Error : " << localError << "  (" << d0 << ", " << d1 << ")" << std::endl;

		if (localError < outlierThreshold)
			++inliers;

		error += localError;
	}

	std::cout << "Inliers     : " << inliers << std::endl;

	return error;
}
#endif

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



ReconstructionDLT Reconstruction3D::solve(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
											const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& K, 
											int maxIterations)
{
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd>	P(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));

	std::vector<ReconstructionDLT> dltArray;

	for (int it = 0; it < maxIterations; ++it)
	{
		std::vector<int> indices = RansacDLT::random4Indices(0, (int)pts.size() - 1);

		// copy points from original array to the 4-array to be used for reconstruction
		std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> points2D;
		for each (auto i in indices)
			points2D.push_back(pts[i]);


		dltArray.push_back(ReconstructionDLT(points2D, K));
		ReconstructionDLT& dlt = dltArray.back();

		dlt.solve();

		std::cout << it << " : " << dlt.error << " , " << dlt.inliersCount << std::endl;
	}


	std::sort(dltArray.begin(), dltArray.end());


	return dltArray.front();
}





ReconstructionDLT::ReconstructionDLT(
		const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
		const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& _K):
		points2D(pts),
		K(_K), 
		inliersCount(0)
{

}

void ReconstructionDLT::solve(double outlierThreshold)
{
	P = std::make_pair(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));

	//
	// Normalize points
	// 
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d> T = DLT::normalizePoints(points2D, points2DNorm);

	//
	// Compute F matrix
	//
	Eigen::MatrixXd Fn = Reconstruction3D::computeF(points2DNorm);
	Fn = Reconstruction3D::applyConstraint(Fn);
	Eigen::MatrixXd F = DLT::denormalizeH(Fn, T);

	//
	// Compute E matrix
	//
	Eigen::MatrixXd E = Reconstruction3D::computeE(K, F);

	//std::vector<Eigen::MatrixXd> P_solutions;
	//Reconstruction3D::computeP(points2DNorm, E, P_solutions);

	P.first = Eigen::MatrixXd::Identity(3, 4);
	//P.first.block(0, 0, 3, 3) = K.first;

	//Eigen::MatrixXd P2 = Reconstruction3D::selectBestP(points2D, K, P_solutions);
	Eigen::MatrixXd P2 = Reconstruction3D::computeP(points2DNorm, E);
	//P.second = K.second * P2;
	P.second = P2;

	//error = Reconstruction3D::computeGeometricError(points2DNorm, P);
	//std::cout << "Geometric Error: " << Reconstruction3D::computeGeometricError(points2DNorm, P) << std::endl;

	inliersCount = 0;
	error = 0;
	for (const auto x : points2DNorm)
	{
		Eigen::Vector3d x0 = x.first.homogeneous();
		Eigen::Vector3d x1 = x.second.homogeneous();

		Eigen::VectorXd X = Triangulation::solve(P, x);

		Eigen::Vector3d xx0 = P.first * X;
		Eigen::Vector3d xx1 = P.second * X;

		xx0 /= xx0[2];
		xx1 /= xx1[2];

		double d0 = Eigen::Vector3d(x0 - xx0).norm();
		double d1 = Eigen::Vector3d(x1 - xx1).norm();

		//std::cout << std::fixed << "x0, x1   : " << x0.transpose() << '\t' << x1.transpose() << std::endl;
		//std::cout << std::fixed << "xx0, xx1 : " << xx0.transpose() << '\t' << xx1.transpose() << std::endl;
		//std::cout << "d0, d1 : " << d0 << ", " << d1 << std::endl;

		//error += d0 * d0 + d1 * d1;
		double localError = d0 + d1;

		if (localError < outlierThreshold)
			++inliersCount;

		error += localError;
	}
	//std::cout << "Sum Error: " << error << std::endl << std::endl;

	P.first.block(0, 0, 3, 3) = K.first;
	P.second = K.second * P2;
}


bool ReconstructionDLT::operator < (ReconstructionDLT const &other)
{
	//	return (this->inliers > other.inliers) && (this->error.first + this->error.second) < (other.error.first + other.error.second);
	//	return (this->error.first + this->error.second) < (other.error.first + other.error.second);
	return (this->inliersCount) > (other.inliersCount) && (this->error) < (other.error);
}
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



Eigen::MatrixXd Reconstruction3D::applyConstraint(const Eigen::MatrixXd& F)
{
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::VectorXd singularValues = svd.singularValues();
	double singularValue = (singularValues[0] + singularValues[1]) * 0.5;
	Eigen::DiagonalMatrix< double, 3, 3 > diagonal(singularValue, singularValue, 0.0);
	//Eigen::DiagonalMatrix< double, 3, 3 > diagonal(singularValues(0), singularValues(1), 0.0);
	Eigen::MatrixXd D = diagonal.toDenseMatrix();
	Eigen::MatrixXd constrainedMat = svd.matrixU() * D * svd.matrixV().transpose();

	//std::cout << std::fixed << "U:" << std::endl << svd.matrixU() << std::endl << std::endl;
	//std::cout << std::fixed << "D:" << std::endl << D << std::endl << std::endl;
	//std::cout << std::fixed << "V:" << std::endl << svd.matrixV() << std::endl << std::endl;


	//std::cout << std::fixed << "Constrained Mat:" << std::endl << constrainedMat << std::endl << std::endl;

	//std::cout << std::fixed 
	//	<< "F restricted determinant: " << constrainedMat.determinant() << std::endl
	//	<< "F determinant           : " << F.determinant() << std::endl << std::endl;

	return constrainedMat;
}


Eigen::MatrixXd Reconstruction3D::computeF(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts)
{
	Eigen::MatrixXd A = Reconstruction3D::buildMatrixA(pts);

	std::cout << std::fixed << "A:" << std::endl << A << std::endl << std::endl;

	//Eigen::MatrixXd kernel = A.fullPivLu().kernel();
	Eigen::JacobiSVD<Eigen::MatrixXd> SVD(A, Eigen::ComputeFullV);
	Eigen::MatrixXd kernel = SVD.matrixV().col(SVD.matrixV().cols() - 1);

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
	//E /= E(2, 2);
	return E;
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


bool Reconstruction3D::checkP(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
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




Eigen::VectorXd Reconstruction3D::triangulation(const Eigen::MatrixXd& P0, const Eigen::MatrixXd& P1, const Eigen::VectorXd& p0, const Eigen::VectorXd& p1)
{
	return Eigen::Vector3d();
}


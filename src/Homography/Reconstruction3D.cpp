#include "Reconstruction3D.h"
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



Eigen::MatrixXd Reconstruction3D::applyRestriction(const Eigen::MatrixXd& F)
{
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::VectorXd singularValues = svd.singularValues();
	double singularValue = (singularValues[0] + singularValues[1]) * 0.5;
	Eigen::DiagonalMatrix< double, 3, 3 > diagonal(singularValue, singularValue, 0.0);
	//Eigen::DiagonalMatrix< double, 3, 3 > diagonal(singularValues(0), singularValues(1), 0.0);
	Eigen::MatrixXd D = diagonal.toDenseMatrix();
	Eigen::MatrixXd F_restricted = svd.matrixU() * D * svd.matrixV().transpose();

	//std::cout << std::fixed << "U:" << std::endl << svd.matrixU() << std::endl << std::endl;
	//std::cout << std::fixed << "D:" << std::endl << D << std::endl << std::endl;
	//std::cout << std::fixed << "V:" << std::endl << svd.matrixV() << std::endl << std::endl;


	//std::cout << std::fixed << "F restricted:" << std::endl << F_restricted << std::endl << std::endl;

	//std::cout << std::fixed 
	//	<< "F restricted determinant: " << F_restricted.determinant() << std::endl
	//	<< "F determinant           : " << F.determinant() << std::endl << std::endl;

	return F_restricted;
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
	E = applyRestriction(E);
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

	//cout << "========== P Computation =============: " << endl
	//	<< "E: " << endl << E << endl << endl
	//	<< "U: " << endl << U << endl << endl
	//	<< "D: " << endl << svd.singularValues() << endl << endl
	//	<< "Vt: " << endl << Vt << endl << endl
	//	<< "W: " << endl << W << endl << endl
	//	<< "u3 :" << u3 << endl << endl;

	Eigen::MatrixXd P0(Eigen::MatrixXd::Identity(3, 4)); // Just K0: no translation or rotation.

	Eigen::MatrixXd P1noT = U * W * Vt; // P without translation part.

	Eigen::MatrixXd P1(3, 4);
	P1.block(0, 0, 3, 3) = P1noT;
	P1.block(0, 3, 3, 1) = u3;

	///cout << "1st sol: " << endl;

	Eigen::MatrixXd bestP1;
	int numCorrectSolutions = 0;
	if (checkP1(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
		//return P1;
	}

	P1.block(0, 3, 3, 1) = -u3;

	//cout << "2nd sol: " << endl;

	if (checkP1(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
		//return P1;
	}

	P1noT = U * W.transpose() * Vt;
	P1.block(0, 0, 3, 3) = P1noT;
	P1.block(0, 3, 3, 1) = u3;
	
	//cout << "3rd sol: " << endl;

	if (checkP1(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
		//return P1;
	}

	P1.block(0, 3, 3, 1) = -u3;

	///cout << "4th sol: " << endl;

	if (checkP1(pts, P0, P1))
	{
		++numCorrectSolutions;
		bestP1 = P1;
		//return P1;
	}

	//cout << "========== P Computation End ============= " << endl << endl
	//	<< "Number of correct solutions: " << numCorrectSolutions << endl << endl;

	if (numCorrectSolutions < 1) 
		std::cerr << "[Error]   : None of the results for P1 was accepted!" << std::endl << std::endl;
	if (numCorrectSolutions > 1) 
		std::cerr << "[Error]   : More than one solution found" << std::endl << std::endl;

	return bestP1;
}


bool Reconstruction3D::checkP1(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
								const Eigen::MatrixXd& P0, 
								const Eigen::MatrixXd& P1)
{
	Eigen::VectorXd p0 = pts[0].first;
	Eigen::VectorXd p1 = pts[0].second;

	std::vector< std::pair< Eigen::VectorXd, Eigen::VectorXd > > pair(1);
	pair[0] = std::pair< Eigen::VectorXd, Eigen::VectorXd >(p0, p1);

#if 0 // TODO
	TriangulationDlt dlt(pair, P0, P1);
	dlt.solve();
	VectorXd point3D = dlt.getPoint3D();

	VectorXd reprojectedImg0 = P0 * point3D;
	VectorXd reprojectedImg1 = P1 * point3D;

	cout << "========== P1 Check =============: " << endl
		<< "P0:" << endl << P0 << endl << endl
		<< "P1:" << endl << P1 << endl << endl
		<< "x: " << endl << p0 << endl << endl
		<< "x': " << endl << p1 << endl << endl
		<< "3d point: " << endl << point3D << endl << endl
		<< "Reprojected x: " << endl << reprojectedImg0 << endl << endl
		<< "Reprojected x': " << endl << reprojectedImg1 << endl << endl
		<< "========== P1 Check End =============: " << endl << endl;

	return reprojectedImg0[2] > 0. && reprojectedImg1[2] > 0.;
#endif
	return false;
}
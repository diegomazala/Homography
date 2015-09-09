#include <QApplication>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "Points.h"
#include "DLT.h"
#include "Reconstruction3D.h"
#include "Triangulation.h"
#include <fstream>




std::pair<Eigen::MatrixXd, Eigen::MatrixXd>					K (Eigen::MatrixXd(3, 3), Eigen::MatrixXd(3, 3));
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>					P (Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>					R (Eigen::MatrixXd(3, 3), Eigen::MatrixXd(3, 3));
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>					Rt(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
std::pair<Eigen::VectorXd, Eigen::VectorXd>					t (Eigen::VectorXd(3),    Eigen::VectorXd(3));

std::vector<Eigen::Vector3d>								Points3D;
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>	Points2D;
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>	Points2DNorm;


void setupMatrices()
{
	//
	// K Matrix
	//
	K.first(0, 0) = K.first(1, 1) = 1.73205;
	K.first(0, 2) = 0;
	K.first(1, 2) = 0;
	K.first(2, 2) = 1.0;
	K.second = K.first;
	

	//
	// R Matrix
	//
	R.first <<
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0;
	//
	// 30 degrees in Y axis
	R.second <<
		0.86603, 0.00000, 0.50000,
		0.00000, 1.00000, 0.00000,
		0.50000, 0.00000, 0.86603;


	//
	// t vector
	//
	t.first << 0.0, 0.0, -2.0;
	t.second << 1.0, 0.0, -2.0;


	//
	// Rt Matrix
	//
	Rt.first.block(0, 0, 3, 3) = R.first;
	Rt.second.block(0, 0, 3, 3) = R.second;
	Rt.first.col(3) = -R.first * t.first;
	Rt.second.col(3) = -R.second * t.second;

	//
	// P Matrix
	//
	P.first = Rt.first;
	P.second = Rt.second;


	std::cout << "K : " << std::endl << K.first << std::endl << std::endl;
	//std::cout << "K2 : " << std::endl << K.second << std::endl << std::endl;
	std::cout << "R1 : " << std::endl << R.first << std::endl << std::endl;
	std::cout << "R2 : " << std::endl << R.second << std::endl << std::endl;
	std::cout << "t1 : " << std::endl << t.first << std::endl << std::endl;
	std::cout << "t2 : " << std::endl << t.second << std::endl << std::endl;
	std::cout << "Rt1: " << std::endl << Rt.first << std::endl << std::endl;
	std::cout << "Rt2: " << std::endl << Rt.second << std::endl << std::endl;
	std::cout << "P1 : " << std::endl << P.first << std::endl << std::endl;
	std::cout << "P2 : " << std::endl << P.second << std::endl << std::endl;
}
	

void generate3DPoints()
{
	Points3D.clear();

	Points3D.push_back(Eigen::Vector3d(-0.5, -0.5, -0.5));
	Points3D.push_back(Eigen::Vector3d(-0.5, 0.5, -0.5));
	Points3D.push_back(Eigen::Vector3d(0.5, 0.5, -0.5));
	Points3D.push_back(Eigen::Vector3d(0.5, -0.5, -0.5));
	Points3D.push_back(Eigen::Vector3d(-0.5, -0.5, 0.5));
	Points3D.push_back(Eigen::Vector3d(-0.5, 0.5, 0.5));
	Points3D.push_back(Eigen::Vector3d(0.5, 0.5, 0.5));
	Points3D.push_back(Eigen::Vector3d(0.5, -0.5, 0.5));
}



void computePX()
{
	Points2D.clear();

	for (auto p3d : Points3D)
	{
		Eigen::VectorXd p0 = P.first * p3d.homogeneous();
		p0 = p0 / p0[2];
		Eigen::VectorXd p1 = P.second * p3d.homogeneous();
		p1 = p1 / p1[2];

		//std::cout
		//	<< "3D Point : " << p3d.transpose() << std::endl
		//	<< "P * X    : " << p0.transpose() << std::endl
		//	<< "P * X    : " << p1.transpose() << std::endl << std::endl;

		Points2D.push_back(std::make_pair(p0.head<2>(), p1.head<2>()));
	}
}


void exportObj(const std::string& filename, const std::vector<Eigen::Vector3d>& points3D)
{
	std::ofstream file;
	file.open(filename);
	for (const auto X : points3D)
	{
		file << "v " << X.transpose() << std::endl;
	}

	file
		<< "f 1 2 3 4" << std::endl
		<< "f 8 7 6 5" << std::endl
		<< "f 5 6 2 1" << std::endl
		<< "f 4 3 7 8" << std::endl
		<< "f 2 6 7 3" << std::endl
		<< "f 5 1 4 8" << std::endl;
	file.close();
}

void exportObj(const std::string& filename, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D)
{
	std::ofstream file;
	file.open(filename);
	for (const auto x : points2D)
	{
		const Eigen::VectorXd& X = Triangulation::solve(P, x);
		file << "v " << X.transpose() << std::endl;
	}

	file
		<< "f 1 2 3 4" << std::endl
		<< "f 8 7 6 5" << std::endl
		<< "f 5 6 2 1" << std::endl
		<< "f 4 3 7 8" << std::endl
		<< "f 2 6 7 3" << std::endl
		<< "f 5 1 4 8" << std::endl;
	file.close();
}



int main(int argc, char* argv[])
{
	setupMatrices();
	generate3DPoints();
	computePX();


	//
	// Normalize points
	// 
	
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d> T = DLT::normalizePoints(Points2D, Points2DNorm);


	///////////////////////////////////////////////////////////////////////////////////////
	//
	// Compute F matrix
	//
	Eigen::MatrixXd Fn = Reconstruction3D::computeF(Points2DNorm);
	
	//std::cout
	//	<< std::endl << std::fixed
	//	<< "F normalized: " << std::endl
	//	<< Fn << std::endl << std::endl;

	Fn = Reconstruction3D::applyConstraint(Fn);

	// 
	// Denormalize F matrix
	//
	Eigen::MatrixXd F = DLT::denormalizeH(Fn, T);

	
	std::cout
		<< std::endl << std::fixed
		<< "F constrained and denormalized: " << std::endl
		<< F << std::endl << std::endl;

	//F /= F(2, 2);
	//std::cout
	//	<< std::endl << std::fixed
	//	<< "F constrained, denormalized and divided by F(2,2): " << std::endl
	//	<< F << std::endl << std::endl;
	//
	/////////////////////////////////////////////////////////////////////////////////////// 


	Eigen::MatrixXd E = Reconstruction3D::computeE(K.first, F);
	std::cout
		<< std::endl << std::fixed
		<< "E: " << std::endl
		<< E << std::endl << std::endl;
	

	///////////////////////////////////////////////////////////////////////////////////////
	//
	// Compute P matrix
	//
	Eigen::MatrixXd Pmat = Reconstruction3D::computeP(Points2DNorm, E);
	//Pmat /= Pmat(2, 2);
	//std::cout
	//	<< std::endl << std::fixed
	//	<< "P0: " << std::endl
	//	<< P.first << std::endl << std::endl;

	//std::cout
	//	<< std::endl << std::fixed
	//	<< "P1: " << std::endl
	//	<< P.second << std::endl << std::endl;

	
	std::cout
		<< std::endl << std::fixed
		<< "Pmat: " << std::endl
		<< Pmat << std::endl << std::endl;


	std::cout << "Error x'Fx=0  : " << Reconstruction3D::computeError(Points2D, F) << std::endl;
	std::cout << "Error x'FnX=0 : " << Reconstruction3D::computeError(Points2D, Fn) << std::endl;
	std::cout << "Error x'EX=0  : " << Reconstruction3D::computeError(Points2D, E) << std::endl << std::endl;



	if (Pmat.isApprox(P.second))
		std::cout << "TEST RESULT: Pmat.isApprox(P) [OK]" << std::endl;
	else
		std::cout << "TEST RESULT: Pmat.isApprox(P) [FAILED]" << std::endl;


	std::cout << std::endl << std::endl;
	for (auto p : Points2DNorm)
	{
		Eigen::VectorXd X = Triangulation::solve(P, p);

		std::cout << std::fixed << X.transpose() << std::endl;
	}
	std::cout << std::endl << std::endl;
	for (auto p : Points2D)
	{
		Eigen::VectorXd X = Triangulation::solve(P, p);

		std::cout << std::fixed << X.transpose() << std::endl;
	}

	exportObj("../../data/Points2DNorm_Src.obj", Points2DNorm);
	exportObj("../../data/Points2D_Src.obj", Points2D);

	P.second = Pmat;

	std::cout << std::endl << std::endl;
	for (auto p : Points2DNorm)
	{
		Eigen::VectorXd X = Triangulation::solve(P, p);

		std::cout << std::fixed << X.transpose() << std::endl;
	}
	std::cout << std::endl << std::endl;
	for (auto p : Points2D)
	{
		Eigen::VectorXd X = Triangulation::solve(P, p);
		std::cout << std::fixed << X.transpose() << std::endl;
	}

	std::cout << std::endl << std::endl;

	exportObj("../../data/Points2DNorm_Pmat.obj", Points2DNorm);
	exportObj("../../data/Points2D_Pmat.obj", Points2D);

	return EXIT_SUCCESS;
}

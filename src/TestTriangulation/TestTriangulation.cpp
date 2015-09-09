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


	//std::cout << "K : " << std::endl << K.first << std::endl << std::endl;
	////std::cout << "K2 : " << std::endl << K.second << std::endl << std::endl;
	//std::cout << "R1 : " << std::endl << R.first << std::endl << std::endl;
	//std::cout << "R2 : " << std::endl << R.second << std::endl << std::endl;
	//std::cout << "t1 : " << std::endl << t.first << std::endl << std::endl;
	//std::cout << "t2 : " << std::endl << t.second << std::endl << std::endl;
	//std::cout << "Rt1: " << std::endl << Rt.first << std::endl << std::endl;
	//std::cout << "Rt2: " << std::endl << Rt.second << std::endl << std::endl;
	//std::cout << "P1 : " << std::endl << P.first << std::endl << std::endl;
	//std::cout << "P2 : " << std::endl << P.second << std::endl << std::endl;
}
	

void generate3DPoints()
{
	Points3D.clear();

#if 0
	Points3D.push_back(Eigen::Vector3d(-0.5, -0.5, -0.5));
	Points3D.push_back(Eigen::Vector3d(-0.5, 0.5, -0.5));
	Points3D.push_back(Eigen::Vector3d(0.5, 0.5, -0.5));
	Points3D.push_back(Eigen::Vector3d(0.5, -0.5, -0.5));
	Points3D.push_back(Eigen::Vector3d(-0.5, -0.5, 0.5));
	Points3D.push_back(Eigen::Vector3d(-0.5, 0.5, 0.5));
	Points3D.push_back(Eigen::Vector3d(0.5, 0.5, 0.5));
	Points3D.push_back(Eigen::Vector3d(0.5, -0.5, 0.5));
#else
	Points3D.push_back(Eigen::Vector3d(-1, -1, -1));
	Points3D.push_back(Eigen::Vector3d(-1, 1, -1));
	Points3D.push_back(Eigen::Vector3d(1, 1, -1));
	Points3D.push_back(Eigen::Vector3d(1, -1, -1));
	Points3D.push_back(Eigen::Vector3d(-1, -1, 1));
	Points3D.push_back(Eigen::Vector3d(-1, 1, 1));
	Points3D.push_back(Eigen::Vector3d(1, 1, 1));
	Points3D.push_back(Eigen::Vector3d(1, -1, 1));
#endif
}




void computePX()
{
	Points2D.clear();

	for (auto p3d : Points3D)
	{
		Eigen::Vector3d p0 = P.first * p3d.homogeneous();
		p0 = p0 / p0[2];
		Eigen::Vector3d p1 = P.second * p3d.homogeneous();
		p1 = p1 / p1[2];

		//std::cout
		//	<< "3D Point : " << p3d.transpose() << std::endl
		//	<< "P * X    : " << p0.transpose() << std::endl
		//	<< "P * X    : " << p1.transpose() << std::endl << std::endl;

		Points2D.push_back(std::make_pair(p0.head<2>(), p1.head<2>()));
	}
}




int main(int argc, char* argv[])
{
	setupMatrices();
	generate3DPoints();
	computePX();

	std::cout << "================ Test Triangulation ==============" << std::endl;
	std::cout << "==                                              ==" << std::endl << std::endl;

	//
	// P Matrix
	//
	P.first = Rt.first;
	P.second = Rt.second;

	std::cout << "P0 : " << std::endl << P.first << std::endl << std::endl;
	std::cout << "P1 : " << std::endl << P.second << std::endl << std::endl;

	int i = 0;
	for (auto x : Points2D)
	{
		Eigen::Vector3d X = Triangulation::solve(P, x).head<3>();
		std::string aprox_str = X.isApprox(Points3D[i]) ? "--[Ok]" : "-- [FAIL]";
		std::cout
			<< std::fixed << aprox_str << std::endl
			<< "2d Correspondence     : " << x.first.transpose() << "   " << x.second.transpose() << std::endl
			<< "Expected 3d point     : " << Points3D[i].transpose() << std::endl
			<< "Triangulated 3d point : " << X.transpose() << std::endl << std::endl;
		++i;
	}

	std::cout << "==                                              ==" << std::endl;
	std::cout << "================ Test Triangulation ==============" << std::endl;


	return EXIT_SUCCESS;
}

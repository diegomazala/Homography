#include <QApplication>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "Points.h"
#include "DLT.h"
#include "Reconstruction3D.h"
#include "Triangulation.h"
#include "ObjHelper.h"
#include <fstream>




std::pair<Eigen::MatrixXd, Eigen::MatrixXd>					P(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
std::pair<Eigen::Matrix3d, Eigen::Matrix3d>					K(Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero());
std::pair<Eigen::Matrix3d, Eigen::Matrix3d>					R(Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero());
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>					Rt(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
std::pair<Eigen::Vector3d, Eigen::Vector3d>					t;

std::vector<Eigen::Vector3d>								Points3D;
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>	Points2D;
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>	Points2DNorm;



void setupMatrices()
{
	//
	// K Matrix
	//
#if 1
	K.first(0, 0) = K.first(1, 1) = 1.73205;
	K.first(0, 2) = 0;
	K.first(1, 2) = 0;
	K.first(2, 2) = 1.0;
	K.second = K.first;
#else	
	#if 0
		K.first(0, 0) = K.first(1, 1) = 114.873 / 0.0130887;
		K.first(0, 2) = 1936;
		K.first(1, 2) = 1296;
		K.first(2, 2) = 1.0;
	#else
		K.second = K.first;
		K.first(0, 0) = K.first(1, 1) = 114.873;
		K.first(0, 2) = 1936 * 0.0130887;
		K.first(1, 2) = 1296 * 0.0130887;
		K.first(2, 2) = 1.0;
	#endif
#endif
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
	//t.first << 0.0, 0.0, -2.0;
	t.first << 0.0, 0.0, 0.0;
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
#if 0
	P.first = K.first * Rt.first;
	P.second = K.second * Rt.second;
#else
	P.first = Rt.first;
	P.second = Rt.second;
#endif


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
	//Eigen::MatrixXd Pmat = Reconstruction3D::computeP(Points2DNorm, E);
	std::vector<Eigen::MatrixXd> P_solutions;
	Reconstruction3D::computeP(Points2DNorm, E, P_solutions);
	//Pmat /= Pmat(2, 2);

	//std::cout
	//	<< std::endl << std::fixed
	//	<< "P0: " << std::endl
	//	<< P.first << std::endl << std::endl;

	//std::cout
	//	<< std::endl << std::fixed
	//	<< "P1: " << std::endl
	//	<< P.second << std::endl << std::endl;

	
	//std::cout
	//	<< std::endl << std::fixed
	//	<< "Pmat: " << std::endl
	//	<< Pmat << std::endl << std::endl;


	std::cout << "[Info]  Error x'Fx=0  : " << Reconstruction3D::computeError(Points2D, F) << std::endl;
	std::cout << "[Info]  Error x'FnX=0 : " << Reconstruction3D::computeError(Points2D, Fn) << std::endl;
	std::cout << "[Info]  Error x'EX=0  : " << Reconstruction3D::computeError(Points2D, E) << std::endl << std::endl;


	P.first = K.first * Rt.first;
	P.second = K.second * Rt.second;



	int i = 0;

	std::cout
		<< std::endl
		<< "[Info]  Exporting : Cube-Points2D_" << i << ".obj" << std::endl
		<< "[Info]  P Solution " << i << std::endl
		<< P.second << std::endl;
	std::string obj_file_name = "../../data/Cube-Points2D_" + std::to_string(i) + ".obj";
	exportCubeObj(obj_file_name, Points2D, P);


	for (auto m : P_solutions)
	{
		++i;

		P.second = K.second * m;

		std::cout
			<< std::endl
			<< "[Info]  Exporting : Cube-Points2D_" << i << ".obj" << std::endl
			<< "[Info]  P Solution " << i << std::endl
			<< P.second << std::endl;
		std::string obj_file_name = "../../data/Cube-Points2D_" + std::to_string(i) + ".obj";
		exportCubeObj(obj_file_name, Points2D, P);
	}


#if 0
	std::cout << "\n\n------------------------------ Computing epipole..." << std::endl;
	P.second = Reconstruction3D::computeP(F);

	P.second /= P.second(2, 2);

	std::cout
		<< std::endl << std::fixed
		<< "P': " << std::endl
		<< P.second << std::endl << std::endl;

	i++;
	std::cout
		<< std::endl
		<< "[Info]  Exporting : Cube-Points2D_" << i << ".obj" << std::endl
		<< "[Info]  P Solution " << i << std::endl
		<< P.second << std::endl;
	obj_file_name = "../../data/Cube-Points2D_" + std::to_string(i) + ".obj";
	exportCubeObj(obj_file_name, Points2D, P);
#endif
	return EXIT_SUCCESS;
}

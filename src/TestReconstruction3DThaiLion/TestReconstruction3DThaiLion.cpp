#include <QApplication>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "Points.h"
#include "DLT.h"
#include "Reconstruction3D.h"
#include "Triangulation.h"
#include "ObjHelper.h"
#include <fstream>




std::pair<Eigen::MatrixXd, Eigen::MatrixXd>					P (Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
std::pair<Eigen::Matrix3d, Eigen::Matrix3d>					K (Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero());
std::pair<Eigen::Matrix3d, Eigen::Matrix3d>					R (Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero());
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>					Rt(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
std::pair<Eigen::Vector3d, Eigen::Vector3d>					t;

std::vector<Eigen::Vector3d>								Points3D;
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>	Points2D;
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>	Points2DNorm;


void setupMatrices()
{

#if 0
	<MLRaster label = "DSC_0176.JPG">
		<VCGCamera	LensDistortion = "0 0"
		PixelSizeMm = "0.0130887 0.0130887"
		TranslationVector = "79.3959 -114.356 -499.541 1"
		CenterPx = "1936 1296"
		RotationMatrix = "0.980106 -0.0199563 0.197469 0 0.0558328 0.982476 -0.177828 0 -0.190459 0.185315 0.964045 0 0 0 0 1 "
		FocalMm = "114.873"
		ViewportPx = "3872 2592" / >
		<Plane semantic = "1" fileName = "DSC_0176.JPG" / >
		< / MLRaster>

		<MLRaster label = "DSC_0179.JPG">
		<VCGCamera	LensDistortion = "0 0"
		PixelSizeMm = "0.0130887 0.0130887"
		TranslationVector = "-227.173 -103.559 -460.851 1"
		CenterPx = "1936 1296"
		RotationMatrix = "0.914099 -0.0148061 -0.40522 0 -0.0540653 0.98596 -0.157987 0 0.401869 0.166324 0.900465 0 0 0 0 1 "
		FocalMm = "114.873"
		ViewportPx = "3872 2592" / >
		<Plane semantic = "1" fileName = "DSC_0179.JPG" / >
		< / MLRaster>
#endif

	//
	// K matrix
	//
#if 1 // pixel unit measure
	K.first(0, 0) = K.first(1, 1) = 114.873 / 0.0130887;
	K.first(0, 2) = 1936;
	K.first(1, 2) = 1296;
	K.first(2, 2) = 1.0;
#else	// mm unit measure
	K.first(0, 0) = K.first(1, 1) = 114.873;
	K.first(0, 2) = 1936 * 0.0130887;
	K.first(1, 2) = 1296 * 0.0130887;
	K.first(2, 2) = 1.0;
#endif
	K.second = K.first;


	//
	// R matrix
	//
	R.first <<
		0.980106, -0.0199563, 0.197469,
		0.0558328, 0.982476, -0.177828,
		-0.190459, 0.185315, 0.964045;

	R.second <<
		0.914099, -0.0148061, -0.40522,
		-0.0540653, 0.98596, -0.157987,
		0.401869, 0.166324, 0.900465;


	//
	// t vector
	//
	t.first << 79.3959, -114.356, -499.541;
	t.second << -227.173, -103.559, -460.851;

	//
	// Rt matrix
	//
	Rt.first.block(0, 0, 3, 3) = R.first;
	Rt.second.block(0, 0, 3, 3) = R.second;
	Rt.first.col(3) = -R.first * t.first;
	Rt.second.col(3) = -R.second * t.second;


	//
	// P matrix
	//
	//P.first = K.first * Rt.first;
	//P.second = K.second * Rt.second;
	P.first = Rt.first;
	P.second = Rt.second;

	std::cout << "K1 : " << std::endl << K.first << std::endl << std::endl;
	std::cout << "K2 : " << std::endl << K.second << std::endl << std::endl;
	std::cout << "R1 : " << std::endl << R.first << std::endl << std::endl;
	std::cout << "R2 : " << std::endl << R.second << std::endl << std::endl;
	std::cout << "t1 : " << std::endl << t.first << std::endl << std::endl;
	std::cout << "t2 : " << std::endl << t.second << std::endl << std::endl;
	std::cout << "Rt1: " << std::endl << Rt.first << std::endl << std::endl;
	std::cout << "Rt2: " << std::endl << Rt.second << std::endl << std::endl;
	std::cout << "P1 : " << std::endl << P.first << std::endl << std::endl;
	std::cout << "P2 : " << std::endl << P.second << std::endl << std::endl;
}
	
void thaiLionPointCorrespondence()
{
//		1600  960   1528  916
//		2164  631   2001  969
//		1924 1224   1620 1204
//		1912 1573   1635 1551
//		2235 1675   1929 1701
//		1623 1993   1378 1932
//		1845 2203   1570 2176
//		2184 2193   1867 2220
}





int main(int argc, char* argv[])
{
	setupMatrices();

	//importExportObjPoints("../../data/thai-lion/thai-lion.obj", "../../data/thai-lion-proj-1.obj", P.first);
	//importExportObjPoints("../../data/thai-lion/thai-lion.obj", "../../data/thai-lion-proj-2.obj", P.second);
	//return 0;

	std::vector<std::string> imageFiles;
	std::string pointsFile;

	for (int i = 0; i < argc; ++i)
	{
		std::string str = argv[i];

		if (str.find(".png") != std::string::npos || str.find(".PNG") != std::string::npos ||
			str.find(".jpg") != std::string::npos || str.find(".JPG") != std::string::npos ||
			str.find(".bmp") != std::string::npos || str.find(".BMP") != std::string::npos)
			imageFiles.push_back(str);

		if (str.find(".txt") != std::string::npos)
			pointsFile = str;
	}
	//
	// Reading points from file
	//
	if (!Points::readFromFile(pointsFile, Points2D))
	{
		std::cout << "[Error] Could not read points from file: " << pointsFile << std::endl;
		return EXIT_FAILURE;
	}
	else
	{
		std::cout << "[Info]  Count of points loaded from file = " << Points2D.size() << std::endl;
	}



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
	Pmat /= Pmat(2, 2);
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




	std::cout << std::endl << "../../data/ThaiLion-Points2DNorm_P.obj" << std::endl;
	for (auto p : Points2DNorm)
	{
		Eigen::VectorXd X = Triangulation::solve(P, p);
		std::cout << std::fixed << X.transpose() << std::endl;
	}
	exportCubeObj("../../data/ThaiLion-Points2DNorm_P.obj", Points2DNorm, P);
	


	std::cout << std::endl << "../../data/ThaiLion-Points2D_P.obj" << std::endl;
	for (auto p : Points2D)
	{
		Eigen::VectorXd X = Triangulation::solve(P, p);
		std::cout << std::fixed << X.transpose() << std::endl;
	}
	exportCubeObj("../../data/ThaiLion-Points2D_P.obj", Points2D, P);



	P.second = Pmat;



	std::cout << std::endl << "../../data/ThaiLion-Points2DNorm_Pmat.obj" << std::endl;
	for (auto p : Points2DNorm)
	{
		Eigen::VectorXd X = Triangulation::solve(P, p);
		std::cout << std::fixed << X.transpose() << std::endl;
	}
	exportCubeObj("../../data/ThaiLion-Points2DNorm_Pmat.obj", Points2DNorm, P);
	


	std::cout << std::endl << "../../data/ThaiLion-Points2D_Pmat.obj" << std::endl;
	for (auto p : Points2D)
	{
		Eigen::VectorXd X = Triangulation::solve(P, p);
		std::cout << std::fixed << X.transpose() << std::endl;
	}
	exportCubeObj("../../data/ThaiLion-Points2D_Pmat.obj", Points2D, P);
	

	return EXIT_SUCCESS;
}

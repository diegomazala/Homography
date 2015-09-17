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
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>	Points2DAll;
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>	Points2DSurf;


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

	//std::cout << "K1 : " << std::endl << K.first << std::endl << std::endl;
	//std::cout << "K2 : " << std::endl << K.second << std::endl << std::endl;
	//std::cout << "R1 : " << std::endl << R.first << std::endl << std::endl;
	//std::cout << "R2 : " << std::endl << R.second << std::endl << std::endl;
	//std::cout << "t1 : " << std::endl << t.first << std::endl << std::endl;
	//std::cout << "t2 : " << std::endl << t.second << std::endl << std::endl;
	//std::cout << "Rt1: " << std::endl << Rt.first << std::endl << std::endl;
	//std::cout << "Rt2: " << std::endl << Rt.second << std::endl << std::endl;
	//std::cout << "P1 : " << std::endl << P.first << std::endl << std::endl;
	//std::cout << "P2 : " << std::endl << P.second << std::endl << std::endl;
}

void thaiLionPointCorrespondence()
{
	Points2D.clear();

	Points2D.push_back(std::make_pair(Eigen::Vector2d(1600, 955), Eigen::Vector2d(1528, 916)));
	Points2D.push_back(std::make_pair(Eigen::Vector2d(2164, 960), Eigen::Vector2d(2001, 969)));
	Points2D.push_back(std::make_pair(Eigen::Vector2d(1924, 1224), Eigen::Vector2d(1620, 1204)));
	Points2D.push_back(std::make_pair(Eigen::Vector2d(1912, 1573), Eigen::Vector2d(1635, 1551)));
	Points2D.push_back(std::make_pair(Eigen::Vector2d(2235, 1675), Eigen::Vector2d(1929, 1701)));
	Points2D.push_back(std::make_pair(Eigen::Vector2d(1623, 1993), Eigen::Vector2d(1378, 1932)));
	Points2D.push_back(std::make_pair(Eigen::Vector2d(1845, 2203), Eigen::Vector2d(1570, 2176)));
	Points2D.push_back(std::make_pair(Eigen::Vector2d(2184, 2193), Eigen::Vector2d(1867, 2220)));
}

void buildCorrespondenceFrom3DPoints(const std::vector<Eigen::Vector3d>& points3D,
	const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& matP,
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D)
{
	points2D.clear();

	assert(matP.first.rows() == 3 && matP.first.cols() == 4
		&& matP.second.rows() == 3 && matP.second.cols() == 4);

	for (auto p3d : points3D)
	{
		Eigen::Vector3d x0 = matP.first * p3d.homogeneous();
		Eigen::Vector3d x1 = matP.second * p3d.homogeneous();

		x0 /= x0[2];
		x1 /= x1[2];

		points2D.push_back(std::make_pair(x0.head<2>(), x1.head<2>()));
	}

}


void exportPointCorrespondence(const std::string& filename, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D)
{
	std::ofstream outFile;
	outFile.open(filename);

	outFile << "POINT_COUNT: " << points2D.size() << std::endl;

	for (auto p : points2D)
		//outFile << std::fixed << '[' << p.first.transpose() << "]\t[" << p.second.transpose() << ']' << std::endl;
		outFile 
		<< std::fixed 
		<< '[' << p.first.x() << ", " << p.first.y()
		<< "]\t[" << p.second.x() << ", " << p.second.y() 
		<< ']' << std::endl;

	outFile.close();
}


void exportPSolutions(const std::vector<Eigen::MatrixXd>& P_solutions, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D)
{
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> P_mat;
	P_mat.first = Eigen::MatrixXd::Identity(3, 4);
	P_mat.first.block(0, 0, 3, 3) = K.first;

	int i = 0;
	for (auto m : P_solutions)
	{
		++i;

		P_mat.second = K.second * m;

		std::cout
			<< std::endl << std::endl
			<< "[Info]  Exporting : ThaiLion_" << i << ".obj ..." << std::endl << std::endl;
		std::string obj_file_name = "../../data/ThaiLion_" + std::to_string(i) + ".obj";
		exportObj(obj_file_name, points2D, P_mat);
	}
}

void test8Points()
{
	setupMatrices();
	thaiLionPointCorrespondence();

	P.first = K.first * Rt.first;
	P.second = K.second * Rt.second;

	readPointsFromObj("../../data/thai-lion.obj", Points3D, 36000);
	buildCorrespondenceFrom3DPoints(Points3D, P, Points2DAll);

	double threshold = 50;
	ReconstructionDLT dlt (Points2D);
	dlt.solve();
	dlt.inliersCount = Reconstruction3D::computeInliers(Points2D, dlt.F, threshold, dlt.error);
	std::cout << "Inliers: " << dlt.inliersCount << " --> " << dlt.error << std::endl;
	
	Eigen::MatrixXd E = Reconstruction3D::computeE(K, dlt.F);


	///////////////////////////////////////////////////////////////////////////////////////
	//
	// Compute P matrix
	//
	P.first = Eigen::MatrixXd::Identity(3, 4);
	P.first.block(0, 0, 3, 3) = K.first;
	//
	
	
	P.second = K.second * Reconstruction3D::computeP(Points2D, E);

	std::string obj_file_name = "../../data/ThaiLion_DLT_test8pts.obj";
	std::cout
		<< std::endl << std::endl
		<< "[Info]  Exporting : " << obj_file_name << std::endl << std::endl;

	exportObj(obj_file_name, Points2DAll, P);
}

void testRansac()
{
	setupMatrices();
	thaiLionPointCorrespondence();
	
	
	//Points2DSurf = Points2D;
	//
	// Reading points from file
	//
	std::string pointsFile = "../../data/ThaiLion-35991-pts.txt";
	if (!Points::readFromFile(pointsFile, Points2DSurf))
	{
		std::cout << "[Error] Could not read points from file: " << pointsFile << std::endl;
		return;
	}
	else
	{
		std::cout << "[Info]  Count of points loaded from file = " << Points2DSurf.size() << std::endl;
	}

	P.first = K.first * Rt.first;
	P.second = K.second * Rt.second;

	readPointsFromObj("../../data/thai-lion.obj", Points3D, 36000);
	buildCorrespondenceFrom3DPoints(Points3D, P, Points2DAll);

	//exportPointCorrespondence("../../data/ThaiLion_All-pts.txt", Points2DAll);
	//exportPointCorrespondence("../../data/ThaiLion_AllSurf-pts.txt", Points2DSurf);
	
	//Points2DSurf = Points2DAll;

	double threshold = 50;

	ReconstructionDLT dlt = Reconstruction3D::solve(Points2DSurf, K, threshold, 500);

	
	std::cout << "Inliers: " << dlt.inliersCount << " --> " << dlt.error << std::endl;

	Eigen::MatrixXd E = Reconstruction3D::computeE(K, dlt.F);

	///////////////////////////////////////////////////////////////////////////////////////
	//
	// Compute P matrix
	//
	P.first = Eigen::MatrixXd::Identity(3, 4);
	P.first.block(0, 0, 3, 3) = K.first;
	//
	std::vector<Eigen::MatrixXd> P_solutions;
	Reconstruction3D::computeP(Points2DSurf, E, P_solutions);

	double outlierThreshold = 20.0;
	int i = 0;
	for (auto m : P_solutions)
	{
		++i;
		P.second = K.second * m;

		std::string obj_file_name = "../../data/ThaiLion_Ransac_" + std::to_string(i) + ".obj";
		std::cout
			<< std::endl << std::endl
			<< "[Info]  Exporting : " << obj_file_name << std::endl << std::endl;

		//exportObj(obj_file_name, Points2DAll, P);
		exportObj(obj_file_name, Points2DSurf, P);
	}
}


int main(int argc, char* argv[])
{
	//test8Points();
	testRansac();
	return 0;

	setupMatrices();


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
	if (!Points::readFromFile("../../data/ThaiLion-35991-pts.txt", Points2DAll))
	{
		std::cout << "[Error] Could not read points from file: ../../data/ThaiLion-35991-pts.txt" << std::endl;
		return EXIT_FAILURE;
	}
	else
	{
		std::cout << "[Info]  Count of points loaded from file = " << Points2DAll.size() << std::endl;
	}



	//
	// Reading points from file
	//
	if (!Points::readFromFile(pointsFile, Points2DSurf))
	{
		std::cout << "[Error] Could not read points from file: " << pointsFile << std::endl;
		return EXIT_FAILURE;
	}
	else
	{
		std::cout << "[Info]  Count of points loaded from file = " << Points2DSurf.size() << std::endl;
	}
	

	//ReconstructionDLT dlt = Reconstruction3D::solve(Points2DSurf, K, 30, 5);
	//std::cout << "Inliers: " << dlt.inliersCount << " --> " << dlt.error << std::endl;

	double threshold = 50;
	ReconstructionDLT dlt(Points2DSurf);
	dlt.solve();
	dlt.inliersCount = Reconstruction3D::computeInliers(Points2DSurf, dlt.F, threshold, dlt.error);
	std::cout << "Inliers: " << dlt.inliersCount << " --> " << dlt.error << std::endl;

	Eigen::MatrixXd E = Reconstruction3D::computeE(K, dlt.F);
	//
	/////////////////////////////////////////////////////////////////////////////////////// 


	for (auto x : Points2DSurf)
	{
		std::cout << x.first.transpose() << "   " << x.second.transpose() << std::endl;
	}

	
	//std::cout << "Error Points2D     x'Ex=0 : " << std::fixed << Reconstruction3D::computeError(Points2D, E) << std::endl;
	//std::cout << "Error Points2DNorm x'Ex=0 : " << std::fixed << Reconstruction3D::computeError(Points2DNorm, E) << std::endl;


#if 1

	///////////////////////////////////////////////////////////////////////////////////////
	//
	// Compute P matrix
	//
	P.first = Eigen::MatrixXd::Identity(3, 4);
	P.first.block(0, 0, 3, 3) = K.first;
	//
	std::vector<Eigen::MatrixXd> P_solutions;
	Reconstruction3D::computeP(Points2DSurf, E, P_solutions);

	//exportPSolutions(P_solutions, Points2DAll);

	double outlierThreshold = 20.0;
	int i = 0;
	for (auto m : P_solutions)
	{
		++i;
		P.second = K.second * m;

		std::string obj_file_name = "../../data/ThaiLion_Ransac_" + std::to_string(i) + ".obj";
		std::cout
			<< std::endl << std::endl
			<< "[Info]  Exporting : " << obj_file_name << std::endl << std::endl;
		
		exportObj(obj_file_name, Points2DAll, P);
	}
#else
	///////////////////////////////////////////////////////////////////////////////////////
	//
	// Compute P matrix
	//
	P.first = Eigen::MatrixXd::Identity(3, 4);
	P.first.block(0, 0, 3, 3) = K.first;
	//
	Eigen::MatrixXd P2 = Reconstruction3D::computeP(Points2DNorm, E);
	P.second = K.second * P2;
	//

	std::string obj_file_name = "../../data/ThaiLion_solution.obj";
	std::cout
		<< std::endl << std::endl
		<< "[Info]  Exporting : " << obj_file_name << " ..." << std::endl << std::endl
		<< P.second << std::endl;
	//exportObj(obj_file_name, Points2D, P);
	exportObj(obj_file_name, Points2DAll, P);
#endif

	return EXIT_SUCCESS;
}

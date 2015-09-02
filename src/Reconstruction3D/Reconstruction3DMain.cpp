#include <QApplication>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "Points.h"
#include "DLT.h"
#include "Reconstruction3D.h"
#include <fstream>


void writeObjFile(const std::string& filename, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const Eigen::Matrix3d& mat)
{
	std::ofstream file;
	file.open(filename);
	
	for (const auto p : pts)
	{
		const auto& x0 = p.first.homogeneous();
		const auto& x1 = p.second.homogeneous();
		auto X1 = mat * x1;
		file << "v " << X1.transpose() << std::endl;
		auto X0 = mat * x0;
		file << "v " << X0.transpose() << std::endl;
	}

	file.close();
}

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		std::cout << "[Error] Missing parameters " << std::endl;
		std::cout << "[Usage] App.exe InputImageFileNameLeft.png InputImageFileNameRight.png PointsFile.txt OutputImageFileName.png" << std::endl;
		std::cout << "[Usage] App.exe ImageLeft.png ImageCenterRight.png ImageRight.png OutputImageFileName.png" << std::endl;
		std::cout << "[Usage] ./Reconstruction3D.exe ../../data/pier/1.jpg ../../data/pier/2.jpg ../../data/pier/3.jpg ../../data/out.png" << std::endl;
		std::cout << "[Usage] ./Reconstruction3D.exe ../../data/pier/1.jpg ../../data/pier/2.jpg ../../data/pier/3.jpg ../../data/out.png" << std::endl;
		std::cout << "[Usage] ./Reconstruction3D.exe ../../data/pier/1.jpg ../../data/pier/2.jpg ../../data/results/pier12.txt ../../data/results/pier12.png" << std::endl;
		return EXIT_FAILURE;
	}

	QApplication app(argc, argv);

	bool generatePoints = true;


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
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d> K(Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero());
	K.first(0, 0) = K.first(1, 1) = 114.873 * 0.0130887;
	K.first(0, 2) = 1936;
	K.first(1, 2) = 1296;
	K.second = K.first;

	//
	// R matrix
	//
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> R(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
	R.first << 0.980106, -0.0199563, 0.197469, 0, 0.0558328, 0.982476, -0.177828, 0, -0.190459, 0.185315, 0.964045, 0, 0, 0, 0, 1;
	R.second << 0.914099, -0.0148061, -0.40522, 0, -0.0540653, 0.98596, -0.157987, 0, 0.401869, 0.166324, 0.900465, 0, 0, 0, 0, 1;


	//
	// t vector
	//
	std::pair<Eigen::Vector3d, Eigen::Vector3d> t;
	t.first << 79.3959, -114.356, -499.541;
	t.second << -227.173, -103.559, -460.851;

	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Rt(R.first, R.second);
	Rt.first.col(3) = t.first;
	Rt.second.col(3) = t.second;



	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> P;
	P.first = K.first * Rt.first;
	P.second = K.second * Rt.second;

	std::cout << "K : " << std::endl << K.first << std::endl << std::endl;
	std::cout << "R : " << std::endl << R.first << std::endl << std::endl;
	std::cout << "t : " << std::endl << t.first << std::endl << std::endl;
	std::cout << "Rt: " << std::endl << Rt.first << std::endl << std::endl;
	std::cout << "P : " << std::endl << P.first << std::endl << std::endl;

	std::pair<QImage, QImage> image;
	QImage outputImage;


	//
	// Reading points from file
	//
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> pointsSrc;
	if (!Points::readFromFile(pointsFile, pointsSrc))
	{
		std::cout << "[Error] Could not read points from file: " << pointsFile << std::endl;
		return EXIT_FAILURE;
	}
	else
	{
		std::cout << "[Info]  Count of points loaded from file = " << pointsSrc.size() << std::endl;
	}



	//
	// Normalize points
	// 
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> pointsNorm;
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d> T = DLT::normalizePoints(pointsSrc, pointsNorm);
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> points(pointsNorm.begin(), pointsNorm.begin() + 8);



	///////////////////////////////////////////////////////////////////////////////////////
	//
	// Compute F matrix
	//
	Eigen::MatrixXd Fn = Reconstruction3D::computeF(points);
	
	std::cout
		<< std::endl << std::fixed
		<< "F normalized: " << std::endl
		<< Fn << std::endl << std::endl;
	Fn = Reconstruction3D::applyConstraint(Fn);

	// 
	// Denormalize F matrix
	//
	Eigen::MatrixXd F = DLT::denormalizeH(Fn, T);
	
	std::cout
		<< std::endl << std::fixed
		<< "F constrained and denormalized: " << std::endl
		<< F << std::endl << std::endl;

	F /= F(2, 2);
	//std::cout
	//	<< std::endl << std::fixed
	//	<< "F restricted, denormalized and normalized by F(2,2): " << std::endl
	//	<< F << std::endl << std::endl;
	//
	/////////////////////////////////////////////////////////////////////////////////////// 

	

	Eigen::MatrixXd E = Reconstruction3D::computeE(K.first, F);



	///////////////////////////////////////////////////////////////////////////////////////
	//
	// Compute P matrix
	//
	Eigen::MatrixXd Pmat = Reconstruction3D::computeP(points, E);

	std::cout
		<< std::endl << std::fixed
		<< "Pmat: " << std::endl
		<< Pmat << std::endl << std::endl;


	writeObjFile("../../data/pts.obj", pointsNorm, Pmat.inverse());
	
	//QImageWidget outputWidget;
	//outputWidget.setImage(outputImage);
	//outputWidget.show();
	//return app.exec();
	return EXIT_SUCCESS;
}

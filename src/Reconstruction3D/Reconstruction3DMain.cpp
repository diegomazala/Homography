#include <QApplication>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "Points.h"
#include "DLT.h"
#include "Reconstruction3D.h"




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

	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> R(Eigen::MatrixXd(3, 4), Eigen::MatrixXd(3, 4));
	R.first << 0.980106, -0.0199563, 0.197469, 0, 0.0558328, 0.982476, -0.177828, 0, -0.190459, 0.185315, 0.964045, 0, 0, 0, 0, 1;
	R.second << 0.914099, -0.0148061, -0.40522, 0, -0.0540653, 0.98596, -0.157987, 0, 0.401869, 0.166324, 0.900465, 0, 0, 0, 0, 1;

	std::pair<Eigen::Vector3d, Eigen::Vector3d> t;
	t.first << 79.3959, -114.356, -499.541;
	t.second << -227.173, -103.559, -460.851;


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
	

	//
	// Compute F matrix
	//
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> points(pointsNorm.begin(), pointsNorm.begin() + 7);
	Eigen::MatrixXd Fn = Reconstruction3D::computeF(points);


	// 
	// Denormalize F matrix
	//
	Eigen::MatrixXd F = DLT::denormalizeH(Fn, T);

	std::cout
		<< std::endl << std::fixed
		<< "F: " << std::endl
		<< F << std::endl << std::endl;

	//QImageWidget outputWidget;
	//outputWidget.setImage(outputImage);
	//outputWidget.show();
	//return app.exec();
	return EXIT_SUCCESS;
}

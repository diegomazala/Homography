#include <QApplication>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "DLT.h"
#include "Points.h"
#include "RansacDLT.h"
#include "GaussNewton.h"



void showImageWidgets(const std::pair<QImage, QImage>& input_image, const QImage& output_image, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points)
{
	std::pair<QImageWidget, QImageWidget> imageWidget;
	std::vector<Eigen::Vector2d> points_left, points_right;
	points_left.reserve(points.size());
	std::for_each(points.begin(), points.end(),
		[&points_left](const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>::value_type& p)
	{ points_left.push_back(p.first); });

	points_right.reserve(points.size());
	std::for_each(points.begin(), points.end(),
		[&points_right](const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>::value_type& p)
	{ points_right.push_back(p.second); });

	imageWidget.first.setImage(input_image.first);
	imageWidget.second.setImage(input_image.second);

	imageWidget.first.setPoints(points_left);
	imageWidget.second.setPoints(points_right);

	imageWidget.first.show();
	imageWidget.second.show();

	QImageWidget outputWidget;
	outputWidget.setImage(output_image);
	outputWidget.show();
}



static int runProgramReadingPointsFile(	const std::pair<std::string, std::string>& inputImageFileName, 
										const std::string& pointsFile, 
										QImage& outputImage, 
										const std::string outputFileName)
{
	std::pair<QImage, QImage> image;

	//
	// Reading points from file
	//
	Points points;
	if (!points.readFromFile(pointsFile))
	{
		std::cout << "[Error] Could not read points from file: " << pointsFile << std::endl;
		return EXIT_FAILURE;
	}
	else
	{
		std::cout << "[Info]  Count of points loaded from file = " << points.count() << std::endl;
	}


	DLT& dltRansac = RansacDLT::solve(points);

	std::cerr
		<< std::fixed << std::endl
		<< "[Info]  DLT RANSAC           : " << std::endl
		<< "[Info]  Projection Error     : " << dltRansac.getError().first << ", " << dltRansac.getError().second << std::endl
		<< "[Info]  Inliers Count        : " << dltRansac.getInliersCount() << " of " << points.count() << std::endl
		<< "[Info]  Inliers Percentage   : " << double(dltRansac.getInliersCount()) / double(points.count()) * 100.0 << "%"
		<< std::endl << std::endl;


	// Recomputing DLT from inliers.
	// This is useful for distribution of the error between the inliers
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> inliers_pts;
	DLT::getInliers(points.getPointArray(), dltRansac.getH(), inliers_pts);
	DLT dltInliers(inliers_pts);
	dltInliers.computeHomography();
	dltInliers.computeGeometricError();
	dltInliers.computeInliers(inliers_pts);

	std::cerr
		<< std::fixed << std::endl
		<< "[Info]  DLT from all inliers : " << std::endl
		<< "[Info]  Projection Error     : " << dltInliers.getError().first << ", " << dltInliers.getError().second << std::endl
		<< "[Info]  Inliers Count        : " << dltInliers.getInliersCount() << " of " << points.count() << std::endl
		<< "[Info]  Inliers Percentage   : " << double(dltInliers.getInliersCount()) / double(points.count()) * 100.0 << "%"
		<< std::endl << std::endl;


	

	DLT& dlt = dltRansac; 
				//dltInliers;

	double original_error = GaussNewton::getSumError(dlt.getPoints(), dlt.getH());
	Eigen::MatrixXd Hgn;
	double gauss_newton_error = GaussNewton::solve(dlt.getPoints(), dlt.getH(), 10, Hgn);


	
	std::cout
		<< std::fixed << std::endl
		<< "[Info]  H Inliers: " << std::endl << dltRansac.getH() << std::endl << std::endl
		<< "[Info]  H Ransac: " << std::endl << dltInliers.getH() << std::endl << std::endl
		<< "[Info]  H Gauss Newton: " << std::endl << Hgn << std::endl << std::endl;


	image.first.load(inputImageFileName.first.c_str());
	image.second.load(inputImageFileName.second.c_str());

	
	projectImages(dltRansac.getH(), image, outputImage);
	outputImage.save("out_dlt_ransac.png");
	projectImages(dltInliers.getH(), image, outputImage);
	outputImage.save("out_dlt_inliers.png");
	projectImages(Hgn, image, outputImage);
	outputImage.save("out_gauss_newton.png");
	

	return EXIT_SUCCESS;
}



static int runProgramGeneratingMatchingPoints(const std::vector<std::string>& imageFiles, QImage& outputImage)
{
	std::pair<QImage, QImage> image;

	// loop throught the image pairs ignoring the last one, which is the output image
	std::pair<std::string, std::string> inputImageFileName;
	inputImageFileName.first = imageFiles[0];
	for (int i = 1; i < imageFiles.size() - 1; ++i)
	{
		inputImageFileName.second = imageFiles[i];
		//
		// generating points for a pair of images
		//
		std::ostringstream cmd;
		cmd << "surf_matcher.exe " << inputImageFileName.first << " " << inputImageFileName.second << " out_surf.png";

		std::cout << std::endl << cmd.str() << std::endl << std::endl;
		std::cout << "[Info]  Wait. Looking for matching points ... " << std::endl << std::endl;

		system(cmd.str().c_str());
		
		std::cout 
			<< std::endl
			<< "[Info]  ... Matching points found." << std::endl;
		//std::cout << "[Info]  Press ESC to close window and continue program" << std::endl;


		const std::string outputImageFileName(imageFiles.back());
		std::string pointsFile = "surf_pts.txt";


		//
		// Reading points from file
		//
		Points points;
		if (!points.readFromFile(pointsFile))
		{
			std::cout << "[Error] Could not read points from file: " << pointsFile << std::endl;
			return EXIT_FAILURE;
		}
		else
		{
			std::cout << "[Info]  Count of points loaded from file = " << points.count() << std::endl;
		}


		DLT& dltRansac = RansacDLT::solve(points);

		std::cerr
			<< std::fixed << std::endl
			<< "[Info]  DLT RANSAC           : " << std::endl
			<< "[Info]  Projection Error     : " << dltRansac.getError().first << ", " << dltRansac.getError().second << std::endl
			<< "[Info]  Inliers Count        : " << dltRansac.getInliersCount() << " of " << points.count() << std::endl
			<< "[Info]  Inliers Percentage   : " << double(dltRansac.getInliersCount()) / double(points.count()) * 100.0 << "%"
			<< std::endl << std::endl;
		
		
		// Recomputing DLT from inliers.
		// This is useful for distribution of the error between the inliers
		std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> inliers_pts;
		DLT::getInliers(points.getPointArray(), dltRansac.getH(), inliers_pts);
		DLT dltInliers(inliers_pts);
		dltInliers.computeHomography();
		dltInliers.computeGeometricError();
		dltInliers.computeInliers(inliers_pts);

		std::cerr
			<< std::fixed << std::endl
			<< "[Info]  DLT from all inliers : " << std::endl
			<< "[Info]  Projection Error     : " << dltInliers.getError().first << ", " << dltInliers.getError().second << std::endl
			<< "[Info]  Inliers Count        : " << dltInliers.getInliersCount() << " of " << points.count() << std::endl
			<< "[Info]  Inliers Percentage   : " << double(dltInliers.getInliersCount()) / double(points.count()) * 100.0 << "%"
			<< std::endl << std::endl;

#if 1

		//
		// Non Linear step : Gauss Newton
		// 
		DLT& dlt = dltRansac;	//dltInliers;

		double original_error = GaussNewton::getSumError(dlt.getPoints(), dlt.getH());
		Eigen::MatrixXd Hgn;
		double gauss_newton_error = GaussNewton::solve(dlt.getPoints(), dlt.getH(), 10, Hgn);



		std::cout
			<< std::fixed << std::endl
			<< "[Info]  H Inliers: " << std::endl << dltRansac.getH() << std::endl << std::endl
			<< "[Info]  H Ransac: " << std::endl << dltInliers.getH() << std::endl << std::endl
			<< "[Info]  H Gauss Newton: " << std::endl << Hgn << std::endl << std::endl;


		image.first.load(inputImageFileName.first.c_str());
		image.second.load(inputImageFileName.second.c_str());

		
		projectImages(dltRansac.getH(), image, outputImage);
		outputImage.save("out_dlt_ransac_" + QString::number(i) + ".png");
		projectImages(dltInliers.getH(), image, outputImage);
		outputImage.save("out_dlt_inliers_" + QString::number(i) + ".png");
		projectImages(Hgn, image, outputImage);
		outputImage.save("out_gauss_newton_" + QString::number(i) + ".png");
#else


		image.first.load(inputImageFileName.first.c_str());
		image.second.load(inputImageFileName.second.c_str());

		//projectImages(dltRansac.getH(), image, outputImage);
		projectImages(dltInliers.getH(), image, outputImage);

		outputImage.save(outputImageFileName.c_str());

		
#endif
		inputImageFileName.first = imageFiles.back();
	}

	
	
	return EXIT_SUCCESS;
}


int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		std::cout << "[Error] Missing parameters " << std::endl;
		std::cout << "[Usage] App.exe InputImageFileNameLeft.png InputImageFileNameRight.png PointsFile.txt OutputImageFileName.png" << std::endl;
		std::cout << "[Usage] App.exe ImageLeft.png ImageCenterRight.png ImageRight.png OutputImageFileName.png" << std::endl;
		std::cout << "[Usage] ./Mosaic.exe ../../data/pier/1.jpg ../../data/pier/2.jpg ../../data/pier/3.jpg ../../data/out.png" << std::endl;
		std::cout << "[Usage] ./Mosaic.exe ../../data/pier/1.jpg ../../data/pier/2.jpg ../../data/results/pier12.txt ../../data/results/pier12.png" << std::endl;
		std::cout << "[Usage] ./Mosaic.exe ../../data/results/pier12.png ../../data/pier/3.jpg ../../data/results/pier123.txt ../../data/results/pier123.png" << std::endl;
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


	QImage outputImage;
	if (pointsFile.empty())
	{
		runProgramGeneratingMatchingPoints(imageFiles, outputImage);
	}
	else
	{
		runProgramReadingPointsFile(std::make_pair(imageFiles[0], imageFiles[1]), pointsFile, outputImage, imageFiles.back());
	}


	QImageWidget outputWidget;
	outputWidget.setImage(outputImage);
	outputWidget.show();


	return app.exec();
}

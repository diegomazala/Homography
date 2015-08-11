#include <QApplication>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "DLT.h"
#include "Points.h"
#include "RansacDLT.h"



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



static int runProgramReadingPointsFile(const std::pair<std::string, std::string>& inputImageFileName, const std::string& pointsFile, QImage& outputImage)
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


	DLT& dltConsensus = RansacDLT::solve(points);

	image.first.load(inputImageFileName.first.c_str());
	image.second.load(inputImageFileName.second.c_str());
	std::cerr
		<< std::fixed
		<< "[Info]  Projection Error  : " << dltConsensus.getError().first << ", " << dltConsensus.getError().second << std::endl
		<< "[Info]  Inliers Count     : " << dltConsensus.getInliersCount() << " of " << points.count() << std::endl
		<< "[Info]  Inliers Percentage: " << double(dltConsensus.getInliersCount()) / double(points.count()) * 100.0 << "%"
		<< std::endl << std::endl;

	projectImages(dltConsensus.getH(), image, outputImage);

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
		std::cout << "[Info]  Wait. Looking for matching points ... " << std::endl;

		system(cmd.str().c_str());
		
		std::cout << "[Info]  Press ESC to close window and continue program" << std::endl;


		std::string outputImageFileName(imageFiles.back());
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


		DLT& dltConsensus = RansacDLT::solve(points);

		image.first.load(inputImageFileName.first.c_str());
		image.second.load(inputImageFileName.second.c_str());
		std::cerr
			<< std::fixed
			<< "[Info]  Projection Error  : " << dltConsensus.getError().first << ", " << dltConsensus.getError().second << std::endl
			<< "[Info]  Inliers Count     : " << dltConsensus.getInliersCount() << " of " << points.count() << std::endl
			<< "[Info]  Inliers Percentage: " << double(dltConsensus.getInliersCount()) / double(points.count()) * 100.0 << "%"
			<< std::endl << std::endl;

		projectImages(dltConsensus.getH(), image, outputImage);

		outputImage.save(outputImageFileName.c_str());

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
		std::cout << "[Usage] ./Homography.exe ../../data/pier/1.jpg ../../data/pier/2.jpg ../../data/pier/3.jpg ../../data/out.png" << std::endl;
		std::cout << "[Usage] ./Homography.exe ../../data/pier/1.jpg ../../data/pier/2.jpg ../../data/results/pier12.txt ../../data/results/pier12.png" << std::endl;
		std::cout << "[Usage] ./Homography.exe ../../data/results/pier12.png ../../data/pier/3.jpg ../../data/results/pier123.txt ../../data/results/pier123.png" << std::endl;
		return EXIT_FAILURE;
	}

	QApplication app(argc, argv);

	bool generatePoints = true;


	std::vector<std::string> imageFiles;
	std::string pointsFile;

	for (int i = 0; i < argc; ++i)
	{
		std::string str = argv[i];

		if (str.find(".png") != std::string::npos || str.find(".jpg") != std::string::npos || str.find(".bmp") != std::string::npos)
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
		runProgramReadingPointsFile(std::make_pair(imageFiles[0], imageFiles[1]), pointsFile, outputImage);
	}


	outputImage.save(imageFiles.back().c_str());
	

	QImageWidget outputWidget;
	outputWidget.setImage(outputImage);
	outputWidget.show();


	return app.exec();
}

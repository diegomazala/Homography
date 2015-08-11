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




int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cout << "[Error] Missing parameters " << std::endl;
		std::cout << "[Usage] App.exe InputImageFileNameLeft.png InputImageFileNameRight.png PointsFile.txt OutputImageFileName.png" << std::endl;
		return EXIT_FAILURE;
	}

	std::pair<std::string, std::string> inputImageFileName(argv[1], argv[2]);
	std::string outputImageFileName = "output.png";
	if (argc > 3)
		outputImageFileName = argv[4];


	//
	// Reading points from file
	//
	Points points;
	if (!points.readFromFile(argv[3]))
	{
		std::cout << "[Error] Could not read points from file: " << argv[1] << std::endl;
		return EXIT_FAILURE;
	}
	else
	{
		std::cout << "[Info]  Count of points loaded from file = " << points.count() << std::endl;
	}


	DLT& dltConsensus = RansacDLT::solve(points);


	//std::cout << std::endl;
	//for (auto it : dltArray)
	//	std::cout << std::fixed << it.getInliersCount() << " : " << it.getError().first + it.getError().second << std::endl;

	Eigen::Matrix3d H;
	//H << 1.027308, -0.004961, -297.475919,
	//	0.066875,     1.014096, -54.126748,
	//	0.000312,     0.000044,    0.878409;

	H << 0.921571, 0.011174, -556.336662,
		0.061009, 0.951503, -76.333429,
		0.000328, 0.000050, 0.738080;


	std::cout << "\nH: " << std::endl << dltConsensus.getH() << std::endl << std::endl;

	//return 0;

#if 1
	QImage outputImage;
	std::pair<QImage, QImage> image;
	image.first.load(inputImageFileName.first.c_str());
	image.second.load(inputImageFileName.second.c_str());
	std::cout 
		<< std::fixed 
		<< "[Info]  Projection Error  : " << dltConsensus.getError().first << ", " << dltConsensus.getError().second << std::endl
		<< "[Info]  Inliers Count     : " << dltConsensus.getInliersCount() << " of " << points.count() << std::endl
		<< "[Info]  Inliers Percentage: " << double(dltConsensus.getInliersCount()) / double(points.count()) * 100.0 << "%"
		<< std::endl << std::endl;

	projectImages(dltConsensus.getH(), image, outputImage);


	//projectImages(H, image, outputImage);
	outputImage.save(outputImageFileName.c_str());

	QApplication app(argc, argv);
	
	QImageWidget outputWidget;
	outputWidget.setImage(outputImage);
	outputWidget.show();

	return app.exec();
#else
	return 0;
#endif
}

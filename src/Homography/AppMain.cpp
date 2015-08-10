#include <QApplication>
#include <random>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "DLT.h"
#include "Points.h"


static std::vector<int> random4Indices(int min, int max)
{
	std::random_device rd;								// obtain a random number from hardware
	std::default_random_engine eng(rd());				// seed the generator
	std::uniform_int_distribution<int> distr(min, max); // define the range

	std::vector<int> indices;
	indices.push_back(distr(eng));
	indices.push_back(distr(eng));
	indices.push_back(distr(eng));
	indices.push_back(distr(eng));


	// Check if there are repeated indices
	auto it = indices.begin();
	while (it != indices.end() - 1)
	{
		if (std::find(it + 1, indices.end(), *it) != indices.end())
			*it = distr(eng);
		else
			++it;
	}
	return indices;
}

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




	int iterations = points.count() * 0.5;	// define the number of iterations as 50% of the number of points
	std::vector<DLT> dltArray;
	dltArray.reserve(iterations);

	int inliersPercentageAchieved = 0;
	int inliersPercentageAccepted = 99;
	//
	// Run the iteration to find the best DLT, i.e, with the biggest number of inliers
	//
	//for (int i = 0; i < iterations; ++i)
	int i = 0;
	while (inliersPercentageAchieved < inliersPercentageAccepted && i++ < points.count())
	{
		dltArray.push_back(DLT());
		DLT& dlt = dltArray.back();

		//
		// Generating 4 indices for points in order to compute homography using these points
		// 
		std::cout << std::endl
			<< "[Info]  Iteration = " << i << std::endl
			<< "[Info]  Random 4 indices = ";
		std::vector<int> indices = random4Indices(0, (int)points.count() - 1);

		// copy points from original array to the 4-array to be used for homography
		for each (auto i in indices)
		{
			dlt.addPoint(points[i]);
			std::cout << i << "  ";
		}

		Eigen::MatrixXd H = dlt.computeHomography();
		dlt.computeError();


		//
		// project all points and count number of inliers and outliers
		//
		int inliers = dlt.computeInliers(points.getPointArray());

		inliersPercentageAchieved = double(inliers) / double(points.count()) * 100;
	}

	std::cout 
		<< std::endl
		<< "[Info]  Number of Iterations : " << i << std::endl
		<< "[Info]  Inliers Acheived     : " << inliersPercentageAchieved << "%" << std::endl << std::endl;

	// Sort the dlts by error in crescent order
	//std::sort(dltArray.begin(), dltArray.end());
	DLT& dltConsensus = dltArray.back();
	
	//std::cout << std::endl;
	//for (auto it : dltArray)
	//	std::cout << std::fixed << it.getInliersCount() << " : " << it.getError().first + it.getError().second << std::endl;

	//Eigen::Matrix3d H;
	//H << 1.027308, -0.004961, -297.475919,
	//	0.066875,     1.014096, -54.126748,
	//	0.000312,     0.000044,    0.878409;


	//std::cout << "\nH: " << std::endl << dltConsensus.getH() << std::endl << std::endl;



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

	QApplication a(argc, argv);
	
	QImageWidget outputWidget;
	outputWidget.setImage(outputImage);
	outputWidget.show();

	return a.exec();
#else
	return 0;
#endif
}

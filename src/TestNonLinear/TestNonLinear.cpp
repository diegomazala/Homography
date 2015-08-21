#include <QApplication>
#include "QImageWidget.h"
#include "ImageHelper.h"
#include "DLT.h"
#include "Points.h"
#include "RansacDLT.h"
#include "GaussNewton.h"




int main(int argc, char* argv[])
{
	QApplication app(argc, argv);

	const std::string outputFileName = "out_test.png";
	QImage outputImage;
	std::pair<std::string, std::string> imageFilename("../../data/pier/1.jpg", "../../data/pier/2.jpg");
	std::pair<QImage, QImage> image;
	if (!image.first.load(imageFilename.first.c_str()))
	{
		std::cerr << "[Error] Could not load image filr: " << imageFilename.first << std::endl;
		return EXIT_FAILURE;
	}

	if (!image.second.load(imageFilename.second.c_str()))
	{
		std::cerr << "[Error] Could not load image filr: " << imageFilename.second << std::endl;
		return EXIT_FAILURE;
	}


	DLT dlt;
	dlt.addPoint(std::make_pair(Eigen::Vector2d(370.361, 198.972), Eigen::Vector2d(81.6814, 171.918)));
	dlt.addPoint(std::make_pair(Eigen::Vector2d(378.627, 206.653), Eigen::Vector2d(89.9805, 179.755)));
	dlt.addPoint(std::make_pair(Eigen::Vector2d(413.459, 214.739), Eigen::Vector2d(124.377, 188.050)));
	dlt.addPoint(std::make_pair(Eigen::Vector2d(407.815, 209.960), Eigen::Vector2d(118.867, 183.199)));


	Eigen::MatrixXd H = dlt.computeHomography();
	//std::cout << "H = \n" << H << std::endl << std::endl;
	
	projectImages(H, image, outputImage);
	outputImage.save("test_out_H.png");

	double sum_error = GaussNewton::getSumError(dlt.getPoints(), H);
	std::cout << "\nDLT H sum_error: " << sum_error << std::endl << std::endl << std::endl;

	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> hp;
	Points::projectPoints(dlt.getPoints(), hp, H);

	//for (int i = 0; i < hp.size(); ++i)
	//{
	//	std::cout << dlt.getPoints()[i].first.transpose() << "  " << hp[i].first.transpose() << "   "
	//		<< dlt.getPoints()[i].second.transpose() << "  " << hp[i].second.transpose() << std::endl;
	//}

	H(0, 2) += 0.1;
	H(1, 2) += -0.1;

	//std::cout << "\ndisturbed H = \n" << H << std::endl << std::endl;
	sum_error = GaussNewton::getSumError(dlt.getPoints(), H);
	std::cout << "\nH disturbed sum_error: " << sum_error << std::endl << std::endl << std::endl;

	projectImages(H, image, outputImage);
	outputImage.save("out_H_disturbed.png");
	//return 0;

	Points::projectPoints(dlt.getPoints(), hp, H);
	//for (int i = 0; i < hp.size(); ++i)
	//{
	//	std::cout << dlt.getPoints()[i].first.transpose() << "  " << hp[i].first.transpose() << "   "
	//		<< dlt.getPoints()[i].second.transpose() << "  " << hp[i].second.transpose() << std::endl;
	//}



	//projectImages(H, image, outputImage);
	//outputImage.save("out_H_Jacobian.png");
	//sum_error = getSumError(dlt.getPoints(), H);
	//std::cout << "\Jacobian H sum_error: " << sum_error << std::endl;
	//return 0;



	const int maxIt = 1;
	bool error_maximized = false;

	//for (int it = 0; it < 3; ++it)
	int it = 0;
	while (it++ < maxIt && !error_maximized)
	{
		std::cout << "\n\n-------- begin loop ------  " << it - 1 << std::endl << std::endl;

		Eigen::MatrixXd dH = GaussNewton::computeDeltaH(dlt.getPoints(), H);

		H = H + dH;

		std::cout << "H = \n" << H << std::endl << std::endl;

		Points::projectPoints(dlt.getPoints(), hp, H);
		double new_sum_error = GaussNewton::getSumError(dlt.getPoints(), hp);

		if (sum_error > new_sum_error)
		{
			std::cout << std::fixed << "\n ** error minimized    : " << sum_error << " -> " << new_sum_error << std::endl;
			sum_error = new_sum_error;
		}
		else
		{
			std::cout << std::fixed << "\n ** error maximized    : " << sum_error << " -> " << new_sum_error << std::endl;
			error_maximized = true;
		}

		std::cout << "-------------- end loop ----------------------" << std::endl << std::endl;

		//if (!error_maximized)
		//{
		//	projectImages(H, image, outputImage);
		//	QString filename("out_H_" + QString::number(it) + ".png");
		//	outputImage.save(filename);
		//}

	}




	//projectImages(H, image, outputImage);
	////outputImage.save(outputFileName.c_str());
	//outputImage.save("out_H.png");


	//projectImages(HJ, image, outputImage);
	////outputImage.save(outputFileName.c_str());
	//outputImage.save("out_HJ.png");



	

	//QImageWidget outputWidget;
	//outputWidget.setImage(outputImage);
	//outputWidget.show();

	std::cout << "\n************** // ********************\n";

	return 0;
	//return app.exec();
}

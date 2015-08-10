#ifndef __IMAGE_HELPER_H__
#define __IMAGE_HELPER_H__

#include <QApplication>
#include <QImage>
#include <iostream>
#include <Eigen/Dense>


static void computImageSize(Eigen::MatrixXd h,
	double in_x, double in_y,
	double in_width, double in_height,
	double& min_x, double& max_x,
	double& min_y, double& max_y)
{
	Eigen::Vector3d p[4] = {
		{ in_x, in_y, 1.0 },
		{ in_width, in_y, 1.0 },
		{ in_width, in_height, 1.0 },
		{ in_x, in_height, 1.0 } };


	Eigen::Vector3d r[4];

	for (int i = 0; i < 4; ++i)
	{
		r[i] = h * p[i];
		r[i] = Eigen::Vector3d(r[i].x() / r[i].z(), r[i].y() / r[i].z(), 1.0);

		std::cout << std::fixed
			<< "p[" << i << "] : (" << p[i].x() << ", " << p[i].y() << ")    "
			<< "r[" << i << "] : (" << r[i].x() << ", " << r[i].y() << ")" << std::endl;
	}

	min_x = r[0].x();
	max_x = r[0].x();

	min_y = r[0].y();
	max_y = r[0].y();

	for (int i = 1; i < 4; ++i)
	{
		if (r[i].x() < min_x)
			min_x = r[i].x();

		if (r[i].y() < min_y)
			min_y = r[i].y();

		if (r[i].x() > max_x)
			max_x = r[i].x();

		if (r[i].y() > max_y)
			max_y = r[i].y();
	}


	std::cout << "\n\nMinMax : (" << min_x << "," << min_y << ") (" << max_x << ", " << max_y << ")" << std::endl;
	std::cout << "Image Size: " << max_x - min_x << ", " << max_y - min_y << std::endl;
}


static void projectImages(const Eigen::MatrixXd H, const std::pair<QImage, QImage>& image, QImage& output)
{
	// composing both images
	Eigen::Vector2d img_size(image.second.width(), image.second.height());
	//Eigen::Vector3d new_size = H.inverse() * img_size.homogeneous();
	//new_size /= new_size[2];

	//std::cout << "original size:  " << img_size.transpose() << std::endl;
	//std::cout << "new size:       " << new_size.transpose() << std::endl;



	double xmin = 0;
	double xmax = 0;
	double ymin = 0;
	double ymax = 0;
	computImageSize(H.inverse(), 0, 0, image.second.width(), image.second.height(), xmin, xmax, ymin, ymax);

	

	double aspect = (xmax - xmin) / (ymax - ymin);

	output = QImage(xmax, ymax - ymin, image.second.format());
	output.fill(Qt::GlobalColor::black);

	std::cout << "Output Size:       " << output.width() << ", " << output.height() << std::endl;

	//double dx = (xmax - xmin) / double(output.width());
	//double dy = (ymax - ymin) / double(output.height());
	double dx = (xmax - xmin) / double(image.first.width());
	double dy = (ymax - ymin) / double(image.first.height());
	std::cout << std::fixed << "dx, dy: " << dx << ", " << dy << std::endl;


	const QImage& input = image.first;
	for (int x = 0; x < input.width(); ++x)
	{
		for (int y = 0; y < input.height(); ++y)
		{
			Eigen::Vector3d p = Eigen::Vector3d(x, y, 1.0);
			//Eigen::Vector3d p = H * Eigen::Vector3d(xmin + x * dx, ymin + y * dy, 1.0);
			p /= p[2];

			if (p.x() > -1 && p.y() > -1
				&& p.x() < image.first.width()
				&& p.y() < image.first.height())
			{
				output.setPixel(x, y, image.first.pixel(p.x(), p.y()));
			}
		}
	}

	for (int x = 0; x < output.width(); ++x)
	{
		for (int y = 0; y < output.height(); ++y)
		{
			//Eigen::Vector3d p = Eigen::Vector3d(xmin + x * dx, ymin + y * dy, 1.0);
			Eigen::Vector3d p = H * Eigen::Vector3d(x, y, 1.0);
			p /= p[2];

			if (p.x() > -1 && p.y() > -1
				&& p.x() < image.first.width()
				&& p.y() < image.first.height())
			{
				output.setPixel(x, y, image.second.pixel(p.x(), p.y()));
			}
		}
	}

#if 0
#if 0
	const QImage& input = image.first;
	for (int x = 0; x < output.width(); ++x)
	{
		for (int y = 0; y < output.height(); ++y)
		{
			Eigen::Vector3d p = H * Eigen::Vector3d(xmin + x * dx, ymin + y * dy, 1.0);
			p /= p[2];

			if (p.x() > -1 && p.y() > -1
				&& p.x() < image.first.width()
				&& p.y() < image.first.height())
			{
				output.setPixel(x, y, image.first.pixel(p.x(), p.y()));
			}


			p = H * Eigen::Vector3d(x, y, 1.0);
			p /= p[2];

			if (p.x() > -1 && p.y() > -1
				&& p.x() < image.second.width()
				&& p.y() < image.second.height())
			{
				output.setPixel(x, y, image.second.pixel(p.x(), p.y()));
			}
		}
	}


#else
	for (int x = 0; x < image.first.width(); ++x)
	{
		for (int y = 0; y < image.first.height(); ++y)
		{
			output.setPixel(x, y, image.first.pixel(x, y));
		}
	}



	for (int x = xmin; x < output.width(); ++x)
	{
		for (int y = ymin; y < output.height(); ++y)
		{
			Eigen::Vector3d p = H * Eigen::Vector3d(x, y, 1.0);
			p /= p[2];

			if (p.x() > -1 && p.y() > -1
				&& p.x() < image.second.width()
				&& p.y() < image.second.height())
			{
				output.setPixel(x, y, image.second.pixel(p.x(), p.y()));
			}
		}
	}
#endif
#endif
}




#endif	// __IMAGE_HELPER_H__
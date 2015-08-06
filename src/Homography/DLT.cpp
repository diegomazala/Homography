#include "DLT.h"
#include <iostream>


DLT::DLT() :isNormalized(true)
{

}

DLT::DLT(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts) : isNormalized(true), points(pts)
{

}

DLT::~DLT()
{
}

void DLT::setNormalized(bool enable)
{
	isNormalized = enable;
}

void DLT::reset()
{
	H = Eigen::Matrix4d::Identity();
	points.clear();
	error = std::make_pair(0.0, 0.0);
}

void DLT::addPoint(const std::pair<Eigen::Vector2d, Eigen::Vector2d>& pt)
{
	points.push_back(pt);
}





Eigen::MatrixXd DLT::computeHomography()
{
	if (isNormalized)
	{
		normalizationTransform = normalizePoints(points, normalizedPoints);
		Eigen::MatrixXd Hn = computeHomography(normalizedPoints);
		this->H = denormalizeH(Hn, normalizationTransform);
	}
	else
	{
		this->H = computeHomography(points);
	}

	return this->H;
}




std::pair<double, double> DLT::computeError()
{
	error = computeError(this->H, this->points);
	return error;
}



Eigen::MatrixXd DLT::denormalizeH(const Eigen::Matrix3d& H, const std::pair<Eigen::Matrix3d, Eigen::Matrix3d>& normalizationTransform)
{
	return normalizationTransform.second.inverse() * H * normalizationTransform.first;
}



std::pair<Eigen::Matrix3d, Eigen::Matrix3d> DLT::buildScaleMatrix(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points,
	const std::pair<Eigen::Matrix3d, Eigen::Matrix3d>& t)
{
	std::pair<double, double> sum_dist(0.0, 0.0);

	// translate points moving the centroid to origin
	// compute the sum of distances
	for (const auto p : points)
	{
		Eigen::Vector3d p1 = t.first * p.first.homogeneous();
		Eigen::Vector3d p2 = t.second * p.second.homogeneous();

		p1 -= Eigen::Vector3d(0, 0, 1);
		p2 -= Eigen::Vector3d(0, 0, 1);

		sum_dist.first += p1.norm();
		sum_dist.second += p2.norm();
	}


	// compute average distance from origin
	std::pair<double, double> average(sum_dist.first / points.size(), sum_dist.second / points.size());
	std::pair<double, double> scale(std::sqrt(2.0) / average.first, std::sqrt(2.0) / average.second);


	// Build T to translate points to origin and scale then to have the average distance from origin equal to sqrt(2.0)
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d>	s(Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Identity());


	// Scale matrix
	s.first(0, 0) = s.first(1, 1) = scale.first;
	s.second(0, 0) = s.second(1, 1) = scale.second;

	return s;
}



std::pair<Eigen::Matrix3d, Eigen::Matrix3d>	DLT::buildTranslationMatrix(const std::pair<Eigen::Vector2d, Eigen::Vector2d>& centroid)
{
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d>	Tt(Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Identity());
	Tt.first(0, 2) = -centroid.first(0);
	Tt.first(1, 2) = -centroid.first(1);
	Tt.second(0, 2) = -centroid.second(0);
	Tt.second(1, 2) = -centroid.second(1);
	return Tt;
}



std::pair<Eigen::Vector2d, Eigen::Vector2d> DLT::findCentroid(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points)
{
	// computing centroid of the points
	int point_count = (int)points.size();
	std::pair<Eigen::Vector2d, Eigen::Vector2d> sum(Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero());
	for (int i = 0; i < point_count; ++i)
	{
		sum.first += points[i].first;
		sum.second += points[i].second;
	}

	std::pair<Eigen::Vector2d, Eigen::Vector2d>	centroid = std::make_pair(sum.first / point_count, sum.second / point_count);

	//std::cout << "[Info]  Centroid: " << std::endl
	//	<< centroid.first.transpose() << std::endl
	//	<< centroid.second.transpose() << std::endl << std::endl;

	return centroid;
}





std::pair<Eigen::Matrix3d, Eigen::Matrix3d> DLT::normalizePoints(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& in_points,
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& out_points)
{
	std::pair<double, double> sum_dist(0.0, 0.0);
	std::pair<Eigen::Vector2d, Eigen::Vector2d>	centroid = findCentroid(in_points);
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d>	Tt = buildTranslationMatrix(centroid);
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d>	Ts = buildScaleMatrix(in_points, Tt);

	std::pair<Eigen::Matrix3d, Eigen::Matrix3d> T;
	T.first = Ts.first * Tt.first;
	T.second = Ts.second * Tt.second;

	sum_dist.first = 0;
	sum_dist.second = 0;

	// normalize point set
	out_points.clear();
	for (const auto p : in_points)
	{
		Eigen::Vector3d p1 = T.first * p.first.homogeneous();
		Eigen::Vector3d p2 = T.second * p.second.homogeneous();

		p1 /= p1[2];
		p2 /= p2[2];

		p1 -= Eigen::Vector3d(0, 0, 1);
		p2 -= Eigen::Vector3d(0, 0, 1);

		out_points.push_back(std::make_pair(Eigen::Vector2d(p1.x(), p1.y()), Eigen::Vector2d(p2.x(), p2.y())));

		sum_dist.first += p1.norm();
		sum_dist.second += p2.norm();

		//std::cout << std::fixed << p1.transpose() << "\t\t" << p2.transpose() << std::endl;
	}

	std::cout << std::endl;


	//std::cout << "[Info]  Average Distance: " << sum_dist.first / in_points.size() << " , " << sum_dist.second / in_points.size() << std::endl;


	//std::cout << "[Info]  T1: " << std::endl << T.first << std::endl << std::endl;
	//std::cout << "[Info]  T2: " << std::endl << T.second << std::endl << std::endl;

	return T;
}


Eigen::MatrixXd DLT::computeHomography(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points)
{
	assert(points.size() == 4);

	Eigen::MatrixXd A(2 * points.size(), 9);	// 8 x 9
	int i = 0, j = 0;
	for (const auto p : points)
	{
		A(i, 0) = -p.first.x();
		A(i, 1) = -p.first.y();
		A(i, 2) = -1.0;
		A(i, 3) = 0.0;
		A(i, 4) = 0.0;
		A(i, 5) = 0.0;
		A(i, 6) = p.first.x() * p.second.x();
		A(i, 7) = p.first.x() * p.second.y();
		A(i, 8) = p.second.x();
		++i;
		A(i, 0) = 0.0;
		A(i, 1) = 0.0;
		A(i, 2) = 0.0;
		A(i, 3) = -p.first.x();
		A(i, 4) = -p.first.y();
		A(i, 5) = -1.0;
		A(i, 6) = p.first.x() * p.second.y();
		A(i, 7) = p.first.y() * p.second.y();
		A(i, 8) = p.second.y();
		++i;
	}

	//std::cout << "A: " << std::endl << A << std::endl << std::endl;

	Eigen::MatrixXd kernel = A.fullPivLu().kernel();

	Eigen::MatrixXd H(3, 3);
	H(0, 0) = kernel(0);
	H(0, 1) = kernel(1);
	H(0, 2) = kernel(2);
	H(1, 0) = kernel(3);
	H(1, 1) = kernel(4);
	H(1, 2) = kernel(5);
	H(2, 0) = kernel(6);
	H(2, 1) = kernel(7);
	H(2, 2) = kernel(8);

	return H;
}




std::pair<double, double> DLT::computeError(const Eigen::MatrixXd H, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points)
{
	std::pair<double, double> error(0.0, 0.0);

	for (auto p : points)
	{
		Eigen::Vector3d HL = H * p.first.homogeneous();
		HL /= HL[2];
		Eigen::Vector3d HR = H.inverse() * p.second.homogeneous();
		HR /= HR[2];

		error.first += std::pow((p.first - Eigen::Vector2d(HR.x(), HR.y())).norm(), 2);
		error.second += std::pow((p.second - Eigen::Vector2d(HL.x(), HL.y())).norm(), 2);
	}


	std::cout << std::fixed << "[Info]  Projection Error  : " << error.first << ", " << error.second << std::endl;

	return error;
}



bool DLT::operator < (DLT const &other)
{
	return (this->error.first + this->error.second) < (other.error.first + other.error.second);
}
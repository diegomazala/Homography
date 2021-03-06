#include "Points.h"
#include <iostream>
#include <fstream>
#include <random>

Points::Points()
{
}


Points::~Points()
{
}

void Points::clear()
{
	points.clear();
}


std::size_t Points::count() const
{
	return points.size();
}

std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& Points::getPointArray()
{
	return this->points;
}

const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& Points::getPointArray() const
{
	return this->points;
}

std::pair<Eigen::Vector2d, Eigen::Vector2d>& Points::operator[](std::size_t index)
{
	return points[index];
};

const std::pair<Eigen::Vector2d, Eigen::Vector2d>& Points::operator[](std::size_t index) const
{
	return points[index];
};


bool Points::readFromFile(const std::string filePath)
{
	return readFromFile(filePath, points);
}




bool Points::readFromFile(const std::string filePath, std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points)
{
	std::ifstream file(filePath);
	if (!file.is_open())
	{
		std::cerr << "Error: Could not open points file: " << filePath << std::endl;
		return false;
	}

	int point_count = 0;
	std::string str;
	char c;
	file >> str >> point_count;

	if (point_count < 1)
	{
		std::cerr << "Error: Could not read points file: " << filePath << ". Point count = " << point_count << std::endl;
		return false;
	}

	points.clear();
	int i = 0;
	while (i++ < point_count && file.good())
	{
		double x0, y0, x1, y1;
		file >> c >> x0 >> c >> y0 >> c >> c >> x1 >> c >> y1 >> c;
		points.push_back(std::make_pair(Eigen::Vector2d(x0, y0), Eigen::Vector2d(x1, y1)));
	}

	if (!file.good())
	{
		std::cerr << "Error: Reading points file failed at line: " << i + 1 << std::endl;
		return false;
	}

	file.close();
	return true;
}



void Points::projectPoints(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts_src,
							std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts_dst,
							const Eigen::Matrix3d& H)
{
	pts_dst.clear();
	for (const auto pt : pts_src)
	{
		const Eigen::Vector2d p1 = pt.first;
		Eigen::Vector3d hp1 = H * p1.homogeneous();
		hp1 /= hp1[2];

		const Eigen::Vector2d p2 = pt.second;
		Eigen::Vector3d hp2 = H.inverse() * p2.homogeneous();
		hp2 /= hp2[2];

		pts_dst.push_back(std::make_pair(Eigen::Vector2d(hp1.x(), hp1.y()), Eigen::Vector2d(hp2.x(), hp2.y())));
	}
}



void Points::projectPoints(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts_src,
							std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& pts_dst,
							const Eigen::Matrix3d& H)
{
	pts_dst.clear();
	for (const auto pt : pts_src)
	{
		const Eigen::Vector2d p1 = pt.first;
		Eigen::Vector3d hp1 = H * p1.homogeneous();

		const Eigen::Vector2d p2 = pt.second;
		Eigen::Vector3d hp2 = H.inverse() * p2.homogeneous();

		pts_dst.push_back(std::make_pair(hp1, hp2));
	}
}




void Points::exportPointCorrespondencies(const std::string& filename, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D)
{
	std::ofstream outFile;
	outFile.open(filename);

	outFile << "POINT_COUNT: " << points2D.size() << std::endl;

	for (auto p : points2D)
	{
		outFile
			<< std::fixed
			<< '[' << p.first.x() << ", " << p.first.y()
			<< "]\t[" << p.second.x() << ", " << p.second.y()
			<< ']' << std::endl;
	}

	outFile.close();
}


bool Points::importPointCorrespondencies(const std::string filePath, std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points)
{
	std::ifstream file(filePath);
	if (!file.is_open())
	{
		std::cerr << "Error: Could not open points file: " << filePath << std::endl;
		return false;
	}

	int point_count = 0;
	std::string str;
	char c;
	file >> str >> point_count;

	if (point_count < 1)
	{
		std::cerr << "Error: Could not read points file: " << filePath << ". Point count = " << point_count << std::endl;
		return false;
	}

	points.clear();
	int i = 0;
	while (i++ < point_count && file.good())
	{
		double x0, y0, x1, y1;
		file >> c >> x0 >> c >> y0 >> c >> c >> x1 >> c >> y1 >> c;
		points.push_back(std::make_pair(Eigen::Vector2d(x0, y0), Eigen::Vector2d(x1, y1)));
	}

	if (!file.good())
	{
		std::cerr << "Error: Reading points file failed at line: " << i + 1 << std::endl;
		return false;
	}

	file.close();
	return true;
}

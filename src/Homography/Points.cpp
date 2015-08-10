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



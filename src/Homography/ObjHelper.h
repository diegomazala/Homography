#ifndef __OBJ_HELPER_H__
#define __OBJ_HELPER_H__

#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include "Triangulation.h"



static bool readPointsFromObj(const std::string& filename, std::vector<Eigen::Vector3d>& points3D, int max_point_count = INT_MAX)
{
	std::ifstream inFile;
	inFile.open(filename);

	if (!inFile.is_open())
	{
		std::cerr << "Error: Could not open obj input file: " << filename << std::endl;
		return false;
	}

	points3D.clear();

	int i = 0;
	while (inFile)
	{
		std::string str;

		if (!std::getline(inFile, str))
		{
			std::cerr << "Error: Problems when reading obj file: " << filename << std::endl;
			return false;
		}

		if (str[0] == 'v')
		{
			std::stringstream ss(str);
			std::vector <std::string> record;

			char c;
			double x, y, z;
			ss >> c >> x >> y >> z;

			Eigen::Vector3d p(x, y, z);
			points3D.push_back(p);
		}

		if (i++ > max_point_count)
			break;
	}

	inFile.close();
	return true;
}

static void exportObj(const std::string& filename, const std::vector<Eigen::Vector3d>& points3D)
{
	std::ofstream file;
	file.open(filename);
	for (const auto X : points3D)
	{
		file << "v " << X.transpose() << std::endl;
	}
	file.close();
}

static void exportObj(const std::string& filename, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D, const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& P)
{
	std::ofstream file;
	file.open(filename);
	for (const auto x : points2D)
	{
		const Eigen::VectorXd& X = Triangulation::solve(P, x);
		file << "v " << X.transpose() << std::endl;
	}
	file.close();
}

static void exportObj(const std::string& filename1, const std::string& filename2, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D)
{
	std::ofstream file;
	file.open(filename1);
	for (const auto x : points2D)
	{
		file << "v " << x.first.homogeneous().transpose() << std::endl;
	}
	file.close();
	file.open(filename2);
	for (const auto x : points2D)
	{
		file << "v " << x.second.homogeneous().transpose() << std::endl;
	}
	file.close();
}



static void exportCubeObj(const std::string& filename, const std::vector<Eigen::Vector3d>& points3D)
{
	std::ofstream file;
	file.open(filename);
	for (const auto X : points3D)
	{
		file << "v " << X.transpose() << std::endl;
	}

	file
		<< "f 1 2 3 4" << std::endl
		<< "f 8 7 6 5" << std::endl
		<< "f 5 6 2 1" << std::endl
		<< "f 4 3 7 8" << std::endl
		<< "f 2 6 7 3" << std::endl
		<< "f 5 1 4 8" << std::endl;
	file.close();
}

static void exportCubeObj(const std::string& filename, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D, const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& P)
{
	std::ofstream file;
	file.open(filename);
	for (const auto x : points2D)
	{
		const Eigen::VectorXd& X = Triangulation::solve(P, x);
		file << "v " << X.transpose() << std::endl;
	}

	file
		<< "f 1 2 3 4" << std::endl
		<< "f 8 7 6 5" << std::endl
		<< "f 5 6 2 1" << std::endl
		<< "f 4 3 7 8" << std::endl
		<< "f 2 6 7 3" << std::endl
		<< "f 5 1 4 8" << std::endl;
	file.close();
}



bool project3DPointsFile(const std::string& filename_in, const std::string& filename_out, const Eigen::MatrixXd& mat, int max_point_count = INT_MAX)
{
	std::ifstream inFile;
	inFile.open(filename_in);

	if (!inFile.is_open())
	{
		std::cerr << "Error: Could not open obj input file: " << filename_in << std::endl;
		return false;
	}

	std::ofstream outFile;
	outFile.open(filename_out);

	assert(mat.rows() == 3 && mat.cols() == 4);

	int i = 0;
	while (inFile)
	{
		std::string str;

		if (!std::getline(inFile, str))
		{
			std::cerr << "Error: Problems when reading obj file: " << filename_in << std::endl;
			return false;
		}


		if (str[0] == 'v')
		{
			std::stringstream ss(str);
			std::vector <std::string> record;

			char c;
			double x, y, z;
			ss >> c >> x >> y >> z;

			Eigen::Vector4d ptX(x, y, z, 1);
			Eigen::Vector3d pt = mat * ptX;
			pt /= pt[2];
			
			outFile << std::fixed << "v " << pt.transpose() << std::endl;
		}

		if (i++ > max_point_count)
			break;

	}

	inFile.close();
	outFile.close();

	return true;
}



#endif //__DLT_H__

#include "RansacDLT.h"
#include <random>
#include <iostream>




DLT RansacDLT::solve(const Points& points)
{
	std::vector<DLT> dltArray;

	int inliersPercentageAchieved = 0;
	int inliersPercentageAccepted = 99;

	//
	// Run the iteration to find the best DLT, i.e, with the biggest number of inliers
	//
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

		inliersPercentageAchieved = (int)(double(inliers) / double(points.count()) * 100);
	}

	std::cout 
		<< std::endl
		<< "[Info]  Number of Iterations : " << i << std::endl
		<< "[Info]  Inliers Acheived     : " << inliersPercentageAchieved << "%" << std::endl << std::endl;

	return dltArray.back();
}




std::vector<int> RansacDLT::random4Indices(int min, int max)
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
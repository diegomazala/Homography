#include "RansacDLT.h"
#include <random>
#include <iostream>




DLT RansacDLT::solve(const Points& points)
{
	double s = 4.0;
	double p = 0.99;
	double N = (double)points.count();
	double T = 0.0;
	double E = (double)points.count();
	std::vector<DLT> dltArray;

	int inliersPercentageAchieved = 0;
	int inliersPercentageAccepted = 99;

	//
	// Run the iteration to find the best DLT, i.e, with the biggest number of inliers
	//
	int i = 0;
	//while (inliersPercentageAchieved < inliersPercentageAccepted && i++ < points.count())
	while (i < N && i++ < points.count())
	{
		dltArray.push_back(DLT());
		DLT& dlt = dltArray.back();

		//
		// Generating 4 indices for points in order to compute homography using these points
		// 
		std::cout << std::endl
			<< "[Info]  Iteration = " << i << std::endl
			<< "[Info]  Random 4 indices = ";
		std::vector<int> indices = randomIndices(4, 0, (int)points.count() - 1);

		// copy points from original array to the 4-array to be used for homography
		for each (auto i in indices)
		{
			dlt.addPoint(points[i]);
			std::cout << i << "  ";
		}

		Eigen::MatrixXd H = dlt.computeHomography();
		dlt.computeGeometricError();


		//
		// project all points and count number of inliers and outliers
		//
		int inliers = dlt.computeInliers(points.getPointArray());

		// Compute error for this iteration and update the global error 
		double Ei = 1.0 - (double(inliers) / double(points.count()));
		if (Ei < E)
			E = Ei;

		N = std::log(1.0 - p) / std::log( 1.0 - std::pow( 1.0 - E, s)  );

		T = (1.0 - E) * points.count();

		inliersPercentageAchieved = (int)(double(inliers) / double(points.count()) * 100);

		std::cout 
			<< std::fixed
			<< "[Info]  Inliers Achieved  : " << inliersPercentageAchieved << "%  ==> " << inliers << " of " << points.count() << std::endl
			<< "[Info]  E, N, T           : " << E << ", " << N << ", " << T << std::endl << std::endl;
	}

	DLT* dltConsensus = &dltArray.back();

	if (inliersPercentageAchieved < inliersPercentageAccepted)
	{
		std::sort(dltArray.begin(), dltArray.end());
		dltConsensus = &dltArray.front();
		inliersPercentageAchieved = (int)(double(dltConsensus->getInliersCount()) / double(points.count()) * 100);
	}

	std::cout
		<< std::endl
		<< "[Info]  -- Consensus:" << std::endl
		<< "[Info]  Number of Iterations : " << i << std::endl
		<< "[Info]  Inliers Achieved     : " << inliersPercentageAchieved << "%" << std::endl
		<< "[Info]  E, N, T              : " << E << ", " << N << ", " << T << std::endl << std::endl;

	return *dltConsensus;
}




std::vector<int> RansacDLT::randomIndices(int count, int min, int max)
{
	std::random_device rd;								// obtain a random number from hardware
	std::default_random_engine eng(rd());				// seed the generator
	std::uniform_int_distribution<int> distr(min, max); // define the range

	std::vector<int> indices;

	for (int i = 0; i < count; ++i)
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
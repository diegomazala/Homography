#include "GaussNewton.h"
#include "Points.h"
#include <iostream>


Eigen::MatrixXd GaussNewton::computeDeltaH(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const Eigen::MatrixXd& H)
{
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> hp;
	Points::projectPoints(pts, hp, H);

	//double sum_error = GaussNewton::getSumError(pts, H);

	Eigen::MatrixXd J = GaussNewton::buildMatrixJ(pts);

	Eigen::MatrixXd A = J.transpose() * J;

	Eigen::VectorXd fx = GaussNewton::buildFx(pts, H);

	//sum_error = GaussNewton::getSumError(pts, H);
	//std::cout << std::endl << "fx" << std::endl << fx << std::endl << std::endl;
	//std::cout << "sum_error: " << sum_error << std::endl;

	Eigen::MatrixXd b = -J.transpose() * fx;
	Eigen::MatrixXd x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

	//std::cout << "fx: " << fx.rows() << ", " << fx.cols() << std::endl;
	//std::cout << "b: " << b.rows() << ", " << b.cols() << std::endl;
	//std::cout << "x: " << x.rows() << ", " << x.cols() << std::endl;

	Eigen::MatrixXd dH(3, 3);
	dH(0, 0) = x(0);
	dH(0, 1) = x(1);
	dH(0, 2) = x(2);
	dH(1, 0) = x(3);
	dH(1, 1) = x(4);
	dH(1, 2) = x(5);
	dH(2, 0) = x(6);
	dH(2, 1) = x(7);
	dH(2, 2) = x(8);

	return dH;
}


double GaussNewton::getSumError(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const Eigen::Matrix3d& H)
{
	double sum_error = 0.0;
	for (const auto pt : pts)
	{
		const Eigen::Vector2d p1 = pt.first;
		Eigen::Vector3d hp1 = H * p1.homogeneous();
		hp1 /= hp1[2];

		const Eigen::Vector2d p2 = pt.second;
		Eigen::Vector3d hp2 = H.inverse() * p2.homogeneous();
		hp2 /= hp2[2];

		double error = std::pow((p1 - Eigen::Vector2d(hp2.x(), hp2.y())).norm(), 2);
		//double error = std::pow((p2 - Eigen::Vector2d(hp1.x(), hp1.y())).norm(), 2);
		//std::cout << std::fixed << "error     : " << error << "   " << p1.transpose() << "   " << hp2.transpose() << std::endl;

		sum_error += error;
	}
	return sum_error;
}


double GaussNewton::getSumError(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& hp)
{
	double sum_error = 0.0;
	for (int i = 0; i < pts.size(); ++i)
	{
		double error = std::pow((pts[i].first - hp[i].second).norm(), 2);
		//std::cout << std::fixed << "error     : " << error << "   " << pts[i].first.transpose() << "   " << hp[i].second.transpose() << std::endl;
		sum_error += error;
	}
	return sum_error;
}



Eigen::MatrixXd GaussNewton::buildMatrixJ(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts)
{
	int i = 0;
	Eigen::MatrixXd J(2 * pts.size(), 9);		// 8 x 9, 2n x 9
	for (const auto p : pts)
	{
		J(i, 0) = p.first.x();
		J(i, 1) = p.first.y();
		J(i, 2) = 1.0;
		J(i, 3) = 0.0;
		J(i, 4) = 0.0;
		J(i, 5) = 0.0;
		J(i, 6) = -p.second.x() * p.first.x();
		J(i, 7) = -p.second.x() * p.first.y();
		J(i, 8) = -p.second.x();
		++i;
		J(i, 0) = 0.0;
		J(i, 1) = 0.0;
		J(i, 2) = 0.0;
		J(i, 3) = p.first.x();
		J(i, 4) = p.first.y();
		J(i, 5) = 1.0;
		J(i, 6) = -p.second.y() * p.first.x();
		J(i, 7) = -p.second.y() * p.first.y();
		J(i, 8) = -p.second.y();
		++i;
	}
	
	return J;
}

Eigen::VectorXd GaussNewton::buildFx(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const Eigen::MatrixXd& H)
{
	Eigen::VectorXd fx(2 * pts.size());	// 8
	int p = 0;
	for (int i = 0; i < pts.size(); ++i)
	{
		const Eigen::Vector2d p1 = pts[i].second;
		Eigen::Vector3d hp3d = H.inverse() * p1.homogeneous();
		hp3d /= hp3d[2];
		Eigen::Vector2d hp(hp3d.x(), hp3d.y());

#if 1
		fx(p++) = std::pow((pts[i].first.x() - hp.x()), 2);
		fx(p++) = std::pow((pts[i].first.y() - hp.y()), 2);

#else
		double error = std::pow((pts[i].first - hp).norm(), 2);
		fx(i) = error;
		//std::cout << std::fixed << "error  fx : " << fx(i) << "   " << pts[i].first.transpose() << "   " << hp.transpose() << std::endl;
#endif
	}

	return fx;
}
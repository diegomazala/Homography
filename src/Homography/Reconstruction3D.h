#ifndef __RECONSTRUCTION_3D_H__
#define __RECONSTRUCTION_3D_H__

#include <Eigen/Dense>
#include <vector>


class Reconstruction3D
{
public:
	
	static Eigen::MatrixXd buildMatrixA(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts);

	static Eigen::MatrixXd applyConstraint(const Eigen::MatrixXd& F);

	static Eigen::MatrixXd computeF(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts);
	
	static Eigen::MatrixXd computeE(const Eigen::MatrixXd& K, const Eigen::MatrixXd& F);
	
	static Eigen::MatrixXd computeP(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
									const Eigen::MatrixXd& E);

	static bool checkP(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
						const Eigen::MatrixXd& P0,
						const Eigen::MatrixXd& P1);

	static Eigen::VectorXd triangulation(
										const Eigen::MatrixXd& P0,
										const Eigen::MatrixXd& P1,
										const Eigen::VectorXd& p0,
										const Eigen::VectorXd& p1);
};



#endif //__RECONSTRUCTION_3D_H__


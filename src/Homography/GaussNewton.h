#ifndef __GAUSS_NEWTON_H__
#define __GAUSS_NEWTON_H__

#include <Eigen/Dense>
#include <vector>


class GaussNewton
{
public:

	static double solve(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, 
						const Eigen::MatrixXd& H, 
						int max_iterations,
						Eigen::MatrixXd& H_out);



	static Eigen::MatrixXd computeDeltaH(	const std::vector<std::pair<Eigen::Vector2d, 
											Eigen::Vector2d>>& pts, 
											const Eigen::MatrixXd& H);


	static double getSumError(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, 
								const Eigen::Matrix3d& H);


	static double getSumError(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, 
								const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& hp);


	static Eigen::MatrixXd buildMatrixJ(const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& pts);

	
	static Eigen::VectorXd buildFx(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, 
									const Eigen::MatrixXd& H);

};



#endif //__GAUSS_NEWTON_H__

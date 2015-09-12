#ifndef __RECONSTRUCTION_3D_H__
#define __RECONSTRUCTION_3D_H__

#include <Eigen/Dense>
#include <vector>


class Reconstruction3D
{
public:
	
	static Eigen::MatrixXd buildMatrixA(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts);

	static Eigen::MatrixXd applyConstraint(const Eigen::MatrixXd& inputMat);

	static Eigen::MatrixXd computeF(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts);

	static Eigen::MatrixXd computeE(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& K, const Eigen::MatrixXd& F);

	static Eigen::MatrixXd computeK(const Eigen::MatrixXd& F);
		
	static Eigen::Matrix3d computeEpipoleMat(const Eigen::Matrix3d& F);
	static Eigen::Vector3d computeEpipole(const Eigen::Matrix3d& F);

	static Eigen::MatrixXd computeP(const Eigen::MatrixXd& F);

	static Eigen::MatrixXd computeP(const Eigen::MatrixXd& F, const Eigen::MatrixXd& eMat);

	static Eigen::MatrixXd computeP(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
		const Eigen::MatrixXd& E);

	static void computeP(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
		const Eigen::MatrixXd& E, std::vector<Eigen::MatrixXd>& P_solutions);

	static bool checkP(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
		const Eigen::MatrixXd& P0,
		const Eigen::MatrixXd& P1);


	static double computeError(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const Eigen::MatrixXd& F);

	static void solveCS(float &c, float &s, float a, float b);
	static void RQdecomposition(Eigen::MatrixXd A, Eigen::Matrix3d &R, Eigen::Matrix3d &Q);
	static void QRdecomposition(Eigen::MatrixXd A, Eigen::Matrix3d &R, Eigen::Matrix3d &Q);
};



#endif //__RECONSTRUCTION_3D_H__


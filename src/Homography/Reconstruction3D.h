#ifndef __RECONSTRUCTION_3D_H__
#define __RECONSTRUCTION_3D_H__

#include <Eigen/Dense>
#include <vector>



struct ReconstructionDLT
{
	ReconstructionDLT(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
		const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& _K);

	void solve(double outlierThreshold = 30.0);

	bool operator < (ReconstructionDLT const &other);

	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> points2D;
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> points2DNorm;
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> K;
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> P;
	double error;
	int inliersCount;
};

class Reconstruction3D
{
public:


	static ReconstructionDLT solve(
								const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts,
								const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& K,
								int maxIterations = 500);

	
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

	static Eigen::MatrixXd selectBestP(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points2D,
										const std::pair<Eigen::Matrix3d, Eigen::Matrix3d>& K,
										const std::vector<Eigen::MatrixXd>& P_solutions);
	
	static double pointLineDistance(Eigen::Vector2d point, Eigen::Vector3d line);

	static double computeError(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const Eigen::MatrixXd& F);
	//static double computeGeometricError(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& P, double outlierThreshold, int& inliers);

	static void solveCS(float &c, float &s, float a, float b);
	static void RQdecomposition(Eigen::MatrixXd A, Eigen::Matrix3d &R, Eigen::Matrix3d &Q);
	static void QRdecomposition(Eigen::MatrixXd A, Eigen::Matrix3d &R, Eigen::Matrix3d &Q);
};





#endif //__RECONSTRUCTION_3D_H__


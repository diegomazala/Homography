#ifndef __TRIANGULATION_H__
#define __TRIANGULATION_H__

#include <Eigen/Dense>
#include <vector>


class Triangulation
{
public:
	
	Triangulation(	const Eigen::MatrixXd& P0,
					const Eigen::MatrixXd& P1,
					const Eigen::VectorXd& x0,
					const Eigen::VectorXd& x1);


	Triangulation(	const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& _P,
					const std::pair<Eigen::VectorXd, Eigen::VectorXd>& _x);


	Eigen::Vector3d solve() const;


	static Eigen::Vector3d solve(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& _P,
								const std::pair<Eigen::VectorXd, Eigen::VectorXd>& _x);


	static Eigen::MatrixXd buildMatrixA(const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& _P,
										const std::pair<Eigen::VectorXd, Eigen::VectorXd>& _x);

private:

	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> P;
	std::pair<Eigen::VectorXd, Eigen::VectorXd> x;

};



#endif //__TRIANGULATION_H__


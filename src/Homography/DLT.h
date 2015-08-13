#ifndef __DLT_H__
#define __DLT_H__

#include <Eigen/Dense>
#include <vector>

class DLT
{
public:

	DLT();
	DLT(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts);
	~DLT();

	const Eigen::MatrixXd& getH() const { return H; };
	const std::pair<double, double>& getError() const { return error; };
	int getInliersCount() const;

	void setNormalized(bool enable);

	void reset();
	void addPoint(const std::pair<Eigen::Vector2d, Eigen::Vector2d>& pt);

	
	
	Eigen::MatrixXd computeHomography();
	std::pair<double, double> computeGeometricError();

	int computeInliers(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts, double distance_threshold = std::sqrt(5.99));

	static Eigen::MatrixXd computeHomography(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points);
	
	static std::pair<Eigen::Matrix3d, Eigen::Matrix3d> normalizePoints(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& in_points,
																		std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& out_points);

	static std::pair<Eigen::Matrix3d, Eigen::Matrix3d> buildScaleMatrix(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points,
																		const std::pair<Eigen::Matrix3d, Eigen::Matrix3d>& t);

	static std::pair<Eigen::Matrix3d, Eigen::Matrix3d> buildTranslationMatrix(const std::pair<Eigen::Vector2d, Eigen::Vector2d>& centroid);

	static std::pair<Eigen::Vector2d, Eigen::Vector2d> findCentroid(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points);

	static Eigen::MatrixXd denormalizeH(const Eigen::Matrix3d& H, const std::pair<Eigen::Matrix3d, Eigen::Matrix3d>& normalizationTransform);

	static std::pair<double, double> computeGeometricError(const Eigen::MatrixXd H, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points);


	bool operator < (DLT const &other);

private:

	bool isNormalized;
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> points;
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> normalizedPoints;
	std::pair<Eigen::Matrix3d, Eigen::Matrix3d> normalizationTransform;
	Eigen::MatrixXd H;
	std::pair<double, double> error;
	int inliers;
};


#endif //__DLT_H__
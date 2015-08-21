#ifndef __POINTS_H__
#define __POINTS_H__

#include <Eigen/Dense>
#include <vector>
#include <string>

class Points
{
public:

	Points();
	~Points();

	bool readFromFile(const std::string filePath);

	void clear();

	std::size_t count() const;

	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& getPointArray();
	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& getPointArray() const;


	static void projectPoints(	const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts_src,
								std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& pts_dst,
								const Eigen::Matrix3d& H);
	
	std::pair<Eigen::Vector2d, Eigen::Vector2d>&		operator[](std::size_t idx);

	const std::pair<Eigen::Vector2d, Eigen::Vector2d>&	operator[](std::size_t idx) const;


	static bool readFromFile(const std::string filePath, std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& points);

private:

	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> points;
};



#endif //__POINTS_H__
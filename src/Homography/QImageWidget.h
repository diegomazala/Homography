

#ifndef QImageWidget_H
#define QImageWidget_H

#include <QLabel>
#include <Eigen/Dense>


class QImageWidget : public QLabel
{
    Q_OBJECT

public:
    QImageWidget(QWidget* parent = nullptr);
    bool loadFile(const QString &);
	void setImage(const QImage& image);

	void setPoints(const std::vector<Eigen::Vector2d> points);

protected:
	virtual void paintEvent(QPaintEvent *);
	virtual void keyReleaseEvent(QKeyEvent *);

	std::vector<QPoint> points;
	int circleSize;
};

#endif

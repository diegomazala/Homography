

#include <QtWidgets>


#include "QImageWidget.h"
#include <iostream>

QImageWidget::QImageWidget(QWidget* parent) : QLabel(parent), circleSize(10)
{
}

void QImageWidget::setPoints(const std::vector<Eigen::Vector2d> pts)
{
	points.clear();
	for (auto p : pts)
		points.push_back(QPoint(p.x(), p.y()));
}

void QImageWidget::setImage(const QImage& image)
{
	if (!image.isNull())
	{
		resize(image.size());
		setPixmap(QPixmap::fromImage(image));
	}
	
}

bool QImageWidget::loadFile(const QString &fileName)
{
    QImage image = QImage(fileName);
    if (image.isNull()) 
	{
        QMessageBox::information(this, QGuiApplication::applicationDisplayName(),
                                 tr("Cannot load %1.").arg(QDir::toNativeSeparators(fileName)));
        setWindowFilePath(QString());
        setPixmap(QPixmap());
        adjustSize();
        return false;
    }

	setImage(image);

    setWindowFilePath(fileName);
    return true;
}



void QImageWidget::paintEvent(QPaintEvent* event)
{
	QLabel::paintEvent(event);

	QPainter painter(this);
	painter.setPen(QColor(255, 255, 0));
	painter.setBrush(Qt::BrushStyle::Dense1Pattern);

	for (auto p : points)
		painter.drawEllipse(p.x() - circleSize / 2, p.y() - circleSize / 2, circleSize, circleSize);
}


void QImageWidget::keyReleaseEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_Q)
		this->close();
}
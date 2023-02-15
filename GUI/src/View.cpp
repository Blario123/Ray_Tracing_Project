#include "View.h"

View::View(QWidget *parent) : QWidget(parent), layout(new QGridLayout), view(new QGraphicsView), scene(new QGraphicsScene), image(new QImage) {
    QSize windowSize = screen()->size();
    scene->addPixmap(QPixmap::fromImage(*image));
	setMinimumSize(windowSize.width()/2, windowSize.height()/2);
	QPalette palette = QPalette(Qt::black);
	setAutoFillBackground(true);
	setPalette(palette);
    view->setScene(scene);
    layout->addWidget(view);
    setLayout(layout);
    emit sizeChanged(width(), height());
}

void View::renderImage() {

}

void View::resizeEvent(QResizeEvent *event) {
    emit sizeChanged(event->size().width(), event->size().height());
    QWidget::resizeEvent(event);
}

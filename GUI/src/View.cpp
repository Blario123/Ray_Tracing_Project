#include "View.h"
#include <QResizeEvent>

View::View(QWidget *parent) : QWidget(parent) {
	setMinimumSize(window()->x()/2, window()->y()/2);
	QPalette palette = QPalette(Qt::black);
	setAutoFillBackground(true);
	setPalette(palette);
}

void View::resizeEvent(QResizeEvent *event) {
    emit sizeChanged(event->size().width(), event->size().height());
//    QResizeEvent(event);
}

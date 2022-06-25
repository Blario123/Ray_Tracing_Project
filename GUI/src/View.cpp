#include "View.h"

View::View(QWidget *parent) : QWidget(parent) {
	setMinimumSize(window()->x()/2, window()->y()/2);
	QPalette palette = QPalette(Qt::black);
	setAutoFillBackground(true);
	setPalette(palette);
}

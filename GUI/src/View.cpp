#include "View.h"

View::View(QWidget *parent) : QWidget(parent) {
	setMinimumSize(1024, 720);
	QPalette palette = QPalette(Qt::black);
	setAutoFillBackground(true);
	setPalette(palette);
}

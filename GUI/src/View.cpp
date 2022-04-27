#include "View.h"

View::View(QWidget *parent) : QWidget(parent) {
	setMinimumSize(100, 100);
	QPalette palette = QPalette(Qt::red);
	setAutoFillBackground(true);
	setPalette(palette);
}

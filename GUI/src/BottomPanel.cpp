#include "BottomPanel.h"

BottomPanel::BottomPanel(QWidget *parent) : QWidget(parent),
											groupBox(new QGroupBox("BottomPanel")),
											layout(new QBoxLayout(QBoxLayout::LeftToRight)) {
	layout->addWidget(groupBox);
	setLayout(layout);
}

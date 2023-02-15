#include "SidePanel.h"

SidePanel::SidePanel(QWidget *parent) : QWidget(parent),
										layout(new QGridLayout),
										groupBox(new QGroupBox("SidePanel")),
										boxLayout(new QGridLayout),
										itemTree(new ItemTree),
										addPushButton(new QPushButton("Add")),
										delPushButton(new QPushButton("Delete")) {
    setMaximumWidth(floor((screen()->size().width() * 0.15625)));
	itemTree->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Preferred);
	for (int i = 0; i < itemTree->columnCount(); i++) {
		itemTree->resizeColumnToContents(i);
	}
	
	boxLayout->addWidget(itemTree, 0, 0, 1, 2);
	boxLayout->addWidget(addPushButton, 1, 0);
	boxLayout->addWidget(delPushButton, 1, 1);
	
	groupBox->setLayout(boxLayout);
	layout->addWidget(groupBox);
	
	connect(addPushButton, &QPushButton::pressed, itemTree, &ItemTree::addPressed);
	connect(delPushButton, &QPushButton::pressed, itemTree, &ItemTree::delPressed);
	setLayout(layout);
}

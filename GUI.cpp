#include <QComboBox>
#include "GUI.h"

GUI::GUI(QWidget *parent) : QWidget(parent),
							gridLayout(new QGridLayout),
							hSplitter(new QSplitter),
							vSplitter(new QSplitter),
							bottomPanel(new BottomPanel),
							fileBar(new FileBar),
							sidePanel(new SidePanel),
							statusBar(new StatusBar),
							view(new View) {
	connect(fileBar->closeAction, &QAction::triggered, this, &QWidget::close);
	connect(fileBar->aboutAction, &QAction::triggered, this, &QApplication::aboutQt);
	connect(fileBar->saveAction, &QAction::triggered, fileBar, &FileBar::save);
	
	hSplitter->addWidget(view);
	hSplitter->addWidget(sidePanel);
	vSplitter->setOrientation(Qt::Vertical);
	vSplitter->addWidget(hSplitter);
	vSplitter->addWidget(bottomPanel);
	
	statusBar->setSizeGripEnabled(false);
	
	gridLayout->setContentsMargins(0, 0, 0, 0);
	
	gridLayout->addWidget(fileBar, 0, 0);
	gridLayout->addWidget(vSplitter, 1, 0);
	gridLayout->addWidget(statusBar, 2, 0);
	
	setLayout(gridLayout);
}

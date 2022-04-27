#ifndef GUI_H
#define GUI_H

#include <QWidget>
#include <QGridLayout>
#include <QLabel>
#include <QApplication>
#include <QGroupBox>
#include <QSplitter>
#include <QSaveFile>
#include <QFileDialog>
#include <QMessageBox>

#include "../GUI/include/BottomPanel.h"
#include "../GUI/include/FileBar.h"
#include "../GUI/include/ItemTree.h"
#include "../GUI/include/SidePanel.h"
#include "../GUI/include/StatusBar.h"
#include "../GUI/include/View.h"

class GUI : public QWidget {
Q_OBJECT
public:
	explicit GUI(QWidget *parent = nullptr);
	
	~GUI() override = default;

private:
	QGridLayout *gridLayout;
	QSplitter *hSplitter;
	QSplitter *vSplitter;
	
	BottomPanel *bottomPanel;
	FileBar *fileBar;
	SidePanel *sidePanel;
	StatusBar *statusBar;
	View *view;
public slots:
signals:
};

#endif //GUI_H

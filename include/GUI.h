#ifndef GUI_H
#define GUI_H

#include <QApplication>
#include <QFileDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QMessageBox>
#include <QSaveFile>
#include <QSplitter>
#include <QWidget>

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
};

#endif //GUI_H

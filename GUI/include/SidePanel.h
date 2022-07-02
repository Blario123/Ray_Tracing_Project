#ifndef SIDEPANEL_H
#define SIDEPANEL_H

#include <QGridLayout>
#include <QGroupBox>
#include <QLineEdit>
#include <QPushButton>
#include <QScreen>
#include <QWidget>

#include "ItemTree.h"

class SidePanel : public QWidget {
Q_OBJECT
public:
	explicit SidePanel(QWidget *parent = nullptr);
	~SidePanel() override = default;
private:
	ItemTree *itemTree;
	QGridLayout *layout;
	QGroupBox *groupBox;
	QGridLayout *boxLayout;
	QPushButton *addPushButton;
	QPushButton *delPushButton;
};

#endif //SIDEPANEL_H

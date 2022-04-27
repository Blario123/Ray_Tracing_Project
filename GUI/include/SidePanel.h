#ifndef SIDEPANEL_H
#define SIDEPANEL_H

#include <QWidget>
#include <QLineEdit>
#include <QGroupBox>
#include <QGridLayout>
#include <QPushButton>
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

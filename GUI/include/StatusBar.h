#ifndef STATUSBAR_H
#define STATUSBAR_H

#include <QStatusBar>

class StatusBar : public QStatusBar {
Q_OBJECT
public:
	explicit StatusBar(QWidget *parent = nullptr);
	~StatusBar() override = default;
};

#endif //STATUSBAR_H

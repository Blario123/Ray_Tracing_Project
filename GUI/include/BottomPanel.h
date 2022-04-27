#ifndef BOTTOMPANEL_H
#define BOTTOMPANEL_H

#include <QWidget>
#include <QBoxLayout>
#include <QGroupBox>

class BottomPanel : public QWidget {
Q_OBJECT
public:
	explicit BottomPanel(QWidget *parent = nullptr);
	
	~BottomPanel() override = default;

private:
	QBoxLayout *layout;
	QGroupBox *groupBox;
};

#endif //BOTTOMPANEL_H

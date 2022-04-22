#ifndef GUI_H
#define GUI_H

#include <QWidget>
#include <QGridLayout>
#include <QLabel>

class GUI : public QWidget {
Q_OBJECT
public:
    explicit GUI(QWidget *parent = nullptr);
    ~GUI() override = default;
};

#endif //GUI_H

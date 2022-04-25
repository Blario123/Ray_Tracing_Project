#ifndef GUI_H
#define GUI_H

#include <QWidget>
#include <QGridLayout>
#include <QLabel>

#include "GUI/BottomPanel.h"
#include "GUI/FileBar.h"
#include "GUI/SidePanel.h"
#include "GUI/StatusBar.h"
#include "GUI/View.h"

class GUI : public QWidget {
Q_OBJECT
public:
    explicit GUI(QWidget *parent = nullptr);
    ~GUI() override = default;
};

#endif //GUI_H

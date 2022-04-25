//
// Created by blario123 on 23/04/2022.
//

#ifndef SIDEPANEL_H
#define SIDEPANEL_H

#include <QWidget>

class SidePanel : public QWidget {
Q_OBJECT
public:
    explicit SidePanel(QWidget *parent = nullptr);
    ~SidePanel() override = default;
};

#endif //SIDEPANEL_H

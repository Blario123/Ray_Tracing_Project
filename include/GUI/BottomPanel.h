//
// Created by blario123 on 23/04/2022.
//

#ifndef BOTTOMPANEL_H
#define BOTTOMPANEL_H

#include <QWidget>

class BottomPanel : public QWidget {
Q_OBJECT
public:
    explicit BottomPanel(QWidget *parent = nullptr);
    ~BottomPanel() override = default;
};

#endif //BOTTOMPANEL_H

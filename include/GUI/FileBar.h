//
// Created by blario123 on 23/04/2022.
//

#ifndef FILEBAR_H
#define FILEBAR_H

#include <QMenuBar>

class FileBar : public QMenuBar {
Q_OBJECT
public:
    explicit FileBar(QWidget *parent = nullptr);
    ~FileBar() override = default;
};

#endif //FILEBAR_H

//
// Created by blario123 on 23/04/2022.
//

#ifndef VIEW_H
#define VIEW_H

#include <QWidget>

class View : public QWidget {
Q_OBJECT
public:
    explicit View(QWidget *parent = nullptr);
    ~View() override = default;
};

#endif //VIEW_H

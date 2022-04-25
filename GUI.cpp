#include "GUI.h"

GUI::GUI(QWidget *parent) : QWidget(parent) {
    auto *gl = new QGridLayout;

    gl->addWidget(l);
    setLayout(gl);
}



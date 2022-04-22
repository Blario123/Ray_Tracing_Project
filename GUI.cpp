#include "GUI.h"

GUI::GUI(QWidget *parent) : QWidget(parent) {
    auto *gl = new QGridLayout;
    auto *l = new QLabel("AAA");
    gl->addWidget(l);
    setLayout(gl);
}



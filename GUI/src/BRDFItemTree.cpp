#include "BRDFItemTree.h"

BRDFItemTree::BRDFItemTree(QWidget *parent) : QTreeWidget(parent) {
    setHeaderLabels(QStringList() << "Name" << "Red" << "Green" << "Blue");
}
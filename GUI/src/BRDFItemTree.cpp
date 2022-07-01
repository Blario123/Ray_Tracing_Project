#include "BRDFItemTree.h"

BRDFItemTree::BRDFItemTree(QWidget *parent) : QTreeWidget(parent) {
    setHeaderLabels(QStringList() << "Name" << "Red" << "Green" << "Blue");
}

void BRDFItemTree::setList(QList<BRDF> &list) {
    QList<BRDF>::iterator i;
    for(i = list.begin(); i < list.end(); i++) {
        QTreeWidgetItem *tempItem = new QTreeWidgetItem;
        tempItem->setText(0, i->name());
        tempItem->setData(1, 0, i->r());
        tempItem->setData(2, 0, i->g());
        tempItem->setData(3, 0, i->b());
        this->addTopLevelItem(tempItem);
    }
}

void BRDFItemTree::addItem(const QTreeWidgetItem &i) {
    addTopLevelItem(i.clone());
}

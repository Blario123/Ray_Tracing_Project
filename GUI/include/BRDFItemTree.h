#ifndef BRDFITEMTREE_H
#define BRDFITEMTREE_H

#include <QListWidgetItem>
#include <QTreeWidget>

#include "BRDF.h"

class BRDFItemTree : public QTreeWidget {
Q_OBJECT
public:
    explicit BRDFItemTree(QWidget *parent = nullptr);
    ~BRDFItemTree() override = default;
    void setList(QList<BRDF> &);
private:
    QList<BRDF> pList;
public slots:
    void addItem(QTreeWidgetItem *i);
    void delItem();
signals:
    void deleteAt(int);
};

#endif //BRDFITEMTREE_H

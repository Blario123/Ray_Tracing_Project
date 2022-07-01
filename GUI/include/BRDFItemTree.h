#ifndef BRDFITEMTREE_H
#define BRDFITEMTREE_H

#include "BRDF.h"
#include <QTreeWidget>

class BRDFItemTree : public QTreeWidget {
Q_OBJECT
public:
    explicit BRDFItemTree(QWidget *parent = nullptr);
    ~BRDFItemTree() override = default;
    void setList(QList<BRDF> &);
private:
    QList<BRDF> pList;
public slots:
    void addItem(const QTreeWidgetItem &);
signals:
};

#endif //BRDFITEMTREE_H

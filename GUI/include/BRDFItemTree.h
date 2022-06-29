#ifndef BRDFITEMTREE_H
#define BRDFITEMTREE_H

#include <QTreeWidget>

class BRDFItemTree : public QTreeWidget {
Q_OBJECT
public:
    explicit BRDFItemTree(QWidget *parent = nullptr);
    ~BRDFItemTree() override = default;
public slots:
signals:
};

#endif //BRDFITEMTREE_H

#ifndef ITEMTREE_H
#define ITEMTREE_H

#include <QComboBox>
#include <QLabel>
#include <QTreeWidget>
#include <QWidget>

class ItemTree : public QTreeWidget {
Q_OBJECT
public:
	explicit ItemTree(QWidget *parent = nullptr);
	~ItemTree() override = default;
public slots:
	void addPressed();
	void delPressed();
	QTreeWidgetItem *createItem(const QString &name);
};

#endif //ITEMTREE_H

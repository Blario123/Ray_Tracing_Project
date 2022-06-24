#include <QLabel>
#include "ItemTree.h"

ItemTree::ItemTree(QWidget *parent) : QTreeWidget(parent) {
	setHeaderLabels({"Item", "Shape", "Name", "x", "y", "z"});

}

void ItemTree::addPressed() {
	addTopLevelItem(createItem("AAA"));
}

void ItemTree::delPressed() {
	takeTopLevelItem(currentIndex().row());
}

QTreeWidgetItem *ItemTree::createItem(const QString &name) {
	auto *item = new QTreeWidgetItem(this);
	auto *cb = new QComboBox;
	cb->addItems(QStringList() << "Sphere" << "Cube" << "Torus" << "Cuboid" << "Reuleaux");
	setItemWidget(item, 1, cb);
    setItemWidget(item, 2, new QLabel(name));
	item->setText(0, QString::number(topLevelItemCount()));
//	item->setText(2, )
//	connect(this, ItemTree::currentItemChanged())
	item->setFlags(item->flags() | Qt::ItemIsEditable);
	return item;
}

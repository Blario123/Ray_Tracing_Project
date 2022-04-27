#include "ItemTree.h"

ItemTree::ItemTree(QWidget *parent) : QTreeWidget(parent) {
	setHeaderLabels({"Item", "Shape", "Name", "x", "y", "z"});
//	header().
// TODO: Setup tree view to have shape, and list and to be able to add items.
}

void ItemTree::addPressed() {
	addTopLevelItem(createItem("AAA"));
}

void ItemTree::delPressed() {
	takeTopLevelItem(currentIndex().row());
}

QTreeWidgetItem *ItemTree::createItem(const QString& name) {
	auto *item = new QTreeWidgetItem(this);
	auto *cb = new QComboBox;
	cb->addItems(QStringList() << "Sphere" << "Cube");
	setItemWidget(item, 1, cb);
	item->setText(0, QString::number(topLevelItemCount()));
//	item->setText(2, )
//	connect(this, ItemTree::currentItemChanged())
	item->setFlags(item->flags() | Qt::ItemIsEditable);
	return item;
}

#ifndef BRDFITEMDIALOG_H
#define BRDFITEMDIALOG_H

#include <QSpinBox>
#include <QLineEdit>
#include <QLabel>
#include <QTreeWidgetItem>
#include <QDialog>
#include <QObject>

class BRDFItemDialog : public QDialog {
    Q_OBJECT
public:
    explicit BRDFItemDialog(QWidget *parent = nullptr) :
        QDialog(parent),
        nameLabel(new QLabel("Name")),
        nameLineEdit(new QLineEdit),
        rLabel(new QLabel("Red")),
        rSpinBox(new QDoubleSpinBox),
        gLabel(new QLabel("Green")),
        gSpinBox(new QDoubleSpinBox),
        bLabel(new QLabel("Blue")),
        bSpinBox(new QDoubleSpinBox) {

    }
    ~BRDFItemDialog() override = default;

private:
    QLabel *nameLabel;
    QLineEdit *nameLineEdit;
    QLabel *rLabel;
    QDoubleSpinBox *rSpinBox;
    QLabel *gLabel;
    QDoubleSpinBox *gSpinBox;
    QLabel *bLabel;
    QDoubleSpinBox *bSpinBox;
    void generateItem() {
        QTreeWidgetItem *item = new QTreeWidgetItem;
        item->setText(0, nameLineEdit->text());
        item->setData(1, 0, rSpinBox->value());
        item->setData(2, 0, gSpinBox->value());
        item->setData(3, 0, bSpinBox->value());
        emit itemAdded(*item);
    }
public slots:
    void accept() override {
        QDialog::accept();
    }
signals:
    void itemAdded(const QTreeWidgetItem &);
};

#endif // BRDFITEMDIALOG_H

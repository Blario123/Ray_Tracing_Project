#ifndef BRDFITEMDIALOG_H
#define BRDFITEMDIALOG_H

#include <QDialog>
#include <QDialogButtonBox>
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QObject>
#include <QPushButton>
#include <QSpinBox>
#include <QTreeWidgetItem>

class BRDFItemDialog : public QDialog {
    Q_OBJECT
public:
    explicit BRDFItemDialog(QWidget *parent = nullptr) :QDialog(parent),
                                                        layout(new QGridLayout),
                                                        nameLabel(new QLabel("Name")),
                                                        nameLineEdit(new QLineEdit),
                                                        rLabel(new QLabel("Red")),
                                                        rSpinBox(new QDoubleSpinBox),
                                                        gLabel(new QLabel("Green")),
                                                        gSpinBox(new QDoubleSpinBox),
                                                        bLabel(new QLabel("Blue")),
                                                        bSpinBox(new QDoubleSpinBox) {
        setWindowTitle("Add BRDF");
        setModal(true);
        auto *buttonBox = new QDialogButtonBox(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);

        layout->addWidget(nameLabel, 0, 0);
        layout->addWidget(nameLineEdit, 0, 1);
        layout->addWidget(rLabel, 1, 0);
        layout->addWidget(rSpinBox, 1, 1);
        layout->addWidget(gLabel, 2, 0);
        layout->addWidget(gSpinBox, 2, 1);
        layout->addWidget(bLabel, 3, 0);
        layout->addWidget(bSpinBox, 3, 1);
        layout->addWidget(buttonBox, 4, 0, 1, 2);

        setLayout(layout);
        connect(buttonBox->button(QDialogButtonBox::Ok), &QPushButton::pressed, this, &BRDFItemDialog::accept);
        connect(buttonBox->button(QDialogButtonBox::Cancel), &QPushButton::pressed, this, &BRDFItemDialog::close);
    }
    ~BRDFItemDialog() override = default;
private:
    QGridLayout *layout;
    QLabel *nameLabel;
    QLineEdit *nameLineEdit;
    QLabel *rLabel;
    QDoubleSpinBox *rSpinBox;
    QLabel *gLabel;
    QDoubleSpinBox *gSpinBox;
    QLabel *bLabel;
    QDoubleSpinBox *bSpinBox;
    bool generateItem() {
        auto *item = new QTreeWidgetItem;
        item->setText(0, nameLineEdit->text());
        item->setData(1, 0, rSpinBox->value());
        item->setData(2, 0, gSpinBox->value());
        item->setData(3, 0, bSpinBox->value());
        if(item->data(0, 0).toString().isEmpty()) {
            return false;
        }
        emit itemAdded(item);
        delete item;
        return true;
    }
public slots:
    void accept() override {
        if(generateItem()) {
            QDialog::accept();
        } else {
            auto *mb = new QMessageBox(QMessageBox::Warning, "Error", "Name cannot be empty.", QMessageBox::Ok, this);
            mb->show();
        }
    }
signals:
    void itemAdded(QTreeWidgetItem *);
};

#endif // BRDFITEMDIALOG_H

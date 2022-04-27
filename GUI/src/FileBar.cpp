#include "FileBar.h"

FileBar::FileBar(QWidget *parent) : QMenuBar(parent),
									saveAction(new QAction("Save")),
									closeAction(new QAction("Close")),
									aboutAction(new QAction("About")) {
	fileMenu = this->addMenu("&File");
	fileMenu->addAction(saveAction);
	fileMenu->addSeparator();
	fileMenu->addAction(closeAction);
	aboutMenu = this->addMenu("&About");
	aboutMenu->addAction(aboutAction);
}

void FileBar::save() {
	QString fileName = QFileDialog::getSaveFileName(this, "Save file.", "", "Portable Network Graphic (*.png);;All Files (*)");
	if(fileName.isEmpty()) {
		return;
	} else {
		QFile file(fileName);
		if(!file.open(QIODevice::WriteOnly)) {
			QMessageBox::information(this, "Unable to open file.", file.errorString());
			return;
		}
		QDataStream out(&file);
		out.setVersion(QDataStream::Qt_6_2);
		// TODO: Find way to data-stream final product through FileBar? Or pull data-stream from a parent?
	}
}

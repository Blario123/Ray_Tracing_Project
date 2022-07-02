#ifndef VIEW_H
#define VIEW_H

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGridLayout>
#include <QImage>
#include <QResizeEvent>
#include <QWidget>

class View : public QWidget {
Q_OBJECT
public:
	explicit View(QWidget *parent = nullptr);
	~View() override = default;
    void resizeEvent(QResizeEvent *) override;
private:
    QGridLayout *layout;
    QGraphicsView *view;
    QGraphicsScene *scene;
    QImage *image;
public slots:
    void renderImage();
signals:
    void sizeChanged(int, int);
};

#endif //VIEW_H

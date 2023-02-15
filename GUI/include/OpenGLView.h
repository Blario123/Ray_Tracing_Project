#ifndef OPENGLVIEW_H
#define OPENGLVIEW_H

#include <QOpenGLWidget>

class OpenGLView : public QOpenGLWidget {
Q_OBJECT
public:
    explicit OpenGLView(QOpenGLWidget *parent = nullptr);
    ~OpenGLView() override = default;
};


#endif // OPENGLVIEW_H
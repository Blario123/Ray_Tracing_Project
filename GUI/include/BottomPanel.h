#ifndef BOTTOMPANEL_H
#define BOTTOMPANEL_H

#include <QWidget>
#include <QBoxLayout>
#include <QGroupBox>
#include <QTabWidget>
#include <QSpinBox>
#include <QLabel>
#include <QCheckBox>
#include <QComboBox>
#include <QPushButton>
#include <QTreeWidgetItem>

#include "BRDF.h"
#include "BRDFItemTree.h"
#include "BRDFItemDialog.h"

class BottomPanel : public QWidget {
Q_OBJECT
public:
	explicit BottomPanel(QWidget *parent = nullptr);
	~BottomPanel() override = default;
private:
    QGridLayout *layout;
    QGridLayout *groupBoxLayout;
	QGroupBox *groupBox;
    QTabWidget *tabWidget;
    QWidget *cameraTab;
    QGridLayout *cameraLayout;
    QWidget *renderTab;
    QGridLayout *renderLayout;
    QWidget *structureTab;
    QGridLayout *structureLayout;
    QWidget *brdfTab;
    QGridLayout *brdfLayout;
    int viewWidth, viewHeight;
    // Camera Tab
    QLabel *cameraCameraPositionLabel;
    QDoubleSpinBox *cameraCameraXPositionSpinBox;
    QDoubleSpinBox *cameraCameraYPositionSpinBox;
    QDoubleSpinBox *cameraCameraZPositionSpinBox;
    QLabel *cameraCameraDirectionLabel;
    QDoubleSpinBox *cameraCameraXDirectionSpinBox;
    QDoubleSpinBox *cameraCameraYDirectionSpinBox;
    QDoubleSpinBox *cameraCameraZDirectionSpinBox;
    QLabel *cameraUpwardDirectionLabel;
    QDoubleSpinBox *cameraUpwardXDirectionSpinBox;
    QDoubleSpinBox *cameraUpwardYDirectionSpinBox;
    QDoubleSpinBox *cameraUpwardZDirectionSpinBox;
    QLabel *cameraHzFOVAngleLabel;
    QDoubleSpinBox *cameraHzFOVAngleSpinBox;
    QLabel *cameraFocalDistanceLabel;
    QDoubleSpinBox *cameraFocalDistanceSpinBox;
    QLabel *cameraLensRadiusLabel;
    QDoubleSpinBox *cameraLensRadiusSpinBox;
    QLabel *cameraResolutionLabel;
    QSpinBox *cameraXResolutionSpinBox;
    QSpinBox *cameraYResolutionSpinBox;
    QLabel *cameraFollowResolutionLabel;
    QCheckBox *cameraFollowResolutionCheckBox;
    // Render Tab
    QLabel *renderNBouncesLabel;
    QSpinBox *renderNBouncesSpinBox;
    QLabel *renderNRandomSamplesLabel;
    QSpinBox *renderNRandomSamplesSpinBox;
    QLabel *renderThresholdIntersectionDistanceLabel;
    QDoubleSpinBox *renderThresholdIntersectionDistanceSpinBox;
    QLabel *renderThresholdNoIntersectionDistanceLabel;
    QDoubleSpinBox *renderThresholdNoIntersectionDistanceSpinBox;
    QLabel *renderFiniteDifferenceSizeLabel;
    QDoubleSpinBox *renderFiniteDifferenceSizeSpinBox;
    QLabel *renderSupersamplingValueLabel;
    QSpinBox *renderSupersamplingValueSpinBox;
    QLabel *renderIsSilentLabel;
    QCheckBox *renderIsSilentCheckBox;
    QLabel *renderRenderModeLabel;
    QComboBox *renderRenderModeComboBox;
    QPushButton *renderRenderPushButton;
    // Structure Tab
    QLabel *structureFloorBRDFLabel;
    QComboBox *structureFloorBRDFComboBox;
    QLabel *structureCeilingBRDFLabel;
    QComboBox *structureCeilingBRDFComboBox;
    QLabel *structureLeftBRDFLabel;
    QComboBox *structureLeftBRDFComboBox;
    QLabel *structureRightBRDFLabel;
    QComboBox *structureRightBRDFComboBox;
    QLabel *structureBackBRDFLabel;
    QComboBox *structureBackBRDFComboBox;
    QLabel *structureSphere1BRDFLabel;
    QComboBox *structureSphere1BRDFComboBox;
    QLabel *structureSphere2BRDFLabel;
    QComboBox *structureSphere2BRDFComboBox;
    // BRDF Tab
    BRDFItemTree *brdfTreeWidget;
    QPushButton *brdfAddBRDFPushButton;
    QPushButton *brdfDelBRDFPushButton;
    bool followResolution = false;
    BRDFItemDialog *brdfItemDialog;
public slots:
    void onRenderPushButtonPressed();
    void onViewSizeChanged(int,int);
    void toggleFollowResolution(bool);
    void onAddPushButtonPressed();
};

#endif //BOTTOMPANEL_H

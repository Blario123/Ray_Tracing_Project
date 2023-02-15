#ifndef BOTTOMPANEL_H
#define BOTTOMPANEL_H

#include <QBoxLayout>
#include <QCheckBox>
#include <QComboBox>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QScreen>
#include <QSpinBox>
#include <QTabWidget>
#include <QTreeWidgetItem>
#include <QWidget>

#include "BRDF.h"
#include "BRDFItemDialog.h"
#include "BRDFItemTree.h"

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
    QList<BRDF> brdfList;
    QStringList brdfNameList;
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
    // BRDF Tab
    BRDFItemTree *brdfTreeWidget;
    QPushButton *brdfAddBRDFPushButton;
    QPushButton *brdfDelBRDFPushButton;
    bool followResolution = false;
    BRDFItemDialog *brdfItemDialog;
public slots:
    void updateList();
    void onRenderPushButtonPressed();
    void onViewSizeChanged(int, int);
    void toggleFollowResolution(bool);
    void onAddPushButtonPressed();
    void addItem(QTreeWidgetItem *);
    void delItem(int);
};

#endif //BOTTOMPANEL_H

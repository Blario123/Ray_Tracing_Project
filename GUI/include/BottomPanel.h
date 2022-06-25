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
    // Structure Tab
    QStringList brdfList;
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
};

#endif //BOTTOMPANEL_H

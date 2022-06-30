#include "BottomPanel.h"

BottomPanel::BottomPanel(QWidget *parent) : QWidget(parent),
											layout(new QGridLayout),
											groupBoxLayout(new QGridLayout),
											groupBox(new QGroupBox("BottomPanel")),
                                            tabWidget(new QTabWidget),
                                            cameraTab(new QWidget),
                                            cameraLayout(new QGridLayout),
                                            renderTab(new QWidget),
                                            renderLayout(new QGridLayout),
                                            structureTab(new QWidget),
                                            structureLayout(new QGridLayout),
                                            brdfTab(new QWidget),
                                            brdfLayout(new QGridLayout),
                                            viewWidth(0),
                                            viewHeight(0),
                                            cameraCameraPositionLabel(new QLabel("Camera position")),
                                            cameraCameraXPositionSpinBox(new QDoubleSpinBox),
                                            cameraCameraYPositionSpinBox(new QDoubleSpinBox),
                                            cameraCameraZPositionSpinBox(new QDoubleSpinBox),
                                            cameraCameraDirectionLabel(new QLabel("Camera direction")),
                                            cameraCameraXDirectionSpinBox(new QDoubleSpinBox),
                                            cameraCameraYDirectionSpinBox(new QDoubleSpinBox),
                                            cameraCameraZDirectionSpinBox(new QDoubleSpinBox),
                                            cameraUpwardDirectionLabel(new QLabel("Upward direction")),
                                            cameraUpwardXDirectionSpinBox(new QDoubleSpinBox),
                                            cameraUpwardYDirectionSpinBox(new QDoubleSpinBox),
                                            cameraUpwardZDirectionSpinBox(new QDoubleSpinBox),
                                            cameraHzFOVAngleLabel(new QLabel("Horizontal FOV angle")),
                                            cameraHzFOVAngleSpinBox(new QDoubleSpinBox),
                                            cameraFocalDistanceLabel(new QLabel("Focal distance")),
                                            cameraFocalDistanceSpinBox(new QDoubleSpinBox),
                                            cameraLensRadiusLabel(new QLabel("Lens radius")),
                                            cameraLensRadiusSpinBox(new QDoubleSpinBox),
                                            cameraResolutionLabel(new QLabel("Camera resolution")),
                                            cameraXResolutionSpinBox(new QSpinBox),
                                            cameraYResolutionSpinBox(new QSpinBox),
                                            cameraFollowResolutionLabel(new QLabel("Follow resolution?")),
                                            cameraFollowResolutionCheckBox(new QCheckBox),
                                            renderNBouncesLabel(new QLabel("Number of bounces")),
                                            renderNBouncesSpinBox(new QSpinBox),
                                            renderNRandomSamplesLabel(new QLabel("Number of random samples")),
                                            renderNRandomSamplesSpinBox(new QSpinBox),
                                            renderThresholdIntersectionDistanceLabel(new QLabel("Threshold intersection distance")),
                                            renderThresholdIntersectionDistanceSpinBox(new QDoubleSpinBox),
                                            renderThresholdNoIntersectionDistanceLabel(new QLabel("Threshold no intersection distance")),
                                            renderThresholdNoIntersectionDistanceSpinBox(new QDoubleSpinBox),
                                            renderFiniteDifferenceSizeLabel(new QLabel("Finite difference size")),
                                            renderFiniteDifferenceSizeSpinBox(new QDoubleSpinBox),
                                            renderSupersamplingValueLabel(new QLabel("Supersampling value")),
                                            renderSupersamplingValueSpinBox(new QSpinBox),
                                            renderIsSilentLabel(new QLabel("Silent render?")),
                                            renderIsSilentCheckBox(new QCheckBox),
                                            renderRenderModeLabel(new QLabel("Render mode")),
                                            renderRenderModeComboBox(new QComboBox),
                                            renderRenderPushButton(new QPushButton("Render")),
                                            structureFloorBRDFLabel(new QLabel("Floor")),
                                            structureFloorBRDFComboBox(new QComboBox),
                                            structureCeilingBRDFLabel(new QLabel("Ceiling")),
                                            structureCeilingBRDFComboBox(new QComboBox),
                                            structureLeftBRDFLabel(new QLabel("Left")),
                                            structureLeftBRDFComboBox(new QComboBox),
                                            structureRightBRDFLabel(new QLabel("Right")),
                                            structureRightBRDFComboBox(new QComboBox),
                                            structureBackBRDFLabel(new QLabel("Back")),
                                            structureBackBRDFComboBox(new QComboBox),
                                            brdfTreeWidget(new BRDFItemTree),
                                            brdfAddBRDFPushButton(new QPushButton("Add")),
                                            brdfDelBRDFPushButton(new QPushButton("Del")) {
    cameraCameraXPositionSpinBox->setValue(0.0);
    cameraCameraYPositionSpinBox->setValue(0.0);
    cameraCameraZPositionSpinBox->setValue(1.0);
    cameraCameraXDirectionSpinBox->setValue(1.0);
    cameraCameraYDirectionSpinBox->setValue(0.0);
    cameraCameraZDirectionSpinBox->setValue(0.0);
    cameraUpwardXDirectionSpinBox->setValue(0.0);
    cameraUpwardYDirectionSpinBox->setValue(0.0);
    cameraUpwardZDirectionSpinBox->setValue(1.0);
    cameraHzFOVAngleSpinBox->setValue(pi / 2.5);
    cameraFocalDistanceSpinBox->setValue(2.0);
    cameraLensRadiusSpinBox->setValue(0.0);
    cameraXResolutionSpinBox->setRange(1, 10000);
    cameraXResolutionSpinBox->setValue(800);
    cameraYResolutionSpinBox->setRange(1, 10000);
    cameraYResolutionSpinBox->setValue(400);

    brdfList << "Red" << "Blue" << "White" << "Pink" << "Purple" << "Yellow" << "Brown" << "Cyan";
    brdfTreeWidget->setFixedHeight(100);

    renderRenderModeComboBox->addItems(QStringList() << "RenderImage" << "RenderImagePerThread" << "RenderImageMultithreaded" << "RenderImageRussian" << "RenderImageRussianPerThread" << "RenderImageRussianMultithreaded");

    structureFloorBRDFComboBox->addItems(brdfList);
    structureFloorBRDFComboBox->setCurrentIndex(2);
    structureCeilingBRDFComboBox->addItems(brdfList);
    structureCeilingBRDFComboBox->setCurrentIndex(2);
    structureLeftBRDFComboBox->addItems(brdfList);
    structureLeftBRDFComboBox->setCurrentIndex(0);
    structureRightBRDFComboBox->addItems(brdfList);
    structureRightBRDFComboBox->setCurrentIndex(1);
    structureBackBRDFComboBox->addItems(brdfList);
    structureBackBRDFComboBox->setCurrentIndex(2);

    renderNBouncesSpinBox->setValue(10);
    renderNRandomSamplesSpinBox->setValue(10);
    renderThresholdIntersectionDistanceSpinBox->setRange(0, 1);
    renderThresholdIntersectionDistanceSpinBox->setDecimals(9);
    renderThresholdIntersectionDistanceSpinBox->setValue(1.0e-8);
    renderThresholdNoIntersectionDistanceSpinBox->setRange(0, 200);
    renderThresholdNoIntersectionDistanceSpinBox->setDecimals(3);
    renderThresholdNoIntersectionDistanceSpinBox->setValue(100.0);
    renderFiniteDifferenceSizeSpinBox->setRange(0, 1);
    renderFiniteDifferenceSizeSpinBox->setDecimals(9);
    renderFiniteDifferenceSizeSpinBox->setValue(1.0e-8);

    cameraLayout->addWidget(cameraCameraPositionLabel, 0, 0);
    cameraLayout->addWidget(cameraCameraXPositionSpinBox, 1, 0);
    cameraLayout->addWidget(cameraCameraYPositionSpinBox, 2, 0);
    cameraLayout->addWidget(cameraCameraZPositionSpinBox, 3, 0);
    cameraLayout->addWidget(cameraCameraDirectionLabel, 0, 1);
    cameraLayout->addWidget(cameraCameraXDirectionSpinBox, 1, 1);
    cameraLayout->addWidget(cameraCameraYDirectionSpinBox, 2, 1);
    cameraLayout->addWidget(cameraCameraZDirectionSpinBox, 3, 1);
    cameraLayout->addWidget(cameraUpwardDirectionLabel, 0, 2);
    cameraLayout->addWidget(cameraUpwardXDirectionSpinBox, 1, 2);
    cameraLayout->addWidget(cameraUpwardYDirectionSpinBox, 2, 2);
    cameraLayout->addWidget(cameraUpwardZDirectionSpinBox, 3, 2);
    cameraLayout->addWidget(cameraHzFOVAngleLabel, 0, 3);
    cameraLayout->addWidget(cameraHzFOVAngleSpinBox, 1, 3);
    cameraLayout->addWidget(cameraFocalDistanceLabel, 0, 4);
    cameraLayout->addWidget(cameraFocalDistanceSpinBox, 1, 4);
    cameraLayout->addWidget(cameraLensRadiusLabel, 2, 4);
    cameraLayout->addWidget(cameraLensRadiusSpinBox, 3, 4);
    cameraLayout->addWidget(cameraResolutionLabel, 0, 5, 1, 2);
    cameraLayout->addWidget(cameraXResolutionSpinBox, 1, 5, 1, 2);
    cameraLayout->addWidget(cameraYResolutionSpinBox, 2, 5, 1, 2);
    cameraLayout->addWidget(cameraFollowResolutionLabel, 3, 5);
    cameraLayout->addWidget(cameraFollowResolutionCheckBox, 3, 6);

    renderLayout->addWidget(renderNBouncesLabel, 0, 0);
    renderLayout->addWidget(renderNBouncesSpinBox, 1, 0);
    renderLayout->addWidget(renderNRandomSamplesLabel, 2, 0);
    renderLayout->addWidget(renderNRandomSamplesSpinBox, 3, 0);
    renderLayout->addWidget(renderThresholdIntersectionDistanceLabel, 0, 1);
    renderLayout->addWidget(renderThresholdIntersectionDistanceSpinBox, 1, 1);
    renderLayout->addWidget(renderThresholdNoIntersectionDistanceLabel, 2, 1);
    renderLayout->addWidget(renderThresholdNoIntersectionDistanceSpinBox, 3, 1);
    renderLayout->addWidget(renderFiniteDifferenceSizeLabel, 0, 2);
    renderLayout->addWidget(renderFiniteDifferenceSizeSpinBox, 1, 2);
    renderLayout->addWidget(renderSupersamplingValueLabel, 2, 2);
    renderLayout->addWidget(renderSupersamplingValueSpinBox, 3, 2);
    renderLayout->addWidget(renderRenderModeLabel, 0, 3, 1, 2);
    renderLayout->addWidget(renderRenderModeComboBox, 1, 3, 1, 2);
    renderLayout->addWidget(renderIsSilentLabel, 2, 3, 1, 2);
    renderLayout->addWidget(renderIsSilentCheckBox, 3, 3);
    renderLayout->addWidget(renderRenderPushButton, 3, 4);

    brdfLayout->addWidget(brdfTreeWidget, 0, 0, 1, 2);
    brdfLayout->addWidget(brdfAddBRDFPushButton, 1, 0);
    brdfLayout->addWidget(brdfDelBRDFPushButton, 1, 1);

    structureLayout->addWidget(structureFloorBRDFLabel, 0, 0);
    structureLayout->addWidget(structureFloorBRDFComboBox, 1, 0);
    structureLayout->addWidget(structureCeilingBRDFLabel, 2, 0);
    structureLayout->addWidget(structureCeilingBRDFComboBox, 3, 0);
    structureLayout->addWidget(structureLeftBRDFLabel, 0, 1);
    structureLayout->addWidget(structureLeftBRDFComboBox, 1, 1);
    structureLayout->addWidget(structureRightBRDFLabel, 2, 1);
    structureLayout->addWidget(structureRightBRDFComboBox, 3, 1);
    structureLayout->addWidget(structureBackBRDFLabel, 0, 2);
    structureLayout->addWidget(structureBackBRDFComboBox, 1, 2);

    cameraTab->setLayout(cameraLayout);
    renderTab->setLayout(renderLayout);
    brdfTab->setLayout(brdfLayout);
    structureTab->setLayout(structureLayout);

    tabWidget->addTab(cameraTab, "Camera");
    tabWidget->addTab(structureTab, "Structure");
    tabWidget->addTab(brdfTab, "BRDF");
    tabWidget->addTab(renderTab, "Render");

	groupBoxLayout->addWidget(tabWidget);
    groupBox->setLayout(groupBoxLayout);

    layout->setContentsMargins(0, 0, 0, 0);

    layout->addWidget(groupBox);
	setLayout(layout);

    setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    connect(renderRenderPushButton, &QPushButton::pressed, this, &BottomPanel::onRenderPushButtonPressed);
    connect(cameraFollowResolutionCheckBox, &QCheckBox::toggled, this, &BottomPanel::toggleFollowResolution);
}

void BottomPanel::onRenderPushButtonPressed() {

}

void BottomPanel::onViewSizeChanged(int w, int h) {
    if(followResolution) {
        cameraXResolutionSpinBox->setValue(w);
        cameraYResolutionSpinBox->setValue(h);
    }
}

void BottomPanel::toggleFollowResolution(bool state) {
    followResolution = state;
}

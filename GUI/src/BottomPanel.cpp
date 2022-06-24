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
                                            renderRenderModeComboBox(new QComboBox) {
    renderRenderModeComboBox->addItems(QStringList() << "RenderImage" << "RenderImagePerThread" << "RenderImageMultithreaded" << "RenderImageRussian" << "RenderImageRussianPerThread" << "RenderImageRussianMultithreaded");

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
    cameraLayout->addWidget(cameraResolutionLabel, 0, 5);
    cameraLayout->addWidget(cameraXResolutionSpinBox, 1, 5);
    cameraLayout->addWidget(cameraYResolutionSpinBox, 2, 5);

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
    renderLayout->addWidget(renderRenderModeLabel, 0, 3);
    renderLayout->addWidget(renderRenderModeComboBox, 1, 3);
    renderLayout->addWidget(renderIsSilentLabel, 2, 3);
    renderLayout->addWidget(renderIsSilentCheckBox, 3, 3);

    cameraTab->setLayout(cameraLayout);
    renderTab->setLayout(renderLayout);

    tabWidget->addTab(cameraTab, "Camera");
    tabWidget->addTab(renderTab, "Render");

	groupBoxLayout->addWidget(tabWidget);
    groupBox->setLayout(groupBoxLayout);

    layout->addWidget(groupBox);
	setLayout(layout);
}

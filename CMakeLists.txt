CMAKE_MINIMUM_REQUIRED(VERSION 3.20)
project(Ray_Tracing_Project)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

if (WIN32)
    set(CMAKE_PREFIX_PATH "C:/Qt/6.2.2/mingw_64/lib/cmake")
else ()
    set(CMAKE_PREFIX_PATH ~/Qt/6.3.0/gcc_64)
endif (WIN32)

find_package(Qt6 COMPONENTS
        Core
        Gui
        Widgets
        REQUIRED)

include_directories(include GUI/include)

set(CPPSOURCE
        main.cpp
        GUI.cpp
        GUI/src/BottomPanel.cpp
        GUI/src/FileBar.cpp
        GUI/src/ItemTree.cpp
        GUI/src/SidePanel.cpp
        GUI/src/StatusBar.cpp
        GUI/src/BRDFItemTree.cpp
        GUI/src/View.cpp)

set(HSOURCE
        include/main.h
        include/GUI.h
        include/Image.h
        include/Observer.h
        include/PhysicalObject.h
        include/Ray.h
        include/RayTracingProject.h
        include/SceneRender.h
        include/SDF.h
        include/stb_image_write.h
        include/Vec3.h
        include/GUI.h
        include/BRDF.h
        GUI/include/BottomPanel.h
        GUI/include/BRDFItemDialog.h
        GUI/include/FileBar.h
        GUI/include/SidePanel.h
        GUI/include/StatusBar.h
        GUI/include/BRDFItemTree.h
        GUI/include/View.h GUI/include/ItemTree.h)

add_executable(${PROJECT_NAME} ${CPPSOURCE} ${HSOURCE})

target_link_libraries(${PROJECT_NAME} Qt::Core Qt::Gui Qt::Widgets)

if (WIN32)
    set(DEBUG_SUFFIX)
    #    if (CMAKE_BUILD_TYPE MATCHES "Debug")
    #        set(DEBUG_SUFFIX "d")
    #    endif ()
    set(QT_INSTALL_PATH "${CMAKE_PREFIX_PATH}")
    if (NOT EXISTS "${QT_INSTALL_PATH}/bin")
        set(QT_INSTALL_PATH "${QT_INSTALL_PATH}/..")
        if (NOT EXISTS "${QT_INSTALL_PATH}/bin")
            set(QT_INSTALL_PATH "${QT_INSTALL_PATH}/..")
        endif ()
    endif ()
    if (EXISTS "${QT_INSTALL_PATH}/plugins/platforms/qwindows${DEBUG_SUFFIX}.dll")
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory
                "$<TARGET_FILE_DIR:${PROJECT_NAME}>/plugins/platforms/")
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                "${QT_INSTALL_PATH}/plugins/platforms/qwindows${DEBUG_SUFFIX}.dll"
                "$<TARGET_FILE_DIR:${PROJECT_NAME}>/plugins/platforms/")
    endif ()
    foreach (QT_LIB Core Gui Widgets)
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                "${QT_INSTALL_PATH}/bin/Qt6${QT_LIB}${DEBUG_SUFFIX}.dll"
                "$<TARGET_FILE_DIR:${PROJECT_NAME}>")
    endforeach (QT_LIB)
endif ()

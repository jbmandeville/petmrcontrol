#-------------------------------------------------
#
# Project created by QtCreator 2021-10-21T10:52:46
#
#-------------------------------------------------

QT       += core gui concurrent

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

LIBS += -L"$$_PRO_FILE_PWD_/3rdparty/fftw-3.3.8/lib" -lfftw3
INCLUDEPATH += 3rdparty/fftw-3.3.8/include

TARGET = petmrControl
TEMPLATE = app

SOURCES += main.cpp\
        ImageIO.cpp \
        anatomy.cpp \
        download.cpp \
        function.cpp \
        petmrMain.cpp \
        utilio.cpp

HEADERS  += petmrMain.h \
    ImageIO.h \
    io.h

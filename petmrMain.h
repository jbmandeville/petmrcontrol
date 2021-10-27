#ifndef PETMRMAIN_H
#define PETMRMAIN_H

#include <QMainWindow>
#include <QWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QProcess>
#include <QComboBox>
#include <QSettings>
#include <QListWidget>
#include <QDebug>

#include "io.h"

////////////////////////////////////////////////////
// Global variables:
// 1) download ID
// 2) download path
// 3) multi-subject template directory
// 4) smoothing
//
// addtional variables needed for anatomy page:
// 1) subject ID (saved a free-surfer directory)
//
// addtional variables needed for EPI page:
// 1) fMRI run list
//
// addtional variables needed for PET page:
//
//

////////////////////////////////////////////////////

class MainWindow : public QMainWindow
{
    Q_OBJECT

private:

    QString _fastmapProcess     = "/Users/jbm/QtApps/build-FM-Desktop_Qt_5_12_2_clang_64bit-Release/FM.app/Contents/MacOS/FM";
    QString _findsessionProcess = "/usr/pubsw/bin/findsession";
    QString _scriptDirectory    = "/homes/1/jbm/script/analyze-fm/";

    enum tabPages
    {
        page_download,
        page_anatomy,
        page_fMRI
    };
    struct FourDFile
    {
        QString name;
        iPoint4D dim;
    };


    QWidget *_centralWidget;
    QTabWidget *_tabs;
    QWidget *_downLoadPage;
    QWidget *_anatomyPage;
    QWidget *_fmriPage;
    QStatusBar *_statusBar;

    // Download page
    QLineEdit *_subjectIDDownload;
    QLabel *_downloadID;
    QLineEdit *_downloadPath;
    QPushButton *_downloadQueryButton;
    QPushButton *_downloadDataButton;

    // Anatomy page
    QComboBox *_anatomyInputDirectoryBox; // "003 004"
    QLineEdit *_subjectIDFreeSurfer;
    QComboBox *_anatomyFileNameBox;       // "raw.nii or "brain.nii"
    QLabel *_freeSurferInputFile;      // "t1/004/raw.nii"
    QLabel *_anatomyInputFile;         // "t1/004/brain.nii"
    QLabel *_anatomyTemplateDirectory; // "t1/004/brain.nii"
    QPushButton *_runFreeSurferButton;
    QPushButton *_selectTemplateDirectoryButton;
    QPushButton *_alignAnatomyToTemplateButton;

    // fMRI page
    QListWidget *_fMRIRunItemsBox;
    QVector<QListWidgetItem> _fMRIRunItems;
    QComboBox *_fMRITemplateDirectoryBox; // an existing directory
    QComboBox *_fMRIFileNameBox;          // "raw.nii or "mc.nii"
    QPushButton *_resliceEPIButton;
    QPushButton *_motionCorrectEPIButton;
    QPushButton *_alignEPIButton;

    // non-GUI variables
    QSettings _savedQSettings;        // needs organization & application name to work (see main.cpp)
    QString _lastTemplateDirectory;
    QVector<FourDFile> _fMRIFiles;
    iPoint4D _dimEPITemplate;

    void createDownloadPage();
    void createAnatomyPage();
    void createfMRIPage();
    void openedAnatomyPage();
    void openedfMRIPage();
    void readQSettings();
    void writeQSettings();
    QString getDimensions(QString fileName, iPoint4D &dim);
    void enableEPIActionButtons();

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
private slots:
    inline void enableGUI(int exitCode, QProcess::ExitStatus exitStatus )
    {
        qInfo() << "exit code" << exitCode << "exit status" << exitStatus;
        _centralWidget->setEnabled(true);
    }
    void queryDownloadPaths();
    void downloadData();
    void aboutApp();
    void exitApp();
    void changedPage(int index);
    void alignAnatomyToTemplate();
    void finishedFMAnatomyAlignment(int exitCode, QProcess::ExitStatus exitStatus );
    void selectTemplateDirectory();

    void resliceEPI();
    //    void motionCorrectEPI();
    void alignEPI();
    void changefMRITemplateDirectory(int indexInBox);
    void finishedFMResliceEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void finishedFMAlignEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void changedfMRIRunCheckBox(QListWidgetItem* item);

};

#endif // PETMRMAIN_H

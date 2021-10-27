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
// 1) download ID   (could be take from log file in principle)
// 2) download path (could be take from log file in principle)
// 3) multi-subject template directory (currently qsettings; should be local)
// 4) smoothing
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
    enum scanCategory
    {
        category_scout,
        category_T1,
        category_EPI,
        category_PET,
        category_MRAC,
        category_VIBE,
        category_UTE,
        category_UNKNOWN
    };
    struct FourDFile
    {
        QString name;
        iPoint4D dim;
    };
    struct programVariables
    {
        QString downloadID;
        QString downloadPath;
        QString subjectID;
    } _savedToDisk;
    struct downloadScan
    {
        QString scanNumber;
        QString sequenceName;
        iPoint4D dim;
        scanCategory category;
        QString categoryName;
        bool selectedForDownload;
    };


    QWidget *_centralWidget;
    QTabWidget *_tabs;
    QWidget *_downLoadPage;
    QWidget *_anatomyPage;
    QWidget *_fmriPage;
    QStatusBar *_statusBar;

    // Download page
    QLineEdit *_subjectIDDownload;
    QComboBox *_downloadIDBox;
    QComboBox *_downloadPathBox;
    QListWidget *_scanItemsBox;
    QVector<QListWidgetItem> _scanItems;
    QPushButton *_querySubjectButton;
    QPushButton *_generateScanListButton;
    QPushButton *_readAvailableScanList;
    QPushButton *_downloadDataButton;

    // Anatomy page
    QStringList _FastmapMSTemplateDirectories;
    QComboBox *_anatomyInputDirectoryBox; // "003 004"
    QLineEdit *_subjectIDFreeSurfer;
    QComboBox *_anatomyFileNameBox;       // "raw.nii or "brain.nii"
    QLabel *_freeSurferInputFile;      // "t1/004/raw.nii"
    QLabel *_anatomyInputFile;         // "t1/004/brain.nii"
    QComboBox *_anatomyTemplateDirectory; // multi-subject template directory
    QPushButton *_runFreeSurferButton;
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
    QVector<downloadScan> _scans;
    QString _lastTemplateDirectory;   // this should be written locally, not read from qsettings
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
    void outputConfigurationFile();
    void getSubjectNameFromFreeDir();

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
private slots:
    inline void enableGUI(int exitCode, QProcess::ExitStatus exitStatus )
    {
        qInfo() << "exit code" << exitCode << "exit status" << exitStatus;
        _centralWidget->setEnabled(true);
    }
    inline void changedDownloadIDBox(int indexInBox)   {_downloadPathBox->setCurrentIndex(indexInBox);}
    inline void changedDownloadPathBox(int indexInBox) {_downloadIDBox->setCurrentIndex(indexInBox);}
    void finishedGeneratingScanList(int exitCode, QProcess::ExitStatus exitStatus);


    void queryDownloadPaths();
    void generateScanList();
    void readAvailableScanList();
    void downloadData();
    void aboutApp();
    void exitApp();
    void changedPage(int index);
    void alignAnatomyToTemplate();
    void finishedFMAnatomyAlignment(int exitCode, QProcess::ExitStatus exitStatus );

    void resliceEPI();
    //    void motionCorrectEPI();
    void alignEPI();
    void changefMRITemplateDirectory(int indexInBox);
    void finishedFMResliceEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void finishedFMAlignEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void changedfMRIRunCheckBox(QListWidgetItem* item);

};

#endif // PETMRMAIN_H

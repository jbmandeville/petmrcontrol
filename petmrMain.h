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
#include <QGroupBox>
#include <QTextBrowser>
#include <QAction>
#include <QDebug>

#include "io.h"

////////////////////////////////////////////////////
// Global variables:
// 1) multi-subject template directory (currently qsettings; should be local)
// 2) smoothing
//
// additional variables needed for EPI page
// - probably not necessary to keep template averaging range (mcTemplate will be kept)
// - get template directory from mcTemplate.nii
// - get fMRI run list from location of mc.nii or align.nii
//
// PET page needs to know
// 1) list of EPI runs with times overlapping with PET
// 2) show PET frames with overlap
//

////////////////////////////////////////////////////

class MainWindow : public QMainWindow
{
    Q_OBJECT

private:

    QString _fastmapProcess     = "/homes/1/jbm/space1/dev/build-FM-Desktop_Qt_5_6_2_GCC_64bit-Release/FM";
    QString _findsessionProcess = "/usr/pubsw/bin/findsession";
    QString _dicomDumpProcess   = "/usr/bin/dcmdump";
    QString _scriptDirectory    = "/homes/1/jbm/script/analyze-fm/";

    enum tabPages_fMRIFiles
    {
        page_download,
        page_anatomy,
        page_fMRI,
        page_PET
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
        QString name;       // "epi/004/raw.nii"
        iPoint4D dim;       // x,y,z,t
        dVector timeTags;   // seconds[dim.t]
        sVector timeText;   // string[dim.t]
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
        QString scanNumberNew;  // "1" --> "001"
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
    QWidget *_petPage;
    QStatusBar *_statusBar;
    QAction *_browserAction;

    QTextBrowser *_outputBrowser;

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
    QGroupBox *_queryDownloadGroupBox;

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
    QGroupBox *_freeSurferGroupBox;
    QGroupBox *_anatomyInputBox;
    QGroupBox *_anatomyAlignmentBox;

    // fMRI page
    QListWidget *_fMRIRunItemBox;
    QVector<QListWidgetItem> _fMRIRunItems;
    QComboBox *_fMRITemplateDirectoryBox; // an existing directory
    QComboBox *_fMRIFileNameBox;          // "raw.nii or "mc.nii"
    QLineEdit *_fMRIMCRange;              // e.g. 1-10
    QPushButton *_resliceEPIButton;
    QPushButton *_motionCorrectEPIButton;
    QPushButton *_alignEPIButton;
    QGroupBox *_fMRIActionBox;

    // PET page
    QComboBox *_petRunBox;
    QListWidget *petFramesBox;
    QVector<QListWidgetItem> _petFrameItems;
    QLabel *_fMRIForPETTemplate;
    QLabel *_fMRIForPETFileName;
    QPushButton *_motionCorrectMatchingMRIButton;
    QPushButton *_motionCorrectPETButton;

    // non-GUI variables
    QSettings _savedQSettings;        // needs organization & application name to work (see main.cpp)
    QVector<downloadScan> _scans;
    QString _lastTemplateDirectory;   // this should be written locally, not read from qsettings
    QVector<FourDFile> _fMRIFiles;
    QVector<FourDFile> _fMRIFilesForPETMC;
    FourDFile _petFile;
    iPoint4D _dimEPITemplate;
    d2Matrix _matchingEPI;  // [_petFile.dim.t][list of pairs (fMRI file,time point)]

    void createDownloadPage();
    void createAnatomyPage();
    void createfMRIPage();
    void createPETPage();
    void openedAnatomyPage();
    void openedfMRIPage();
    void openedPETPage();
    void readQSettings();
    void writeQSettings();
    void enableEPIActionButtons();
    void outputDownloadList();
    void getSubjectNameFromFreeDir();
    void readAvailableScanList();
    void readUnpackLog();
    bool enableDownloadData();
    void reformatAcquisitionTimes(downloadScan scan);
    void updateFileNameBox();
    void findPETandFMRIOverlap();

    QString getDimensions(QString fileName, iPoint4D &dim);
    void getTimeTags(QString fileName, dVector &timeTags, sVector &timeText );
    QString twoDigits(short time);
    dPoint2D petFrameTime(int iFrame);
    void writeJipCommandFileForMatchingMRI();


public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
private slots:
    inline void enableGUI(int exitCode, QProcess::ExitStatus exitStatus )
    {
        qInfo() << "exit code" << exitCode << "exit status" << exitStatus;
        _centralWidget->setEnabled(true);
        showBrowser(false);
    }
    inline void outputToBrowser()
    {
        QProcess *process = qobject_cast<QProcess*>(sender());
        _outputBrowser->append(process->readAllStandardOutput());
    }
    inline void showBrowser(bool show)
    {
        FUNC_ENTER << show;
        if ( show )
            _outputBrowser->show();
        else
            _outputBrowser->hide();
    }

    inline void changedDownloadIDBox(int indexInBox)   {_downloadPathBox->setCurrentIndex(indexInBox);}
    inline void changedDownloadPathBox(int indexInBox) {_downloadIDBox->setCurrentIndex(indexInBox);}
    void finishedGeneratingScanList(int exitCode, QProcess::ExitStatus exitStatus);
    void finishedDownloadData(int exitCode, QProcess::ExitStatus exitStatus);

    void queryDownloadPaths();
    void generateScanList();
    void downloadData();
    void finishedDownloadDataNew();
    void aboutApp();
    void exitApp();
    void changedPage(int index);
    void alignAnatomyToTemplate();
    void finishedFMAnatomyAlignment(int exitCode, QProcess::ExitStatus exitStatus );

    void resliceEPI();
    void motionCorrectEPI();
    void alignEPI();
    void changefMRITemplateDirectory(int indexInBox);
    void changedfMRIFileName(int indexInBox);
    void finishedFMResliceEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void finishedFMAlignEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void changedfMRIRunCheckBox(QListWidgetItem* item);
    void changedDownloadScanCheckBox(QListWidgetItem *item);
    
    void finishedMotionCorrectEPI(int exitCode, QProcess::ExitStatus exitStatus);
    
    void updatePETRunBox(int indexInBox);
    void changedPETFrameSelection(QListWidgetItem *item);
    void motionCorrectMatchingMRI();
//    void motionCorrectPET();
    void finishedMotionCorrectMatchingMRI(int exitCode, QProcess::ExitStatus exitStatus);

};

#endif // PETMRMAIN_H

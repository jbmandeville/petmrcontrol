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
#include <QRadioButton>
#include <QDebug>

#include "io.h"

////////////////////////////////////////////////////
// Global variables:
// 1) multi-subject template directory (currently qsettings; should be local)
// 2) smoothing
//
// TO DO:
// 1) PET page: need to make file type variable (in case "raw.nii" is deleted)
// 2) EPI page: after reslice, update file list (now have to leave page and go back)
// 3) anatomy page: update file names after runFreeSurfer
// 4)

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
        page_PET,
        page_clean
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
        bool existsOnDisk;
        bool selectedForDownload;
        bool selectedForCleaning;
    };
    struct savedSettings
    {
        // Image Window
        QByteArray imageWindowGeometry;
        QByteArray browserWindowGeometry;
        QString lastTemplateDirectory;   // this should be written locally, not read from qsettings
    };

    QWidget *_centralWidget;
    QTabWidget *_tabs;
    QWidget *_downLoadPage;
    QWidget *_anatomyPage;
    QWidget *_fmriPage;
    QWidget *_petPage;
    QWidget *_cleanPage;
    QStatusBar *_statusBar;
    QAction *_browserAction;

    QTextBrowser *_outputBrowser;

    // Radiobuttons
    QRadioButton *_radioButtonHumanBay7;
    QRadioButton *_radioButtonNHPBay6;
    QRadioButton *_radioButtonNHPBay7;

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
    QGroupBox *_freeSurferGroupBox;
    QStringList _FastmapMSTemplateDirectories;
    QComboBox *_anatomyDirBox; // "003 004"
    QLineEdit *_subjectIDFreeSurfer;
    QComboBox *_anatomyFileNameBox;       // "raw.nii or "brain.nii"
    QComboBox *_anatomyTemplateDirectory; // multi-subject template directory
    QPushButton *_runFreeSurferButton;
    QPushButton *_alignAnatomyButton;
    QPushButton *_extractFreeSurferOverlaysButton;

    // fMRI page
    QListWidget *_fMRIRunItemBox;
    QVector<QListWidgetItem> _fMRIRunItems;
    QComboBox *_fMRITemplateDirBox; // an existing directory
    QComboBox *_fMRIFileNameBox;          // "raw.nii or "mc.nii"
    QLineEdit *_fMRIMCRange;              // e.g. 1-10
    QPushButton *_resliceEPIButton;
    QPushButton *_motionCorrectEPIButton;
    QPushButton *_alignEPIButton;

    // PET page
    QComboBox *_petDirBox;
    QListWidget *_petFramesBox;
    QVector<QListWidgetItem> _petFrameItems;
    QLabel *_fMRIForPETTemplate;
    QLabel *_fMRIForPETFileName;
    QString _anatomyFileNameForPETReslice;
    QString _alignFileNameForPETRegistration;
    QPushButton *_motionCorrectMatchingMRIButton;
    QPushButton *_motionCorrectPETButton;
    QPushButton *_reslicePETButton;
    QPushButton *_alignPETButton;
    QPushButton *_analyzeTAC;

    // clean page
    QListWidget *_cleanScanTypesBox;
    QVector<QListWidgetItem> _cleanScanTypeItems;
    QPushButton *_cleanDICOMs;
    QPushButton *_cleanNII_auxilliary;  // MR*.nii, test.nii, matchingMRI-rs.nii, matchingMRI-mc.nii
    QPushButton *_cleanNII_mc;          // mc.nii: this needs to be kept in order to change smoothing
    QPushButton *_cleanNII_raw;         // raw.nii
    QLabel *_totalSizeSubDirs;

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
    void createCleanPage();
    void openedAnatomyPage();
    void openedfMRIPage();
    void openedPETPage();
    void openedCleanPage();
    void readQSettings();
    void writeQSettings();
    void outputDownloadList();
    void getSubjectNameFromFreeDir();
    void readAvailableScanList();
    void readUnpackLog();
    bool enableDownloadData();
    void reformatAcquisitionTimes(downloadScan scan);
    void updateFileNameBox();
    void findPETandFMRIOverlap();
    void setupScanTypes();
    void updateCleaningList();
    void updateAnatomyFileName();
    void installSRTMAnalysis();
    void writeFramesTable(QString fileName);
    void writeTimeModelFile(QString directoryName);
    void writeGLMFile(QString directoryName);

    QString getDimensions(QString fileName, iPoint4D &dim);
    void getTimeTags(QString fileName, dVector &timeTags, sVector &timeText );
    QString twoDigits(short time);
    dPoint2D petFrameTime(int iFrame);
    void writeJipCommandFileForMatchingMRI();
    bool getPETMCInterpolationRequired();

    inline bool anatomyFileExists(QString fileName) {return anatomyFileExists(_anatomyDirBox->currentText(),fileName);}
    inline bool petFileExists(QString fileName)     {return petFileExists(_petDirBox->currentText(),fileName);}
    inline bool epiFileExists(QString fileName)     {return epiFileExists(_fMRITemplateDirBox->currentText(),fileName);}
    bool anatomyFileExists(QString dirName, QString fileName);
    bool epiFileExists(QString dirName, QString fileName);
    bool petFileExists(QString dirName, QString fileName);

    void enableEPIActionButtons();
    void enableAnatomyActionButtons();
    void enablePETActionButtons();

    void findDICOMs(bool remove);
    void findAuxFiles(bool remove);
    void findMCFiles(bool remove);
    void findRawFiles(bool remove);
    void findAllFiles();

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
private slots:
    void dataOriginChanged();

    inline void enableGUI(int exitCode, QProcess::ExitStatus exitStatus )
    {
        qInfo() << "finished";
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
    inline void changedAnatomyFileName(int indexInBox) {enableAnatomyActionButtons();}
    void changedAnatomyDirName(int indexInBox);

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
    void extractFreeSurferOverlays();
    void finishedExtractOverlays(int exitCode, QProcess::ExitStatus exitStatus );
    void runFreeSurfer();
    void finishedRunFreeSurfer(int exitCode, QProcess::ExitStatus exitStatus );
    void displayAnatomy();

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
    void displayEPI();

    void updatePETRunBox(int indexInBox);
    void changedPETFrameSelection(QListWidgetItem *item);
    void motionCorrectMatchingMRI();
    void applyMotionCorrectionToPET();
    void reslicePET();
    void alignPET();
    void analyzeTAC();
    void finishedMotionCorrectMatchingMRI(int exitCode, QProcess::ExitStatus exitStatus);
    void finishedApplyingMCToPET(int exitCode, QProcess::ExitStatus exitStatus);
    void finishedFMReslicePET(int exitCode, QProcess::ExitStatus exitStatus );
    void finishedFMAlignPET(int exitCode, QProcess::ExitStatus exitStatus );
    void finishedFMAnalyzeTAC(int exitCode, QProcess::ExitStatus exitStatus );

    void changedScanTypeCheckBox(QListWidgetItem *item);
    inline void cleanDICOMFiles()     {findDICOMs(true);         findDICOMs(false);         findAllFiles();}
    inline void cleanMCFiles()        {findMCFiles(true);        findMCFiles(false);        findAllFiles();}
    inline void cleanAuxNIIFiles()    {findAuxFiles(true);       findAuxFiles(false);       findAllFiles();}
    inline void cleanRawNIIFiles()    {findRawFiles(true);       findRawFiles(false);       findAllFiles();}

};

#endif // PETMRMAIN_H

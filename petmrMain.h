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
#include <QStatusBar>

#include "io.h"

////////////////////////////////////////////////////
// Global variables:
//
// TO DO:

////////////////////////////////////////////////////

class MainWindow : public QMainWindow
{
    Q_OBJECT

private:

//    QString _fastmapProcess     = "/homes/1/jbm/space1/dev/build-FM-Desktop_Qt_5_6_2_GCC_64bit-Release/FM";
    QString _fastmapProcess     = "/usr/pubsw/packages/jip/bin/Linux-x86_64/FM";
    QString _findsessionProcess = "/usr/pubsw/bin/findsession";
    QString _dicomDumpProcess   = "/usr/bin/dcmdump";
    QString _scriptDirectory    = "/space/deltabp/1/users/public/script/analyze-petmr/";

    enum tabPages
    {
        page_download,
        page_anatomy,
        page_fMRI,
        page_PET,
        page_clean,
        page_size
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
        bool isAntomicalT1=false;
        bool isMRAC=false;
        bool existsOnDisk;
        bool selectedForDownload;
        bool selectedForCleaning;
    };
    struct savedSettings
    {
        // Image Window
        QByteArray imageWindowGeometry;
        QByteArray browserWindowGeometry;
        double fmSmoothing=0.;
    };

    // time tag corrections
    int _EPITimeCorrection=0;
    int _PETTimeCorrection=0;

    QWidget *_centralWidget;
    QTabWidget *_tabs;
    QWidget *_downLoadPage;
    QWidget *_anatomyPage;
    QWidget *_fmriPage;
    QWidget *_petPage;
    QWidget *_cleanPage;
    QStatusBar *_statusBar;
    QVector<QTextEdit *> _noteBox;
    QAction *_showNotesAction;

    QTextBrowser *_outputBrowser;
    QWidget *_helpTool;
    QTextBrowser *_helpBrowser;
    int _helpPageIndex=0;

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
    QPushButton *_setupTransferDirectory;
    QGroupBox *_queryDownloadGroupBox;

    // Anatomy page
    QGroupBox *_freeSurferGroupBox;
    QStringList _FastmapMSTemplateDirectories;
    QComboBox *_anatomyDirBox; // "003 004"
    QLineEdit *_subjectIDFreeSurfer;
    QComboBox *_anatomyFileBox;       // "raw.nii or "brain.nii"
    QComboBox *_anatomyTemplateDirectory; // multi-subject template directory
    QPushButton *_runFreeSurferButton;
    QPushButton *_alignAnatomyButton;
    QPushButton *_extractFreeSurferOverlaysButton;
    QLineEdit *_smoothingAnatomy;

    // fMRI page
    QListWidget *_fMRIRunItemBox;
    QVector<QListWidgetItem> _fMRIRunItems;
    QComboBox *_fMRITemplateDirBox; // an existing directory
    QComboBox *_fMRIFileBox;          // "raw.nii or "mrangeLabelc.nii"
    QLineEdit *_fMRIMCRange;              // e.g. 1-10
    QPushButton *_doEverythingEPIButton;
    QPushButton *_resliceEPIButton;
    QPushButton *_motionCorrectEPIButton;
    QPushButton *_alignEPIButton;
    QLineEdit *_smoothingfMRI;

    // PET page
    QComboBox *_petDirBox;
    QComboBox *_petFileBox;
    QListWidget *_petFramesBox;
    QVector<QListWidgetItem> _petFrameItems;
    QLabel *_fMRIForPETTemplate;
    QLabel *_fMRIForPETFileName;
    QString _alignFileNameForPETRegistration;
    QLineEdit *_smoothingPET;
    QPushButton *_doEverythingPETButton;
    QPushButton *_motionCorrectMatchingMRIButton;
    QPushButton *_motionCorrectPETButton;
    QPushButton *_alignPETButton;
    QPushButton *_analyzeTAC;
    QGroupBox *_mcPETBox;

    // clean page
    QLineEdit *_correctEPITimeTags;
    QLineEdit *_correctPETTimeTags;
    QListWidget *_cleanScanTypesBox;
    QVector<QListWidgetItem> _cleanScanTypeItems;
    QPushButton *_cleanDICOMs;
    QPushButton *_cleanNII_auxilliary;  // MR*.nii, test.nii, matchingMRI-rs.nii, matchingMRI-mc.nii
    QPushButton *_cleanNII_mc;          // mc.nii: this needs to be kept in order to change smoothing
    QPushButton *_cleanNII_raw;         // raw.nii
    QLabel *_totalSizeSubDirs;

    // non-GUI variables
    QVector<downloadScan> _scans;
    QVector<FourDFile> _fMRIFiles;
    QVector<FourDFile> _fMRIFilesForPETMC;
    FourDFile _petFile;
    iPoint4D _dimEPITemplate;
    i2Matrix _matchingEPI;  // [_petFile.dim.t][list of pairs (fMRI file,time point)]
    savedSettings _savedSettings;
    QSettings _savedQSettings;

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
    QString readFileTextArgument(QString fileName, QString parameterName);
    bool enableDownloadData();
    bool enableTransferDirectory();
    void reformatAcquisitionTimes(downloadScan scan);
    void readSubjectVariables();
    void readSmoothing(int which);
    void writeAllNotes();
    void loadNotes();
    void loadHelp(int whichTab);
    int whichTabName(QString name);
    void setupScanTypes();
    void copyDICOMs(bool T1, QString destinationDir);

    void updateAnatomyFileName();
    void setTemplate();

    void populateEPIFileNameBox();
    void setDefaultIndexEPIFileNameBox();
    void updateEPIFileNameBox(QString fileName);
    void writeJipCommandFileForMCAveraging();
    QString readJipCommandFileForMCAveraging();

    void updatePETFileNameBox(QString fileName);
    void updatePETFileNameBox();
    void findPETandFMRIOverlap();
    void installSRTMAnalysis();
    void writeFramesTable(QString fileName);
    void writeTimeModelFile(QString directoryName);
    void writeGLMFile(QString directoryName);

    void updateCleaningList();

    QString getDimensions(QString fileName, iPoint4D &dim);
    void getTimeTags(QString fileName, int correction, dVector &timeTags, sVector &timeText );
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
    QString alignCOMFileName(int which);

    inline void spawnProcess(QProcess *process, QString exe, QStringList arguments,
                             QString message, QString browserTitle )
    {
        _statusBar->showMessage(message);
        _centralWidget->setEnabled(false);
        if ( !browserTitle.isEmpty() )
        {
            _outputBrowser->setWindowTitle("Query Progress");
            _outputBrowser->show();
        }
        qInfo() <<  exe << arguments;
        process->start(exe,arguments);
    }

    inline void finishedProcess()
    {
        _statusBar->clearMessage();
        showOutputBrowser(false);
        _centralWidget->setEnabled(true);
    }

    void enableEPIActionButtons();
    void enableAnatomyActionButtons();

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
        finishedProcess();
    }
    inline void outputToBrowser()
    {
        QProcess *process = qobject_cast<QProcess*>(sender());
        _outputBrowser->append(process->readAllStandardOutput());
    }
    inline void showOutputBrowser(bool show)
    {
        if ( show ) _outputBrowser->show();
        else        _outputBrowser->hide();
    }
    inline void showHelpBrowser(bool show)
    {
        if ( show ) _helpTool->show();
        else        _helpTool->hide();
    }
    void showNone();
    void showNotes(bool show);
    void helpGoBackward();
    void helpGoForward();

    inline void changedAnatomyFileName(int indexInBox) {enableAnatomyActionButtons();}
    void changedAnatomyDirName(int indexInBox);

    inline void changedDownloadIDBox(int indexInBox)   {_downloadPathBox->setCurrentIndex(indexInBox);}
    inline void changedDownloadPathBox(int indexInBox) {_downloadIDBox->setCurrentIndex(indexInBox);}
    void finishedGeneratingScanList(int exitCode, QProcess::ExitStatus exitStatus);
    void finishedDownloadData(int exitCode, QProcess::ExitStatus exitStatus);

    void queryDownloadPaths();
    void generateScanList();
    void downloadData();
    void setupTransferDirectory();
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
    void writeSubjectVariables();

    void doEverthingEPI();
    void resliceEPI();
    void motionCorrectEPI();
    void alignEPI();
    void changefMRITemplateDirectory(int indexInBox);
    void changedfMRIFileName(int indexInBox);
    void finishedFMResliceEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void finishedLinkResliceEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void finishedFMAlignEPI(int exitCode, QProcess::ExitStatus exitStatus );
    void changedfMRIRunCheckBox(QListWidgetItem* item);
    void changedDownloadScanCheckBox(QListWidgetItem *item);
    void changedSmoothingAnatomy();
    void changedSmoothingfMRI();
    void changedSmoothingPET();
    void changedEPITimeTagCorrection();
    void changedPETTimeTagCorrection();

    void finishedMotionCorrectEPI(int exitCode, QProcess::ExitStatus exitStatus);
    void displayEPI();

    void updatePETDirBox(int indexInBox);
    void changedPETFrameSelection(QListWidgetItem *item);
    void motionCorrectMatchingMRI();
    void applyMotionCorrectionToPET();
    void doEverthingPET();
    void alignPET();
    void analyzeTAC();
    void enablePETActionButtons();
    void finishedMotionCorrectMatchingMRI(int exitCode, QProcess::ExitStatus exitStatus);
    void finishedApplyingMCToPET(int exitCode, QProcess::ExitStatus exitStatus);
    void finishedFMAlignPET(int exitCode, QProcess::ExitStatus exitStatus );
    void finishedFMAnalyzeTAC(int exitCode, QProcess::ExitStatus exitStatus );

    void changedScanTypeCheckBox(QListWidgetItem *item);
    inline void cleanDICOMFiles()     {findDICOMs(true);         findDICOMs(false);         findAllFiles();}
    inline void cleanMCFiles()        {findMCFiles(true);        findMCFiles(false);        findAllFiles();}
    inline void cleanAuxNIIFiles()    {findAuxFiles(true);       findAuxFiles(false);       findAllFiles();}
    inline void cleanRawNIIFiles()    {findRawFiles(true);       findRawFiles(false);       findAllFiles();}

};

#endif // PETMRMAIN_H

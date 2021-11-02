#include <QtWidgets>
#include "petmrMain.h"

void MainWindow::createAnatomyPage()
{
    _anatomyPage = new QWidget();
    auto *setupFreeSurferLayout = new QGridLayout();

    ////////////////////////////////////////////////
    // input directory
    ////////////////////////////////////////////////
    auto *anatDirLabel = new QLabel("Sub-directory on 't1'");
    auto *inputLayout = new QGridLayout();
    _anatomyInputDirectoryBox = new QComboBox();
    inputLayout->addWidget(anatDirLabel,0,0);
    inputLayout->addWidget(_anatomyInputDirectoryBox,0,1);
    _anatomyInputBox = new QGroupBox("Input directory for anatomy");
    _anatomyInputBox->setLayout(inputLayout);

    ////////////////////////////////////////////////
    // freeSurfer
    ////////////////////////////////////////////////
    auto *subjectIDLabel     = new QLabel("Subject ID",_anatomyPage);
    auto *freeInputFileLabel = new QLabel("Input file",_anatomyPage);
    _subjectIDFreeSurfer     = new QLineEdit("?");
    _freeSurferInputFile     = new QLabel();
    setupFreeSurferLayout->addWidget(subjectIDLabel,0,0);
    setupFreeSurferLayout->addWidget(_subjectIDFreeSurfer,0,1);
    setupFreeSurferLayout->addWidget(freeInputFileLabel,1,0);
    setupFreeSurferLayout->addWidget(_freeSurferInputFile,1,1);

    _runFreeSurferButton       = new QPushButton("Run FreeSurfer",_downLoadPage);
//    connect(_runFreeSurferButton, SIGNAL(pressed()), this, SLOT(runFreeSurfer()));

    auto *freeSurferLayout = new QVBoxLayout();
    freeSurferLayout->addLayout(setupFreeSurferLayout);
    freeSurferLayout->addWidget(_runFreeSurferButton);
    freeSurferLayout->setSpacing(0);

    _freeSurferGroupBox = new QGroupBox("Run FreeSurfer to delineate anatomy (takes HOURS)");
    _freeSurferGroupBox->setLayout(freeSurferLayout);

    ////////////////////////////////////////////////
    // anatomy registration
    ////////////////////////////////////////////////
    auto *anatInputFileLabel  = new QLabel("Input file: ",_anatomyPage);
    auto *templateDirLabel    = new QLabel("Template directory: ",_anatomyPage);
    _anatomyFileNameBox       = new QComboBox();
    _anatomyInputFile         = new QLabel("Input file = ?");
    _anatomyTemplateDirectory = new QComboBox();
    int iSelection = 0;
    for (int jList=0; jList<_FastmapMSTemplateDirectories.count(); jList+=2)
    {
        _anatomyTemplateDirectory->addItem(_FastmapMSTemplateDirectories.at(jList));
        if ( !_lastTemplateDirectory.compare(_FastmapMSTemplateDirectories.at(jList)) )
            iSelection = jList/2;
    }
    _anatomyTemplateDirectory->setCurrentIndex(iSelection);

    auto *anatomyFileLayout = new QGridLayout();
    anatomyFileLayout->addWidget(anatInputFileLabel,0,0);
    anatomyFileLayout->addWidget(_anatomyFileNameBox,0,1);
    anatomyFileLayout->addWidget(templateDirLabel,1,0);
    anatomyFileLayout->addWidget(_anatomyTemplateDirectory,1,1);
    anatomyFileLayout->setSpacing(0);

    auto *anatomyAlignmentLayout = new QVBoxLayout();
    _alignAnatomyToTemplateButton  = new QPushButton("Align to template",_anatomyPage);
    connect(_alignAnatomyToTemplateButton, SIGNAL(pressed()), this, SLOT(alignAnatomyToTemplate()));

    anatomyAlignmentLayout->addLayout(anatomyFileLayout);
    anatomyAlignmentLayout->addWidget(_alignAnatomyToTemplateButton);

    _anatomyAlignmentBox = new QGroupBox("Use fastmap to align anatomy to a multi-subject template");
    _anatomyAlignmentBox->setLayout(anatomyAlignmentLayout);

    QString freeDir = "free";
    QFileInfo checkDir(freeDir);
    if (checkDir.exists() && checkDir.isDir())
        getSubjectNameFromFreeDir();    // get subject name

    ////////////////////////////////////////////////
    // full page layout
    ////////////////////////////////////////////////
    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(_anatomyInputBox);
    pageLayout->addWidget(_freeSurferGroupBox);
    pageLayout->addWidget(_anatomyAlignmentBox);
    pageLayout->setSpacing(0);
    _anatomyPage->setLayout(pageLayout);
}

void MainWindow::getSubjectNameFromFreeDir()
{
    // Read the time model file
    QDir const freeDir("./free");
    QStringList const folderList = freeDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot | QDir::NoSymLinks);
    FUNC_INFO << folderList;
    if ( folderList.count() == 1 )
    {
        _subjectIDDownload->setText(folderList.at(0));
        _subjectIDFreeSurfer->setText(folderList.at(0));
        _freeSurferGroupBox->setEnabled(false);
        _queryDownloadGroupBox->setEnabled(false);
        _anatomyAlignmentBox->setStyleSheet("background-color:lightYellow;");
        _freeSurferGroupBox->setStyleSheet("background-color:lightYellow;");
        setWindowTitle(QString("petmrcontrol: %1").arg(_subjectIDFreeSurfer->text()));
    }
}

void MainWindow::openedAnatomyPage()
{
    FUNC_ENTER;

    QDir const anatomyTopDir("./t1");
    if (!anatomyTopDir.exists())
    {
        _anatomyInputDirectoryBox->clear();
        return;
    }
    QStringList const folderList = anatomyTopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    qInfo() << folderList;
    _anatomyInputDirectoryBox->clear();
    for (int jList=0; jList<folderList.size(); jList++)
        _anatomyInputDirectoryBox->addItem(folderList.at(jList));
    _anatomyInputDirectoryBox->setCurrentIndex(_anatomyInputDirectoryBox->count()-1);

    // For the current selection, show all NIFTI files
    QString path = "./t1/" + _anatomyInputDirectoryBox->currentText();
    QDir anatomyDir(path);
    anatomyDir.setNameFilters(QStringList()<<"*.nii");
    QStringList fileList = anatomyDir.entryList();

    _anatomyFileNameBox->clear();
    int indexRaw=-1;  int indexBrain=-1;
    for (int jList=0; jList<fileList.size(); jList++)
    {
        _anatomyFileNameBox->addItem(fileList.at(jList));
        if ( fileList.at(jList) == "brain.nii") indexBrain = jList;
        if ( fileList.at(jList) == "raw.nii")   indexRaw   = jList;
    }
    if ( indexBrain >= 0 )
    {
        _anatomyFileNameBox->setCurrentIndex(indexBrain);
        _anatomyInputBox->setEnabled(false);
        _anatomyAlignmentBox->setEnabled(false);
        _anatomyInputBox->setStyleSheet("background-color:lightYellow;");
        _anatomyAlignmentBox->setStyleSheet("background-color:lightYellow;");
    }
    else if ( indexRaw >= 0 )
        _anatomyFileNameBox->setCurrentIndex(indexRaw);
    else
        _anatomyFileNameBox->setCurrentIndex(0);

    FUNC_INFO << 2;
    QString fileName = "t1/" + _anatomyInputDirectoryBox->currentText() + "/raw.nii";
    _freeSurferInputFile->setText(fileName);
    FUNC_INFO << 3;
    fileName = "t1/" + _anatomyInputDirectoryBox->currentText() + "/" + _anatomyFileNameBox->currentText();
    _anatomyInputFile->setText(fileName);

    // Enable/disable:  // if anatomy file was found, it can be used for either freeSurfer or alignment
    bool enable = _anatomyFileNameBox->count() > 0;
    _runFreeSurferButton->setEnabled(enable);
    fileName = _anatomyTemplateDirectory->currentText();
//    enable &= fileName.compare("unknown",Qt::CaseInsensitive);
    _alignAnatomyToTemplateButton->setEnabled(enable);
}

void MainWindow::alignAnatomyToTemplate()
{
    FUNC_ENTER;
    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedFMAnatomyAlignment(int, QProcess::ExitStatus)));

    QStringList arguments;
    arguments.append(_anatomyInputFile->text());
    arguments.append("-O");
    arguments.append("alignment");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());
    qInfo() << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    FUNC_EXIT;
}

void MainWindow::finishedFMAnatomyAlignment(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "exit code" << exitCode << "exit status" << exitStatus;
    _centralWidget->setEnabled(true);
}

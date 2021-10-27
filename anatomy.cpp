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
    auto *inputBox = new QGroupBox("Input directory for anatomy");
    inputBox->setLayout(inputLayout);

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

    auto *freeBox = new QGroupBox("Run FreeSurfer to delineate anatomy (takes HOURS)");
    freeBox->setLayout(freeSurferLayout);

    QString freeDir = "free";
    QFileInfo checkDir(freeDir);
    if (checkDir.exists() && checkDir.isDir())
    {
        // get subject name
        getSubjectNameFromFreeDir();
        freeBox->setEnabled(false);
    }

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

    auto *anatomyLayout = new QVBoxLayout();
    _alignAnatomyToTemplateButton  = new QPushButton("Align to template",_anatomyPage);
    connect(_alignAnatomyToTemplateButton, SIGNAL(pressed()), this, SLOT(alignAnatomyToTemplate()));

    anatomyLayout->addLayout(anatomyFileLayout);
    anatomyLayout->addWidget(_alignAnatomyToTemplateButton);

    auto *anatomyBox = new QGroupBox("Use fastmap to align anatomy to a multi-subject template");
    anatomyBox->setLayout(anatomyLayout);

    ////////////////////////////////////////////////
    // full page layout
    ////////////////////////////////////////////////
    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(inputBox);
    pageLayout->addWidget(freeBox);
    pageLayout->addWidget(anatomyBox);
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
        _subjectIDFreeSurfer->setText(folderList.at(0));
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
    // Save the old list
    QStringList oldFolderList;
    for (int jList=0; jList<_anatomyInputDirectoryBox->count(); jList++)
        oldFolderList.append(_anatomyInputDirectoryBox->itemText(jList));
    QStringList const folderList = anatomyTopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);

    // Did something change?
    bool changeOnDisk = _anatomyInputDirectoryBox->count() != folderList.size();
    if ( !changeOnDisk )
    {
        for (int jList=0; jList<_anatomyInputDirectoryBox->count(); jList++)
            changeOnDisk |= _anatomyInputDirectoryBox->itemText(jList) != folderList.at(jList);
    }
    if ( changeOnDisk )
    { // Something changed
        qInfo() << folderList;
        _anatomyInputDirectoryBox->clear();
        for (int jList=0; jList<folderList.size(); jList++)
            _anatomyInputDirectoryBox->addItem(folderList.at(jList));
        _anatomyInputDirectoryBox->setCurrentIndex(_anatomyInputDirectoryBox->count()-1);
    }

    FUNC_INFO << 1;
    // For the current selection, show all NIFTI files
    QString path = "./t1/" + _anatomyInputDirectoryBox->currentText();
    QDir anatomyDir(path);
    anatomyDir.setNameFilters(QStringList()<<"*.nii");
    QStringList fileList = anatomyDir.entryList();
    // Did something change?
    changeOnDisk = _anatomyFileNameBox->count() != fileList.size();
    if ( !changeOnDisk )
    {
        for (int jList=0; jList<_anatomyFileNameBox->count(); jList++)
            changeOnDisk |= _anatomyFileNameBox->itemText(jList) != fileList.at(jList);
    }
    if ( changeOnDisk )
    { // Something changed
        qInfo() << fileList;
        _anatomyFileNameBox->clear();
        int indexRaw=-1;  int indexBrain=-1;
        for (int jList=0; jList<fileList.size(); jList++)
        {
            _anatomyFileNameBox->addItem(fileList.at(jList));
            if ( fileList.at(jList) == "brain.nii") indexBrain = jList;
            if ( fileList.at(jList) == "raw.nii")   indexRaw   = jList;
        }
        if ( indexBrain >= 0 )
            _anatomyFileNameBox->setCurrentIndex(indexBrain);
        else if ( indexRaw >= 0 )
            _anatomyFileNameBox->setCurrentIndex(indexRaw);
        else
            _anatomyFileNameBox->setCurrentIndex(0);
    }
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
    auto *view = new QTextBrowser;
    QObject::connect(process, &QProcess::readyReadStandardOutput, [process,view]()
    {
        auto output=process->readAllStandardOutput();
        view->append(output);
    });
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

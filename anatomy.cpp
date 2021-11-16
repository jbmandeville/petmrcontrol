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
    auto *anatInputFileLabel  = new QLabel("Input file: ",_anatomyPage);
    auto *inputLayout = new QGridLayout();
    _anatomyInputDirectoryBox = new QComboBox();
    _anatomyFileNameBox  = new QComboBox();
    connect(_anatomyInputDirectoryBox, SIGNAL(currentIndexChanged(int)),this, SLOT(changedAnatomyDirName(int)));
    connect(_anatomyFileNameBox, SIGNAL(activated(int)),this, SLOT(changedAnatomyFileName(int)));

    inputLayout->addWidget(anatDirLabel,0,0);
    inputLayout->addWidget(_anatomyInputDirectoryBox,0,1);
    inputLayout->addWidget(anatInputFileLabel,1,0);
    inputLayout->addWidget(_anatomyFileNameBox,1,1);

    auto *anatomyInputBox = new QGroupBox("Input directory for anatomy");
    anatomyInputBox->setLayout(inputLayout);
    anatomyInputBox->setStyleSheet("border: 1px dotted gray");

    ////////////////////////////////////////////////
    // freeSurfer
    ////////////////////////////////////////////////
    auto *subjectIDLabel     = new QLabel("Subject ID",_anatomyPage);
    _subjectIDFreeSurfer     = new QLineEdit("");
    setupFreeSurferLayout->addWidget(subjectIDLabel,0,0);
    setupFreeSurferLayout->addWidget(_subjectIDFreeSurfer,0,1);

    _runFreeSurferButton       = new QPushButton("Run FreeSurfer (raw --> brain)",_anatomyPage);
    connect(_runFreeSurferButton, SIGNAL(pressed()), this, SLOT(runFreeSurfer()));

    auto *freeSurferLayout = new QVBoxLayout();
    freeSurferLayout->addLayout(setupFreeSurferLayout);
    freeSurferLayout->addWidget(_runFreeSurferButton);
    freeSurferLayout->setSpacing(0);

    _freeSurferGroupBox = new QGroupBox("Run FreeSurfer to delineate anatomy (TAKES HOURS)");
    _freeSurferGroupBox->setLayout(freeSurferLayout);

    ////////////////////////////////////////////////
    // anatomy registration
    ////////////////////////////////////////////////
    auto *templateDirLabel    = new QLabel("Template directory: ",_anatomyPage);
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
    anatomyFileLayout->addWidget(templateDirLabel,0,0);
    anatomyFileLayout->addWidget(_anatomyTemplateDirectory,0,1);
    anatomyFileLayout->setSpacing(0);

    auto *anatomyAlignmentLayout = new QVBoxLayout();
    _alignAnatomyButton  = new QPushButton("Align to template (raw/brain --> align)",_anatomyPage);
    connect(_alignAnatomyButton, SIGNAL(pressed()), this, SLOT(alignAnatomyToTemplate()));
    _extractFreeSurferOverlaysButton  = new QPushButton("extract freeSurfer ROIs from atlas file)",_anatomyPage);
    connect(_extractFreeSurferOverlaysButton, SIGNAL(pressed()), this, SLOT(extractFreeSurferOverlays()));

    anatomyAlignmentLayout->addLayout(anatomyFileLayout);
    anatomyAlignmentLayout->addWidget(_alignAnatomyButton);
    anatomyAlignmentLayout->addWidget(_extractFreeSurferOverlaysButton);

    auto *anatomyAlignmentBox = new QGroupBox("Use fastmap to align anatomy to a multi-subject template");
    anatomyAlignmentBox->setLayout(anatomyAlignmentLayout);
    anatomyAlignmentBox->setStyleSheet("border: 1px dotted gray");

    QString freeDir = "free";
    QFileInfo checkDir(freeDir);
    if (checkDir.exists() && checkDir.isDir())
        getSubjectNameFromFreeDir();    // get subject name

    ////////////////////////////////////////////////
    // full page layout
    ////////////////////////////////////////////////
    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(anatomyInputBox);
    pageLayout->addWidget(_freeSurferGroupBox);
    pageLayout->addWidget(anatomyAlignmentBox);
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
        QDir thisDir = QDir::currentPath();
        QStringList subDirs = thisDir.absolutePath().split("/");
        FUNC_INFO << subDirs;
        int nList = qMin(3,subDirs.count());
        FUNC_INFO << "nList" << nList;
        QString list;  list.append("[..]/");
        for (int jList=nList; jList>0; jList--)
        {
            list.append(subDirs.at(subDirs.count()-jList));
            if ( jList != 1 ) list.append("/");
        }
        _queryDownloadGroupBox->setEnabled(false);
        setWindowTitle(QString("subject %1 @ %2").arg(_subjectIDFreeSurfer->text()).arg(list));
    }
    enableAnatomyActionButtons();
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
    FUNC_INFO << folderList;
    _anatomyInputDirectoryBox->clear();
    for (int jList=0; jList<folderList.size(); jList++)
        _anatomyInputDirectoryBox->addItem(folderList.at(jList));
    _anatomyInputDirectoryBox->setCurrentIndex(_anatomyInputDirectoryBox->count()-1);

//    changedAnatomyDirName(_anatomyInputDirectoryBox->currentIndex());
}

void MainWindow::changedAnatomyDirName(int indexInBox)
{
    updateAnatomyFileName();
}

void MainWindow::updateAnatomyFileName()
{
    // For the current selection, show all NIFTI files
    QString path = "./t1/" + _anatomyInputDirectoryBox->currentText();
//    QString path = "./t1/" + _anatomyInputDirectoryBox->itemText(indexInBox);
    QDir anatomyDir(path);
    anatomyDir.setNameFilters(QStringList()<<"*.nii");
    QStringList fileList = anatomyDir.entryList();

    _anatomyFileNameBox->clear();
    int indexRaw=-1;  int indexBrain=-1;  int indexAlign=-1;
    for (int jList=0; jList<fileList.size(); jList++)
    {
        _anatomyFileNameBox->addItem(fileList.at(jList));
        if ( fileList.at(jList) == "align.nii") indexAlign = jList;
        if ( fileList.at(jList) == "brain.nii") indexBrain = jList;
        if ( fileList.at(jList) == "raw.nii")   indexRaw   = jList;
    }
    if ( indexAlign >= 0 )
        _anatomyFileNameBox->setCurrentIndex(indexAlign);
    else if ( indexBrain >= 0)
        _anatomyFileNameBox->setCurrentIndex(indexBrain);
    else if ( indexRaw >= 0)
        _anatomyFileNameBox->setCurrentIndex(indexRaw);
    else if ( _anatomyFileNameBox->count() > 0 )
        _anatomyFileNameBox->setCurrentIndex(0);
    // Enable/disable:  // if anatomy file was found, it can be used for either freeSurfer or alignment
    enableAnatomyActionButtons();
}

void MainWindow::extractFreeSurferOverlays()
{
    FUNC_ENTER << _anatomyInputDirectoryBox->currentText();
    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedExtractOverlays(int, QProcess::ExitStatus)));

    QString brainName = "t1/" + _anatomyInputDirectoryBox->currentText() + "/brain.nii";
    QString atlasName = "t1/" + _anatomyInputDirectoryBox->currentText() + "/atlas.nii";
    QString comName   = "t1/" + _anatomyInputDirectoryBox->currentText() + "/align.com";
    QString outputDir = "t1/" + _anatomyInputDirectoryBox->currentText() + "/templateOverlaysFromFreeSurfer";

    QStringList arguments;
    arguments.append(brainName);
    arguments.append("-A");
    arguments.append(atlasName);
    arguments.append("-I");
    arguments.append(comName);
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());
    arguments.append("--output-file");
    arguments.append(outputDir);
    arguments.append("--preprocess");
    arguments.append("register-overlays-forward");
    arguments.append("--quit");
    qInfo() << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    FUNC_EXIT;
}

void MainWindow::finishedExtractOverlays(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: alignment";
    qInfo() << "exit code" << exitCode << "exit status" << exitStatus;
    _centralWidget->setEnabled(true);
}

void MainWindow::alignAnatomyToTemplate()
{
    FUNC_ENTER;
    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedFMAnatomyAlignment(int, QProcess::ExitStatus)));

    QString inputFileName = "t1/" + _anatomyInputDirectoryBox->currentText() + "/" + _anatomyFileNameBox->currentText();

    QStringList arguments;
    arguments.append(inputFileName);
    arguments.append("-O");
    arguments.append("alignment");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());
    arguments.append("--output-file");
    arguments.append("align");
    qInfo() << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    FUNC_EXIT;
}

void MainWindow::finishedFMAnatomyAlignment(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: alignment";
    qInfo() << "exit code" << exitCode << "exit status" << exitStatus;
    updateAnatomyFileName();
    _centralWidget->setEnabled(true);
}

bool MainWindow::anatomyFileExists(QString fileName)
{
    bool fileExists = false;
    QString fullName = "t1/" + _anatomyInputDirectoryBox->currentText() + "/" + fileName;
    FUNC_INFO << "check file" << fullName;
    QFileInfo checkFile(fullName);
    if ( checkFile.exists() && checkFile.isFile() )
         fileExists = true;
    return fileExists;
}

void MainWindow::enableAnatomyActionButtons()
{
    QString fileName = _anatomyFileNameBox->currentText();
    bool fileIsRaw = ! fileName.compare("raw.nii");
    bool fileIsBrain = ! fileName.compare("brain.nii");

    // Check for various files
    bool rawExists      = anatomyFileExists("raw.nii");
    bool brainExists    = anatomyFileExists("brain.nii");
    bool alignExists    = anatomyFileExists("align.nii");
    bool atlasExists    = anatomyFileExists("atlas.nii");
    bool alignComExists = anatomyFileExists("align.com");

    fileName = "t1/" + _anatomyInputDirectoryBox->currentText() + "/templateOverlaysFromFreeSurfer";
    QFileInfo checkOverlays(fileName);
    bool overlaysExist = checkOverlays.exists() && checkOverlays.isDir();

    QString freeDir = "free";
    QFileInfo checkDir(freeDir);
    bool freeExists = checkDir.exists() && checkDir.isDir();

    bool subjectIDEmpty = _subjectIDFreeSurfer->text().isEmpty();

    // enable
    _runFreeSurferButton->setEnabled(fileIsRaw  && rawExists);
    _alignAnatomyButton->setEnabled((fileIsBrain && brainExists) || (fileIsRaw  && rawExists));
    _extractFreeSurferOverlaysButton->setEnabled(brainExists && atlasExists && alignComExists && !subjectIDEmpty);

    FUNC_INFO << "colorize" << alignExists << freeExists << overlaysExist;

    // colorize
    if ( alignExists )
        _alignAnatomyButton->setStyleSheet("background-color:lightYellow;");
    if ( freeExists )
        _runFreeSurferButton->setStyleSheet("background-color:lightYellow;");
    if ( overlaysExist )
        _extractFreeSurferOverlaysButton->setStyleSheet("background-color:lightYellow;");
}

void MainWindow::runFreeSurfer()
{
    FUNC_ENTER;
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, "Subject name OK?", "Run freeSurfer (this takes hours)?",
                                  QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::Yes)
    {
        auto *process = new QProcess;
        _centralWidget->setEnabled(false);
        _outputBrowser->setWindowTitle("Run freeSurfer: recon-all");
        showBrowser(true);
        QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
        connect(process, SIGNAL(finished(int, QProcess::ExitStatus)),
                this, SLOT(finishedRunFreeSurfer(int, QProcess::ExitStatus)));

        QString exe = _scriptDirectory + "runFreeSurfer.csh";
        QStringList arguments;
        arguments.append(_anatomyInputDirectoryBox->currentText());
        arguments.append(_subjectIDFreeSurfer->text());
        qInfo() << exe << arguments;
        process->start(exe,arguments);
    }
}

void MainWindow::finishedRunFreeSurfer(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: freeSurfer";
    updateAnatomyFileName();
    _centralWidget->setEnabled(true);
    showBrowser(false);
}

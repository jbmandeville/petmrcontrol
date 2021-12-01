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
    auto *anatInputFileLabel  = new QLabel("file name ",_anatomyPage);
    _anatomyDirBox = new QComboBox();
    _anatomyFileBox  = new QComboBox();
    connect(_anatomyDirBox, SIGNAL(currentIndexChanged(int)),this, SLOT(changedAnatomyDirName(int)));
    connect(_anatomyFileBox, SIGNAL(activated(int)),this, SLOT(changedAnatomyFileName(int)));

    auto *displayButton = new QPushButton("display file (fastmap)");
    connect(displayButton, SIGNAL(pressed()), this, SLOT(displayAnatomy()));
    displayButton->setToolTip("Display file listed above");

    auto *displayLayout = new QVBoxLayout();
    displayLayout->addWidget(displayButton);

    auto *inputFileLayout = new QGridLayout();
    inputFileLayout->addWidget(anatDirLabel,0,0);
    inputFileLayout->addWidget(_anatomyDirBox,0,1);
    inputFileLayout->addWidget(anatInputFileLabel,1,0);
    inputFileLayout->addWidget(_anatomyFileBox,1,1);

    auto *inputLayout = new QVBoxLayout();
    inputLayout->addLayout(inputFileLayout);
    inputLayout->addLayout(displayLayout);

    auto *anatomyInputBox = new QGroupBox("Input directory for anatomy");
    anatomyInputBox->setLayout(inputLayout);
//    anatomyInputBox->setStyleSheet("border: 1px dotted gray");

    ////////////////////////////////////////////////
    // freeSurfer
    ////////////////////////////////////////////////
    auto *subjectIDLabel     = new QLabel("Subject ID",_anatomyPage);
    _subjectIDFreeSurfer     = new QLineEdit("");
    setupFreeSurferLayout->addWidget(subjectIDLabel,0,0);
    setupFreeSurferLayout->addWidget(_subjectIDFreeSurfer,0,1);

    _runFreeSurferButton       = new QPushButton("Run FreeSurfer (raw --> brain)",_anatomyPage);
    _runFreeSurferButton->setToolTip("Run freeSurfer recon-all");
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
    for (int jList=0; jList<_FastmapMSTemplateDirectories.count(); jList+=2)
        _anatomyTemplateDirectory->addItem(_FastmapMSTemplateDirectories.at(jList));
    _anatomyTemplateDirectory->setCurrentIndex(0);
    connect(_anatomyTemplateDirectory, SIGNAL(activated(int)),this, SLOT(writeSubjectVariables()));

    auto *anatomyFileLayout = new QGridLayout();
    anatomyFileLayout->addWidget(templateDirLabel,0,0);
    anatomyFileLayout->addWidget(_anatomyTemplateDirectory,0,1);
    anatomyFileLayout->setSpacing(0);

    _smoothingAnatomy = new QLineEdit("0.");
    connect(_smoothingAnatomy, SIGNAL(editingFinished()), this, SLOT(changedSmoothingAnatomy()));

    auto *anatomyAlignmentLayout = new QVBoxLayout();
    _alignAnatomyButton  = new QPushButton("Align to template (raw/brain --> align)",_anatomyPage);
    connect(_alignAnatomyButton, SIGNAL(pressed()), this, SLOT(alignAnatomyToTemplate()));
    _alignAnatomyButton->setToolTip("Align to the multi-subject template" + _anatomyTemplateDirectory->currentText());

    _extractFreeSurferOverlaysButton  = new QPushButton("extract freeSurfer ROIs from atlas file)",_anatomyPage);
    connect(_extractFreeSurferOverlaysButton, SIGNAL(pressed()), this, SLOT(extractFreeSurferOverlays()));
    _extractFreeSurferOverlaysButton->setToolTip("Move freeSurfer segmented ROIs to template space.");

    anatomyAlignmentLayout->addLayout(anatomyFileLayout);
    anatomyAlignmentLayout->addWidget(_alignAnatomyButton);
    anatomyAlignmentLayout->addWidget(_extractFreeSurferOverlaysButton);

    auto *anatomyAlignmentBox = new QGroupBox("Use fastmap to align anatomy to a multi-subject template");
    anatomyAlignmentBox->setLayout(anatomyAlignmentLayout);
//    anatomyAlignmentBox->setStyleSheet("border: 1px dotted gray");

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
        int nList = qMin(3,subDirs.count());
        QString list;  list.append("[..]/");
        for (int jList=nList; jList>0; jList--)
        {
            list.append(subDirs.at(subDirs.count()-jList));
            if ( jList != 1 ) list.append("/");
        }
        _queryDownloadGroupBox->setEnabled(false);
        setWindowTitle(QString("subject %1 @ %2").arg(_subjectIDFreeSurfer->text()).arg(list));
    }
    else
    { // freeSurfer directory does not exist: just use directory
        QDir thisDir = QDir::currentPath();
        QStringList subDirs = thisDir.absolutePath().split("/");
        int nList = qMin(3,subDirs.count());
        QString list;  list.append("[..]/");
        for (int jList=nList; jList>0; jList--)
        {
            list.append(subDirs.at(subDirs.count()-jList));
            if ( jList != 1 ) list.append("/");
        }
        setWindowTitle(QString("%1").arg(list));
    }
    enableAnatomyActionButtons();
}

void MainWindow::openedAnatomyPage()
{
    FUNC_ENTER;

    QDir const anatomyTopDir("./t1");
    if (!anatomyTopDir.exists())
    {
        _anatomyDirBox->clear();
        return;
    }
    QStringList const folderList = anatomyTopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    FUNC_INFO << folderList;
    _anatomyDirBox->clear();
    for (int jList=0; jList<folderList.size(); jList++)
        _anatomyDirBox->addItem(folderList.at(jList));
    _anatomyDirBox->setCurrentIndex(_anatomyDirBox->count()-1);

    readSmoothing(0);
    readSmoothing(1);
    readSmoothing(2);

//    changedAnatomyDirName(_anatomyDirBox->currentIndex());
}

void MainWindow::changedAnatomyDirName(int indexInBox)
{
    updateAnatomyFileName();
}

void MainWindow::updateAnatomyFileName()
{
    // For the current selection, show all NIFTI files
    QString path = "./t1/" + _anatomyDirBox->currentText();
//    QString path = "./t1/" + _anatomyDirBox->itemText(indexInBox);
    QDir anatomyDir(path);
    anatomyDir.setNameFilters(QStringList()<<"raw.nii"<<"brain.nii"<<"align.nii");
    QStringList fileList = anatomyDir.entryList();

    _anatomyFileBox->clear();
    int indexRaw=-1;  int indexBrain=-1;  int indexAlign=-1;
    for (int jList=0; jList<fileList.size(); jList++)
    {
        _anatomyFileBox->addItem(fileList.at(jList));
        if ( fileList.at(jList) == "align.nii") indexAlign = jList;
        if ( fileList.at(jList) == "brain.nii") indexBrain = jList;
        if ( fileList.at(jList) == "raw.nii")   indexRaw   = jList;
    }
    if ( indexAlign >= 0 )
        _anatomyFileBox->setCurrentIndex(indexAlign);
    else if ( indexBrain >= 0)
        _anatomyFileBox->setCurrentIndex(indexBrain);
    else if ( indexRaw >= 0)
        _anatomyFileBox->setCurrentIndex(indexRaw);
    else if ( _anatomyFileBox->count() > 0 )
        _anatomyFileBox->setCurrentIndex(0);
    // Enable/disable:  // if anatomy file was found, it can be used for either freeSurfer or alignment
    enableAnatomyActionButtons();
}

void MainWindow::extractFreeSurferOverlays()
{
    FUNC_ENTER << _anatomyDirBox->currentText();
    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedExtractOverlays(int, QProcess::ExitStatus)));

    QString brainName = "t1/" + _anatomyDirBox->currentText() + "/brain.nii";
    QString atlasName = "t1/" + _anatomyDirBox->currentText() + "/atlas.nii";
    QString comName   = "t1/" + _anatomyDirBox->currentText() + "/align.com";
    QString outputDir = "t1/" + _anatomyDirBox->currentText() + "/templateOverlaysFromFreeSurfer";

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

    QString inputFileName = "t1/" + _anatomyDirBox->currentText() + "/" + _anatomyFileBox->currentText();

    QStringList arguments;
    arguments.append(inputFileName);
    arguments.append("-O");
    arguments.append("alignment");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());
    arguments.append("--output-file");
    arguments.append("align");

    if ( anatomyFileExists("align.com") )
    {
        QString comName = "t1/" + _anatomyDirBox->currentText() + "/align.com";
        arguments.append("-I");
        arguments.append(comName);
    }


    qInfo() << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    QMessageBox msgBox;
    QString line1 = "1) Check/adjust alignment (calculator button)\n";
    QString line2 = "2) Run images through alignment pipeline (run button)\n";
    QString line3 = "3) Save images (next to run button)\n";
    QString line4 = "4) quit gracefully (control-q or menu exit)";
    QString text = line1 + line2 + line3 + line4;
    msgBox.setText(text);
    msgBox.setIcon(QMessageBox::Information);
    msgBox.exec();

    FUNC_EXIT;
}

void MainWindow::finishedFMAnatomyAlignment(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: alignment";
    qInfo() << "exit code" << exitCode << "exit status" << exitStatus;
    updateAnatomyFileName();
    _centralWidget->setEnabled(true);
}

bool MainWindow::anatomyFileExists(QString dirName, QString fileName)
{
    bool fileExists = false;
    QString fullName = "t1/" + dirName + "/" + fileName;
    QFileInfo checkFile(fullName);
    if ( checkFile.exists() && checkFile.isFile() )
         fileExists = true;
    return fileExists;
}

void MainWindow::enableAnatomyActionButtons()
{
    QString fileName = _anatomyFileBox->currentText();
    bool fileIsRaw = ! fileName.compare("raw.nii");
    bool fileIsBrain = ! fileName.compare("brain.nii");

    // Check for various files
    bool rawExists      = anatomyFileExists("raw.nii");
    bool brainExists    = anatomyFileExists("brain.nii");
    bool alignExists    = anatomyFileExists("align.nii");
    bool atlasExists    = anatomyFileExists("atlas.nii");
    bool alignComExists = anatomyFileExists("align.com");

    fileName = "t1/" + _anatomyDirBox->currentText() + "/templateOverlaysFromFreeSurfer";
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
        arguments.append(_anatomyDirBox->currentText());
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

void MainWindow::displayAnatomy()
{
    FUNC_ENTER;
    auto *process = new QProcess;

    QString rootFileName = _anatomyFileBox->currentText();
    QString inputFileName = "t1/" + _anatomyDirBox->currentText() + "/" + rootFileName;
    QStringList arguments;
    arguments.append(inputFileName);

    if ( !rootFileName.compare("align.nii") )
    {
        arguments.append("-T");
        arguments.append(_anatomyTemplateDirectory->currentText());

        QString ovlListName = "t1/" + _anatomyDirBox->currentText() + "/templateOverlaysFromFreeSurfer/overlay-list.dat";
        QFileInfo checkFile(ovlListName);
        if ( checkFile.exists() && checkFile.isFile() )
        {
            arguments.append("-o");
            arguments.append(ovlListName);
        }
    }
    else if ( !rootFileName.compare("brain.nii") )
    {
        if ( anatomyFileExists("atlas.nii") )
        {
            QString atlasFileName = "t1/" + _anatomyDirBox->currentText() + "/atlas.nii";
            arguments.append("-A");
            arguments.append(atlasFileName);
        }
    }
    qInfo() << _fastmapProcess << arguments;
    process->startDetached(_fastmapProcess,arguments);

    FUNC_EXIT;
}

#include <QtWidgets>
#include "petmrMain.h"

void MainWindow::createfMRIPage()
{
    _fmriPage = new QWidget();

    _fMRIRunItemBox = new QListWidget();
    _fMRIRunItemBox->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::MinimumExpanding);
    _fMRIRunItemBox->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    connect(_fMRIRunItemBox, SIGNAL(itemChanged(QListWidgetItem*)),this, SLOT(changedfMRIRunCheckBox(QListWidgetItem*)));
//    _fMRIRunItems = new QVector<QListWidgetItem *>;

    auto *runsLayout = new QVBoxLayout();
    runsLayout->addWidget(_fMRIRunItemBox);

    auto *templateLabel = new QLabel("template directory");
    auto *fileLabel     = new QLabel("file name(s)");
    auto *rangeLabel    = new QLabel("averaging range for template");
    _fMRITemplateDirBox = new QComboBox();
    _fMRIFileBox        = new QComboBox();
    _fMRIMCRange        = new QLineEdit("1-10");
    connect(_fMRITemplateDirBox, SIGNAL(activated(int)),this, SLOT(changefMRITemplateDirectory(int)));
    connect(_fMRIFileBox,        SIGNAL(activated(int)),this, SLOT(changedfMRIFileName(int)));

    auto *fileLayout = new QGridLayout();
    fileLayout->addWidget(templateLabel,0,0);
    fileLayout->addWidget(_fMRITemplateDirBox,0,1);
    fileLayout->addWidget(rangeLabel,1,0);
    fileLayout->addWidget(_fMRIMCRange,1,1);
    fileLayout->addWidget(fileLabel,2,0);
    fileLayout->addWidget(_fMRIFileBox,2,1);

    auto *displayButton = new QPushButton("display file(s) with fastmap");
    connect(displayButton, SIGNAL(pressed()), this, SLOT(displayEPI()));
    auto *displayLayout = new QVBoxLayout();
    displayLayout->addWidget(displayButton);

    auto *setupLayout = new QVBoxLayout();
    setupLayout->addLayout(runsLayout);
    setupLayout->addLayout(fileLayout);
    setupLayout->addLayout(displayLayout);

    auto *runsBox = new QGroupBox("List of EPI runs to include in analysis");
    runsBox->setLayout(setupLayout);
//    runsBox->setStyleSheet("border: 1px dotted gray");

    _doEverythingEPIButton  = new QPushButton("do everything",_fmriPage);
    _doEverythingEPIButton->setEnabled(false);
    _doEverythingEPIButton->setCheckable(true);

    _resliceEPIButton       = new QPushButton("reslice runs as necessary (raw --> reslice)",_fmriPage);
    _motionCorrectEPIButton = new QPushButton("motion-correct runs (raw/reslice -->mc)",_fmriPage);
    _alignEPIButton         = new QPushButton("Align to template (raw/reslice/mc -> align)",_fmriPage);
    connect(_doEverythingEPIButton,  SIGNAL(clicked(bool)), this, SLOT(doEverthingEPI()));
    connect(_resliceEPIButton,       SIGNAL(pressed()), this, SLOT(resliceEPI()));
    connect(_motionCorrectEPIButton, SIGNAL(pressed()), this, SLOT(motionCorrectEPI()));
    connect(_alignEPIButton,         SIGNAL(pressed()), this, SLOT(alignEPI()));

    auto *actionLayout = new QVBoxLayout();
    actionLayout->addWidget(_doEverythingEPIButton);
    actionLayout->addWidget(_resliceEPIButton);
    actionLayout->addWidget(_motionCorrectEPIButton);
    actionLayout->addWidget(_alignEPIButton);

    auto *fMRIActionBox = new QGroupBox("Process EPI runs");
    fMRIActionBox->setLayout(actionLayout);
    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(runsBox);
    pageLayout->addWidget(fMRIActionBox);
    _fmriPage->setLayout(pageLayout);
}

void MainWindow::doEverthingEPI()
{
    FUNC_ENTER;
    _resliceEPIButton->setStyleSheet("background-color:white;");
    _motionCorrectEPIButton->setStyleSheet("background-color:white;");
    _alignEPIButton->setStyleSheet("background-color:white;");
    _doEverythingEPIButton->setStyleSheet("background-color:white;");

    FUNC_INFO << _resliceEPIButton->isEnabled();
    if ( _resliceEPIButton->isEnabled() )
        resliceEPI();
    else
        motionCorrectEPI();
}

void MainWindow::changedfMRIRunCheckBox(QListWidgetItem *item)
{
    int iSelected=-1;
    for ( int jItem=0; jItem<_fMRIRunItems.size(); jItem++)
    {
        if ( item == &_fMRIRunItems.at(jItem) )
            iSelected = jItem;
    }
    FUNC_INFO << iSelected;
    if ( !_fMRIRunItems[iSelected].checkState() &&
         iSelected == _fMRITemplateDirBox->currentIndex() )
    { // unselected the template directory, so pick the first selection
        iSelected = -1;  int numberChecked=0;
        for ( int jItem=0; jItem<_fMRIRunItems.size(); jItem++)
        {
            if ( _fMRIRunItems[jItem].checkState() )
            {
                numberChecked++;
                if ( iSelected < 0 ) iSelected = jItem;
            }
        }
        FUNC_INFO << "numberChecked" << numberChecked;
        iSelected = qMax(iSelected,0);
        _fMRITemplateDirBox->setCurrentIndex(iSelected);
        changefMRITemplateDirectory(iSelected);
    }
    enableEPIActionButtons();
}

void MainWindow::changefMRITemplateDirectory(int indexInBox)
{// fMRI template directory = "004", "005", ...
    FUNC_ENTER << indexInBox << _fMRIFiles.size();
    _dimEPITemplate = _fMRIFiles[indexInBox].dim;
    FUNC_INFO << "_dimEPITemplate3" << _dimEPITemplate.x << _dimEPITemplate.y << _dimEPITemplate.z;
    for (int jList=0; jList<_fMRIRunItems.size(); jList++)
    {
        if ( jList == indexInBox && _fMRIRunItems[jList].checkState() )
            _fMRIRunItems[jList].setBackgroundColor(Qt::cyan);
        else
            _fMRIRunItems[jList].setBackgroundColor(Qt::white);
    }
}

void MainWindow::openedfMRIPage()
{
    FUNC_ENTER;
    openedAnatomyPage();

    qInfo() << "query EPI data";
    QDir const fMRITopDir("./epi");
    if (!fMRITopDir.exists())
    {
        _fMRIRunItems.clear();
        return;
    }
    FUNC_INFO << 1;
    if (_fMRITemplateDirBox->count() == 0 )
    {
        QStringList const folderList = fMRITopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
        for (int jList=0; jList<folderList.size(); jList++)
        {
            if ( epiFileExists(folderList.at(jList),"raw.nii")     ||
                 epiFileExists(folderList.at(jList),"reslice.nii") ||
                 epiFileExists(folderList.at(jList),"align.nii")  )
                _fMRITemplateDirBox->addItem(folderList.at(jList));
        }
    }

    // Find the template directory
    int iTemplate=-1;
    for (int jList=0; jList<_fMRITemplateDirBox->count(); jList++)
    {
        if (epiFileExists(_fMRITemplateDirBox->itemText(jList),"mcTemplate.nii"))
        {
            iTemplate = jList;
            FUNC_INFO << "found template" << iTemplate;
        }
    }
    iTemplate = qMax(0,iTemplate);

    // Given the fMRI template directory, pick a file type
    // For the first selection, show all NIFTI files
    _fMRITemplateDirBox->setCurrentIndex(iTemplate);
    changefMRITemplateDirectory(_fMRITemplateDirBox->currentIndex());
    populateEPIFileNameBox();

    enableEPIActionButtons();

    FUNC_EXIT;
}

void MainWindow::updateEPIFileNameBox(QString fileName)
{
    FUNC_ENTER << fileName;
    populateEPIFileNameBox();
    int foundIndex=-1;
    for (int jList=0; jList<_fMRIFileBox->count(); jList++)
    {
        if ( !_fMRIFileBox->itemText(jList).compare(fileName) )
            foundIndex = jList;
    }
    if ( foundIndex >= 0 ) changedfMRIFileName(foundIndex);

    FUNC_EXIT << "foundIndex" << foundIndex;
}

void MainWindow::populateEPIFileNameBox()
{
    QString path = "./epi/" + _fMRITemplateDirBox->currentText();
    FUNC_ENTER << path;
    QDir templateDir(path);
    templateDir.setNameFilters(QStringList()<<"raw.nii"<<"reslice.nii"<<"mc.nii"<<"align.nii");
    QStringList fileList = templateDir.entryList();
    FUNC_INFO << fileList;
    _fMRIFileBox->clear();
    for (int jList=0; jList<fileList.size(); jList++)
        _fMRIFileBox->addItem(fileList.at(jList));

    setDefaultIndexEPIFileNameBox();
}

void MainWindow::setDefaultIndexEPIFileNameBox()
{
    int indexRaw=-1;  int indexReslice=-1;  int indexMC=-1;  int indexAlign=-1;
    for (int jList=0; jList<_fMRIFileBox->count(); jList++)
    {
        QString fileName = _fMRIFileBox->itemText(jList);
        if ( fileName == "raw.nii")     indexRaw     = jList;
        if ( fileName == "reslice.nii") indexReslice = jList;
        if ( fileName == "mc.nii")      indexMC      = jList;
        if ( fileName == "align.nii")   indexAlign   = jList;
    }
    if ( indexAlign >= 0 )
        _fMRIFileBox->setCurrentIndex(indexAlign);
    else if ( indexMC >= 0 )
        _fMRIFileBox->setCurrentIndex(indexMC);
    else if ( indexReslice >= 0 )
        _fMRIFileBox->setCurrentIndex(indexReslice);
    else if ( indexRaw >= 0 )
        _fMRIFileBox->setCurrentIndex(indexRaw);
    else
        _fMRIFileBox->setCurrentIndex(0);
    changedfMRIFileName(_fMRIFileBox->currentIndex());
}

void MainWindow::changedfMRIFileName(int indexInBox)
{ // fMRI file name = "raw.nii", "reslice.nii", "mc.nii", "align.nii"
    FUNC_ENTER << indexInBox;
    _fMRIFileBox->setCurrentIndex(indexInBox);
    QDir const fMRITopDir("./epi/");
    if (!fMRITopDir.exists())
    {
        _fMRIRunItems.clear();
        return;
    }
    QStringList const folderList = fMRITopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);

    QString path = fMRITopDir.dirName() + "/" + _fMRITemplateDirBox->currentText();
    QString fileName = path + "/" + _fMRIFileBox->currentText();
    getDimensions(fileName, _dimEPITemplate);
    FUNC_INFO << "_dimEPITemplate1" << _dimEPITemplate.x << _dimEPITemplate.y << _dimEPITemplate.z;
    FUNC_INFO << "from file" << fileName;

    // Count files for allocation
    int nFiles=0;
    for (int jList=0; jList<folderList.size(); jList++)
    {
        FourDFile file;
        file.name        = "epi/" + folderList.at(jList) + "/" + _fMRIFileBox->currentText();
        QString tagsName = "epi/" + folderList.at(jList) + "/time-tags.txt";
        getDimensions(file.name, file.dim);
        getTimeTags(tagsName,file.timeTags,file.timeText);
        if ( file.timeTags.size() != 0 )
            nFiles++;
    }

    FUNC_INFO << "resize to" << nFiles;
    _fMRIRunItems.resize(nFiles);
    _fMRIFiles.resize(nFiles);
    nFiles=0;
    for (int jList=0; jList<folderList.size(); jList++)
    {
        FourDFile file;
        file.name        = "epi/" + folderList.at(jList) + "/" + _fMRIFileBox->currentText();
        QString tagsName = "epi/" + folderList.at(jList) + "/time-tags.txt";
        QString text = getDimensions(file.name, file.dim);
        getTimeTags(tagsName,file.timeTags,file.timeText);
        FUNC_INFO << "jList" << jList << "file" << file.name;
        if ( file.timeTags.size() != 0 )
        {
            FUNC_INFO << "add fmri file" << folderList.at(jList) << text;
            QListWidgetItem item;

            item.setText(folderList.at(jList) + " :    " + text);
            bool hide = file.dim.z == 0;
            item.setHidden(hide);
            if ( hide )
            {
                item.setCheckState(Qt::Unchecked);
                item.setFlags(Qt::ItemIsUserCheckable);
            }
            else
            {
                item.setCheckState(Qt::Checked);
                item.setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
            }
            _fMRIRunItems[nFiles] = item;
            _fMRIFiles[nFiles] = file;
            _fMRIRunItemBox->addItem(&_fMRIRunItems[nFiles]);
            nFiles++;
            FUNC_INFO << "added fmri file" << folderList.at(jList) << text;
        }
    }

    FUNC_INFO << "_fMRIFiles size" << _fMRIFiles.size();
    FUNC_INFO << "_dimEPITemplate2" << _dimEPITemplate.x << _dimEPITemplate.y << _dimEPITemplate.z;
    changefMRITemplateDirectory(_fMRITemplateDirBox->currentIndex());

    enableEPIActionButtons();
    FUNC_EXIT;
}

void MainWindow::enableEPIActionButtons()
{
    FUNC_ENTER;
    bool allSameDimension = true;
    for (int jList=0; jList<_fMRIFiles.size(); jList++)
    {
        bool includeFile = _fMRIRunItems[jList].checkState();
        if ( includeFile )        {
            FourDFile file = _fMRIFiles[jList];
            allSameDimension &= file.dim.x == _dimEPITemplate.x;
            allSameDimension &= file.dim.y == _dimEPITemplate.y;
            allSameDimension &= file.dim.z == _dimEPITemplate.z;
        }
    }
    _resliceEPIButton->setEnabled(!allSameDimension);
    bool selectRaw     = !_fMRIFileBox->currentText().compare("raw.nii");
    bool selectReslice = !_fMRIFileBox->currentText().compare("reslice.nii");
    bool selectMC      = !_fMRIFileBox->currentText().compare("mc.nii");
    _motionCorrectEPIButton->setEnabled(allSameDimension && (selectRaw || selectReslice));
    _alignEPIButton->setEnabled(allSameDimension && (selectRaw || selectReslice || selectMC));
    _doEverythingEPIButton->setEnabled(selectRaw);

    if ( epiFileExists("reslice.nii") )
        _resliceEPIButton->setStyleSheet("background-color:lightYellow;");
    if ( epiFileExists("mc.nii") )
        _motionCorrectEPIButton->setStyleSheet("background-color:lightYellow;");
    if ( epiFileExists("align.nii") )
    {
        _alignEPIButton->setStyleSheet("background-color:lightYellow;");
        _doEverythingEPIButton->setStyleSheet("background-color:lightYellow;");
    }


    FUNC_EXIT << allSameDimension;
}

void MainWindow::resliceEPI()
{
    FUNC_ENTER;
    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)),
            this, SLOT(finishedFMResliceEPI(int, QProcess::ExitStatus)));

    QStringList arguments;
    arguments.append("-O");
    arguments.append("reslice");
    arguments.append("-a");
    QString templateFileName = "./epi/" + _fMRITemplateDirBox->currentText()
            + "/" + _fMRIFileBox->currentText();
    arguments.append(templateFileName);
    arguments.append("--preprocess");
    arguments.append("reslice");
    arguments.append("--output-file");
    arguments.append("reslice");
    arguments.append("--quit");
    for (int jFile=0; jFile<_fMRIFiles.size(); jFile++)
    {
        bool includeFile = _fMRIRunItems[jFile].checkState();
        if ( includeFile )
        {
            FUNC_INFO << "jFile" << jFile << "name" << _fMRIFiles[jFile].name << "args" << arguments;
            iPoint4D dim = _fMRIFiles[jFile].dim;
            bool sameDimensions = dim.x == _dimEPITemplate.x;
            sameDimensions     &= dim.y == _dimEPITemplate.y;
            sameDimensions     &= dim.z == _dimEPITemplate.z;
            if ( !sameDimensions )
                arguments.append(_fMRIFiles[jFile].name);
        }
    }
    FUNC_INFO << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    FUNC_EXIT;
}
void MainWindow::finishedFMResliceEPI(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: reslice EPI";
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;

    QString exe = _scriptDirectory + "linkNonReslicedRaws.csh";
    QStringList arguments;
    for (int jFile=0; jFile<_fMRIFiles.size(); jFile++)
    {
        bool includeFile = _fMRIRunItems[jFile].checkState();
        if ( includeFile )
        {
            FUNC_INFO << "jFile" << jFile << "name" << _fMRIFiles[jFile].name << "args" << arguments;
            iPoint4D dim = _fMRIFiles[jFile].dim;
            bool sameDimensions = dim.x == _dimEPITemplate.x;
            sameDimensions     &= dim.y == _dimEPITemplate.y;
            sameDimensions     &= dim.z == _dimEPITemplate.z;
            if ( sameDimensions )
                arguments.append(_fMRIFiles[jFile].name);
        }
    }

    qInfo() <<  exe << arguments;
    auto *process = new QProcess;
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)),
            this, SLOT(finishedLinkResliceEPI(int, QProcess::ExitStatus)));

    process->start(exe,arguments);
}

void MainWindow::finishedLinkResliceEPI(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: link EPI";
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;

    _centralWidget->setEnabled(true);
    updateEPIFileNameBox("reslice.nii");

    FUNC_INFO << "next up?" << _doEverythingEPIButton->isChecked();
    if ( _doEverythingEPIButton->isChecked() )
    {
        _resliceEPIButton->setStyleSheet("background-color:lightYellow;");
        _motionCorrectEPIButton->setStyleSheet("background-color:white;");
        _alignEPIButton->setStyleSheet("background-color:white;");
        enableEPIActionButtons();
        motionCorrectEPI();
    }
    else
        enableEPIActionButtons();
}

bool MainWindow::epiFileExists(QString dirName, QString fileName)
{
    bool fileExists = false;
    QString fullName = "epi/" + dirName + "/" + fileName;
    FUNC_INFO << "check file" << fullName;
    QFileInfo checkFile(fullName);
    if ( checkFile.exists() && checkFile.isFile() )
         fileExists = true;
    return fileExists;
}

void MainWindow::alignEPI()
{
    FUNC_ENTER;
    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)),
            this, SLOT(finishedFMAlignEPI(int, QProcess::ExitStatus)));

    bool alignEPIComExists = epiFileExists("align.com");
    bool alignT1ComExists = anatomyFileExists("align.com");

    QStringList arguments;
    arguments.append("-O");
    arguments.append("align");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());
    arguments.append("--output-file");
    arguments.append("align");
    // include individual anatomy file
    QString individualAnatomy = "t1/" + _anatomyDirBox->currentText() + "/align.nii";
    QFileInfo checkFile(individualAnatomy);
    if ( checkFile.exists() && checkFile.isFile() )
    {
        arguments.append("-a");
        arguments.append(individualAnatomy);
    }
    for (int jFile=0; jFile<_fMRIFiles.size(); jFile++)
    {
        bool includeFile = _fMRIRunItems[jFile].checkState();
        if ( includeFile )
            arguments.append(_fMRIFiles[jFile].name);
    }
    if ( alignEPIComExists )
    {
        QString comName = "epi/" + _fMRITemplateDirBox->currentText() + "/align.com";
        arguments.append("-I");
        arguments.append(comName);
    }
    else if ( alignT1ComExists )
    {
        QString comName = "t1/" + _anatomyDirBox->currentText() + "/align.com";
        arguments.append("-I");
        arguments.append(comName);
    }
    FUNC_INFO << _fastmapProcess << arguments;
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
void MainWindow::finishedFMAlignEPI(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: align EPI";
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    updateEPIFileNameBox("align.nii");
    _centralWidget->setEnabled(true);
    if ( _doEverythingEPIButton->isChecked() )
    {
        _doEverythingEPIButton->setStyleSheet("background-color:lightYellow;");
        _doEverythingEPIButton->setChecked(false);
    }
}

void MainWindow::motionCorrectEPI()
{
    auto *process = new QProcess;
    _outputBrowser->setWindowTitle("Motion-correct EPI");
    showBrowser(true);
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedMotionCorrectEPI(int, QProcess::ExitStatus)));

    process->setProcessChannelMode(QProcess::MergedChannels);

    QString exe = _scriptDirectory + "motionCorrectEPI.csh";
    QStringList arguments;
    // arguments: [MC template dir] [template range (e.g. 1-10)] [list of files: "epi/011/reslice.nii ..."]
    arguments.append(_fMRITemplateDirBox->currentText());
    arguments.append(_fMRIMCRange->text());
    // list of EPI directories to motion-correct ...
    for (int jFile=0; jFile<_fMRIFiles.size(); jFile++)
    {
        bool includeFile = _fMRIRunItems[jFile].checkState();
        if ( includeFile )
            arguments.append(_fMRIFiles[jFile].name);
    }
    qInfo() <<  exe << arguments;
    process->start(exe,arguments);
    _centralWidget->setEnabled(false);
}
void MainWindow::finishedMotionCorrectEPI(int exitCode, QProcess::ExitStatus exitStatus)
{
    qInfo() << "finished: motion-correct EPI";
    updateEPIFileNameBox("mc.nii");
    _centralWidget->setEnabled(true);
    showBrowser(false);
    if ( _doEverythingEPIButton->isChecked() )
    {
        _motionCorrectEPIButton->setStyleSheet("background-color:lightYellow;");
        _alignEPIButton->setStyleSheet("background-color:white;");
        alignEPI();
    }
    else
        enableEPIActionButtons();
}

void MainWindow::displayEPI()
{
    FUNC_ENTER;
    auto *process = new QProcess;

    QString rootFileName = _fMRIFileBox->currentText();
    QStringList arguments;
    for (int jFile=0; jFile<_fMRIFiles.size(); jFile++)
    {
        bool includeFile = _fMRIRunItems[jFile].checkState();
        if ( includeFile )
            arguments.append(_fMRIFiles[jFile].name);
    }

    if ( !rootFileName.compare("align.nii") )
    {
        QString ovlListName = "t1/" + _anatomyDirBox->currentText() + "/templateOverlaysFromFreeSurfer/overlay-list.dat";
        FUNC_INFO << "check file" << ovlListName;
        QFileInfo checkFile(ovlListName);
        if ( checkFile.exists() && checkFile.isFile() )
        {
            arguments.append("-o");
            arguments.append(ovlListName);
        }
        else
            FUNC_INFO << "does not exist:" << ovlListName;
    }

    qInfo() << _fastmapProcess << arguments;
    process->startDetached(_fastmapProcess,arguments);

    FUNC_EXIT;
}

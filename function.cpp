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
    _fMRIFileNameBox          = new QComboBox();
    _fMRIMCRange              = new QLineEdit("1-10");
    connect(_fMRITemplateDirBox, SIGNAL(activated(int)),this, SLOT(changefMRITemplateDirectory(int)));
    connect(_fMRIFileNameBox,          SIGNAL(activated(int)),this, SLOT(changedfMRIFileName(int)));

    auto *fileLayout = new QGridLayout();
    fileLayout->addWidget(templateLabel,0,0);
    fileLayout->addWidget(_fMRITemplateDirBox,0,1);
    fileLayout->addWidget(rangeLabel,1,0);
    fileLayout->addWidget(_fMRIMCRange,1,1);
    fileLayout->addWidget(fileLabel,2,0);
    fileLayout->addWidget(_fMRIFileNameBox,2,1);

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

    _resliceEPIButton       = new QPushButton("reslice runs as necessary (raw --> reslice)",_fmriPage);
    _motionCorrectEPIButton = new QPushButton("motion-correct runs (raw/reslice -->mc)",_fmriPage);
    _alignEPIButton         = new QPushButton("Align to template (raw/reslice/mc -> align)",_fmriPage);
    connect(_resliceEPIButton,       SIGNAL(pressed()), this, SLOT(resliceEPI()));
    connect(_motionCorrectEPIButton, SIGNAL(pressed()), this, SLOT(motionCorrectEPI()));
    connect(_alignEPIButton,         SIGNAL(pressed()), this, SLOT(alignEPI()));

    auto *actionLayout = new QVBoxLayout();
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
        if (epiFileExists("mcTemplate.nii"))
        {
            iTemplate = jList;
            FUNC_INFO << "found template" << iTemplate;
        }
    }
    iTemplate = qMax(0,iTemplate);

    // Given the fMRI template directory, pick a file type
    // For the first selection, show all NIFTI files
    _fMRITemplateDirBox->setCurrentIndex(iTemplate);
    updateFileNameBox();
    changefMRITemplateDirectory(_fMRITemplateDirBox->currentIndex());

    enableEPIActionButtons();

    FUNC_EXIT;
}

void MainWindow::updateFileNameBox()
{
    QString path = "./epi/" + _fMRITemplateDirBox->currentText();
    FUNC_ENTER << path;
    QDir templateDir(path);
    templateDir.setNameFilters(QStringList()<<"*.nii");
    QStringList fileList = templateDir.entryList();
    FUNC_INFO << fileList;
    _fMRIFileNameBox->clear();
    int indexRaw=-1;  int indexReslice=-1;  int indexMC=-1;  int indexAlign=-1;
    for (int jList=0; jList<fileList.size(); jList++)
    {
        if ( fileList.at(jList).compare("mcTemplate.nii") )  // don't add template to list
            _fMRIFileNameBox->addItem(fileList.at(jList));
        if ( fileList.at(jList) == "raw.nii")     indexRaw     = jList;
        if ( fileList.at(jList) == "reslice.nii") indexReslice = jList;
        if ( fileList.at(jList) == "mc.nii")      indexMC      = jList;
        if ( fileList.at(jList) == "align.nii")   indexAlign   = jList;
    }
    if ( indexAlign >= 0 )
        _fMRIFileNameBox->setCurrentIndex(indexAlign);
    else if ( indexMC >= 0 )
        _fMRIFileNameBox->setCurrentIndex(indexMC);
    else if ( indexReslice >= 0 )
        _fMRIFileNameBox->setCurrentIndex(indexReslice);
    else if ( indexRaw >= 0 )
        _fMRIFileNameBox->setCurrentIndex(indexRaw);
    else
        _fMRIFileNameBox->setCurrentIndex(0);

    FUNC_INFO << "indexRaw indexReslice indexMC indexAlign" << indexRaw << indexReslice << indexMC << indexAlign;

    if ( indexReslice >= 0 )
        _resliceEPIButton->setStyleSheet("background-color:lightYellow;");
    if ( indexMC >= 0 )
        _motionCorrectEPIButton->setStyleSheet("background-color:lightYellow;");
    if ( indexAlign >= 0 )
        _alignEPIButton->setStyleSheet("background-color:lightYellow;");

    FUNC_INFO << 2;
    changedfMRIFileName(_fMRIFileNameBox->currentIndex());
}

void MainWindow::changedfMRIFileName(int indexInBox)
{ // fMRI file name = "raw.nii", "reslice.nii", "mc.nii", "align.nii"
    FUNC_ENTER;
    QDir const fMRITopDir("./epi/");
    if (!fMRITopDir.exists())
    {
        _fMRIRunItems.clear();
        return;
    }
    FUNC_INFO << 1;
    QStringList const folderList = fMRITopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);

    QString path = fMRITopDir.dirName() + "/" + _fMRITemplateDirBox->currentText();
    QString fileName = path + "/" + _fMRIFileNameBox->currentText();
    getDimensions(fileName, _dimEPITemplate);
    FUNC_INFO << "_dimEPITemplate1" << _dimEPITemplate.x << _dimEPITemplate.y << _dimEPITemplate.z;
    FUNC_INFO << "from file" << fileName;

    // Count files for allocation
    int nFiles=0;
    for (int jList=0; jList<folderList.size(); jList++)
    {
        FourDFile file;
        file.name        = "epi/" + folderList.at(jList) + "/" + _fMRIFileNameBox->currentText();
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
        file.name        = "epi/" + folderList.at(jList) + "/" + _fMRIFileNameBox->currentText();
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
            FUNC_INFO << "here 1";
            _fMRIRunItems[nFiles] = item;
            FUNC_INFO << "here 2";
            _fMRIFiles[nFiles] = file;
            FUNC_INFO << "here 3";
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
    bool selectRaw     = !_fMRIFileNameBox->currentText().compare("raw.nii");
    bool selectReslice = !_fMRIFileNameBox->currentText().compare("reslice.nii");
    bool selectMC      = !_fMRIFileNameBox->currentText().compare("mc.nii");
    _motionCorrectEPIButton->setEnabled(allSameDimension && (selectRaw || selectReslice));
    _alignEPIButton->setEnabled(allSameDimension && (selectRaw || selectReslice || selectMC));

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
            + "/" + _fMRIFileNameBox->currentText();
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
    process->startDetached(exe,arguments);
    _centralWidget->setEnabled(true);
    updateFileNameBox();

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

    bool alignComExists = epiFileExists("align.com");

    QStringList arguments;
    arguments.append("-O");
    arguments.append("align");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());
    // include individual anatomy file
    QString individualAnatomy = "t1/" + _anatomyDirBox->currentText() + "/align.nii";
    QFileInfo checkFile(individualAnatomy);
    if ( checkFile.exists() && checkFile.isFile() )
    {
        arguments.append("-a");
        arguments.append(individualAnatomy);
    }
    arguments.append("--output-file");
    arguments.append("align");
    for (int jFile=0; jFile<_fMRIFiles.size(); jFile++)
    {
        bool includeFile = _fMRIRunItems[jFile].checkState();
        if ( includeFile )
            arguments.append(_fMRIFiles[jFile].name);
    }
    if ( alignComExists )
    {
        QString comName = "epi/" + _fMRITemplateDirBox->currentText() + "/align.com";
        arguments.append("-I");
        arguments.append(comName);
    }
    FUNC_INFO << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    FUNC_EXIT;
}
void MainWindow::finishedFMAlignEPI(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: align EPI";
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    updateFileNameBox();
    _centralWidget->setEnabled(true);
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
    updateFileNameBox();
    _centralWidget->setEnabled(true);
    showBrowser(false);
}

void MainWindow::displayEPI()
{
    FUNC_ENTER;
    auto *process = new QProcess;

    QString rootFileName = _fMRIFileNameBox->currentText();
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

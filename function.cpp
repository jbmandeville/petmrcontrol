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

    auto *dirLayout = new QVBoxLayout();
    dirLayout->addWidget(_fMRIRunItemBox);

    auto *fileLabel     = new QLabel("file name(s)");
    _fMRIFileBox        = new QComboBox();
    connect(_fMRIFileBox,  SIGNAL(activated(int)),this, SLOT(changedfMRIFileName(int)));

    auto *fileLayout = new QGridLayout();
    fileLayout->addWidget(fileLabel,0,0);
    fileLayout->addWidget(_fMRIFileBox,0,1);

    auto *displayButton = new QPushButton("display file(s) with fastmap");
    connect(displayButton, SIGNAL(pressed()), this, SLOT(displayEPI()));

    auto *runsLayout = new QVBoxLayout();
    runsLayout->addLayout(dirLayout);
    runsLayout->addLayout(fileLayout);
    runsLayout->addWidget(displayButton);
    auto *runsBox = new QGroupBox("List of EPI runs to include in analysis");
    runsBox->setLayout(runsLayout);

    auto *templateLabel = new QLabel("directory");
    auto *rangeLabel    = new QLabel("averaging range");
    _fMRITemplateDirBox = new QComboBox();
    connect(_fMRITemplateDirBox, SIGNAL(activated(int)),this, SLOT(changefMRITemplateDirectory(int)));
    _fMRIMCRange        = new QLineEdit("1-10");
    _fMRITemplateDirBox->setMaximumWidth(150);
    _fMRIMCRange->setMaximumWidth(150);

    auto *mcLayout = new QGridLayout();
    mcLayout->addWidget(templateLabel,0,0);
    mcLayout->addWidget(_fMRITemplateDirBox,0,1);
    mcLayout->addWidget(rangeLabel,0,2);
    mcLayout->addWidget(_fMRIMCRange,0,3);
    auto *mcBox = new QGroupBox("Motion-correction template");
    mcBox->setLayout(mcLayout);

    auto *smoothingLabel = new QLabel("post-alignment smoothing width");
    _smoothingfMRI = new QLineEdit("0.");
    _smoothingfMRI->setMaximumWidth(150);
    connect(_smoothingfMRI, SIGNAL(editingFinished()), this, SLOT(changedSmoothingfMRI()));
    auto *alignLayout = new QGridLayout();
    alignLayout->addWidget(smoothingLabel,0,0);
    alignLayout->addWidget(_smoothingfMRI,0,1);
    auto *alignBox = new QGroupBox("Alignment to template space");
    alignBox->setLayout(alignLayout);

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

    auto *actionBox = new QGroupBox("Process EPI runs");
    actionBox->setLayout(actionLayout);
    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(runsBox);
    pageLayout->addWidget(mcBox);
    pageLayout->addWidget(alignBox);
    pageLayout->addWidget(actionBox);
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
    FUNC_ENTER << _fMRIRunItems.size();
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

    for (int jFile=0; jFile<_fMRIRunItems.size(); jFile++)
    {
        FUNC_INFO << "check state" << jFile << _fMRIRunItems[jFile].checkState();
    }
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

    QString range = readJipCommandFileForMCAveraging();
    FUNC_INFO << "range =" << range;
    _fMRIMCRange->setText(range);

    readSmoothing(0);
    readSmoothing(1);
    readSmoothing(2);

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
    // Count files for allocation
    int nFiles=0;
    for (int jList=0; jList<folderList.size(); jList++)
    {
        FourDFile file;
        QString tagsName = "epi/" + folderList.at(jList) + "/time-tags.txt";
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
        getTimeTags(tagsName,file.timeTags,file.timeText);
        QString text = getDimensions(file.name, file.dim);
        FUNC_INFO << "jList" << jList << "file" << file.name;
        if ( !folderList.at(jList).compare(_fMRITemplateDirBox->currentText()) )
        {
            _dimEPITemplate.x = file.dim.x;
            _dimEPITemplate.y = file.dim.y;
            _dimEPITemplate.z = file.dim.z;
        }
        if ( file.timeTags.size() != 0 )
        {
            _fMRIFiles[nFiles] = file;
            FUNC_INFO << "pet tags size" << _petFile.timeTags.size();
            if ( _petFile.timeTags.size() > 0 )
            {
                double timeInPet = (file.timeTags.at(0) - _petFile.timeTags.at(0))/60.;
                QString number;  number.setNum(timeInPet,'g',4);
                text = text + ", time in PET = " + number + " min";
            }

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
    auto *process = new QProcess();
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
    QString message = "Reslice EPI; this can take tens of minutes";
    spawnProcess(process,_fastmapProcess,arguments,message,"");

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

    auto *process = new QProcess();
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)),
            this, SLOT(finishedLinkResliceEPI(int, QProcess::ExitStatus)));

    QString message = "Link non-resliced raw.nii files to `reslice.nii' for consistency.";
    spawnProcess(process,exe,arguments,message,"");
}

void MainWindow::finishedLinkResliceEPI(int exitCode, QProcess::ExitStatus exitStatus )
{
    qInfo() << "finished: link EPI";
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;

    finishedProcess();
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
    auto *process = new QProcess();
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
    FUNC_INFO << "set args" << _fMRIRunItems.size();
    for (int jFile=0; jFile<_fMRIRunItems.size(); jFile++)
    {
        FUNC_INFO << "check state" << jFile << _fMRIRunItems[jFile].checkState();
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

    QString message = "Align EPI to template space; this requires interaction (potential tweaking)";
    spawnProcess(process,_fastmapProcess,arguments,message,"");

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
    finishedProcess();
    if ( _doEverythingEPIButton->isChecked() )
    {
        _doEverythingEPIButton->setStyleSheet("background-color:lightYellow;");
        _doEverythingEPIButton->setChecked(false);
    }
}

void MainWindow::writeJipCommandFileForMCAveraging()
{
    FUNC_ENTER;
    QString fileName = "epi/" + _fMRITemplateDirBox->currentText() + "/makeMCTemplate.com";
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);

    out << "# create template file for motion-correction" << "\n";
    QString inputFileName  = "epi/" + _fMRITemplateDirBox->currentText() + "/" + _fMRIFileBox->currentText();
    QString outputFileName = "epi/" + _fMRITemplateDirBox->currentText() + "/mcTemplate.nii";
    out << "read " << inputFileName << ">" << _fMRIMCRange->text() << " all\n";
    out << "average av\n";
    out << "write " << outputFileName << " av\n";
    out << "bye\n";
    file.close();
}

QString MainWindow::readJipCommandFileForMCAveraging()
{
    QString range = "1-10";  // default

    QString fileName = "epi/" + _fMRITemplateDirBox->currentText() + "/makeMCTemplate.com";
    QString argument = readFileTextArgument(fileName, "read");

    FUNC_INFO << "read subject" << argument;
    QStringList subList = argument.split(QRegularExpression(">"));
    FUNC_INFO << "subList" << subList;
    if ( subList.size() > 1 ) range = subList.at(1);
    return range;
}

void MainWindow::motionCorrectEPI()
{
    writeJipCommandFileForMCAveraging();

    auto *process = new QProcess();
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedMotionCorrectEPI(int, QProcess::ExitStatus)));

    process->setProcessChannelMode(QProcess::MergedChannels);

    QString exe = _scriptDirectory + "motionCorrectEPI.csh";
    QStringList arguments;
    // arguments: [MC template dir] [template range (e.g. 1-10)] [list of files: "epi/011/reslice.nii ..."]
    arguments.append(_fMRITemplateDirBox->currentText());
    QString comFileName = "epi/" + _fMRITemplateDirBox->currentText() + "/makeMCTemplate.com";
    arguments.append(comFileName);
    // list of EPI directories to motion-correct ...
    for (int jFile=0; jFile<_fMRIFiles.size(); jFile++)
    {
        bool includeFile = _fMRIRunItems[jFile].checkState();
        if ( includeFile )
            arguments.append(_fMRIFiles[jFile].name);
    }
    QString message = "Motion-correct EPI; this takes about 10 minutes";
    spawnProcess(process,exe,arguments,message,"Motion-correct EPI");
}
void MainWindow::finishedMotionCorrectEPI(int exitCode, QProcess::ExitStatus exitStatus)
{
    qInfo() << "finished: motion-correct EPI";
    updateEPIFileNameBox("mc.nii");
    finishedProcess();
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
    auto *process = new QProcess();

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

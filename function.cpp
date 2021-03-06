#include <QtWidgets>
#include "petmrMain.h"
#include "ImageIO.h"

void MainWindow::createfMRIPage()
{
    _fmriPage = new QWidget();

    _fMRIRunItemsBox = new QListWidget();
    _fMRIRunItemsBox->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::MinimumExpanding);
    _fMRIRunItemsBox->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    connect(_fMRIRunItemsBox, SIGNAL(itemClicked(QListWidgetItem*)),
            this, SLOT(changedfMRIRunCheckBox(QListWidgetItem*)));

    auto *runsLayout = new QVBoxLayout();
    runsLayout->addWidget(_fMRIRunItemsBox);

    auto *templateLabel = new QLabel("template directory");
    auto *fileLabel     = new QLabel("file name");
    _fMRITemplateDirectoryBox = new QComboBox();
    _fMRIFileNameBox          = new QComboBox();
    connect(_fMRITemplateDirectoryBox, SIGNAL(activated(int)),
            this, SLOT(changefMRITemplateDirectory(int)));

    auto *fileLayout = new QGridLayout();
    fileLayout->addWidget(templateLabel,0,0);
    fileLayout->addWidget(_fMRITemplateDirectoryBox,0,1);
    fileLayout->addWidget(fileLabel,1,0);
    fileLayout->addWidget(_fMRIFileNameBox,1,1);

    auto *setupLayout = new QVBoxLayout();
    setupLayout->addLayout(runsLayout);
    setupLayout->addLayout(fileLayout);

    auto *runsBox = new QGroupBox("List of EPI runs to include in analysis");
    runsBox->setLayout(setupLayout);

    _resliceEPIButton       = new QPushButton("reslice runs as necessary",_anatomyPage);
    _motionCorrectEPIButton = new QPushButton("motion-correct runs",_anatomyPage);
    _alignEPIButton         = new QPushButton("Align to template",_anatomyPage);
    connect(_resliceEPIButton,       SIGNAL(pressed()), this, SLOT(resliceEPI()));
//    connect(_motionCorrectEPIButton, SIGNAL(pressed()), this, SLOT(motionCorrectEPI()));
    connect(_alignEPIButton,         SIGNAL(pressed()), this, SLOT(alignEPI()));

    auto *actionLayout = new QVBoxLayout();
    actionLayout->addWidget(_resliceEPIButton);
    actionLayout->addWidget(_motionCorrectEPIButton);
    actionLayout->addWidget(_alignEPIButton);

    auto *actionBox = new QGroupBox("Process EPI runs");
    actionBox->setLayout(actionLayout);

    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(runsBox);
    pageLayout->addWidget(actionBox);
    _fmriPage->setLayout(pageLayout);
}

void MainWindow::changedfMRIRunCheckBox(QListWidgetItem* item)
{
    FUNC_ENTER;
    int iSelected=-1;
    for ( int jItem=0; jItem<_fMRIRunItems.size(); jItem++)
    {
        if ( item == &_fMRIRunItems.at(jItem) )
            iSelected = jItem;
    }
    if ( !_fMRIRunItems[iSelected].checkState() &&
         iSelected == _fMRITemplateDirectoryBox->currentIndex() )
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
        _fMRITemplateDirectoryBox->setCurrentIndex(iSelected);
        changefMRITemplateDirectory(iSelected);
    }
    enableEPIActionButtons();
}

void MainWindow::changefMRITemplateDirectory(int indexInBox)
{
    FUNC_ENTER << indexInBox << _fMRIFiles.size();
    _dimEPITemplate = _fMRIFiles[indexInBox].dim;
    FUNC_INFO << "_dimEPITemplate" << _dimEPITemplate.x << _dimEPITemplate.y << _dimEPITemplate.z;
    for (int jList=0; jList<_fMRIRunItems.size(); jList++)
    {
        if ( jList == indexInBox && _fMRIRunItems[jList].checkState() )
            _fMRIRunItems[jList].setBackgroundColor(Qt::yellow);
        else
            _fMRIRunItems[jList].setBackgroundColor(Qt::white);
    }
}

void MainWindow::openedfMRIPage()
{
    FUNC_ENTER;

    QDir const fMRITopDir("./epi");
    if (!fMRITopDir.exists())
    {
        _fMRIRunItems.clear();
        return;
    }

    FUNC_INFO << 1;
    QStringList const folderList = fMRITopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);

    // fMRI template directory combo box
    _fMRITemplateDirectoryBox->clear();
    for (int jList=0; jList<folderList.size(); jList++)
        _fMRITemplateDirectoryBox->addItem(folderList.at(jList));

    // Given the fMRI template directory, pick a file type
    // For the first selection, show all NIFTI files
    QString path = "./epi/" + _fMRITemplateDirectoryBox->currentText();
    QDir templateDir(path);
    FUNC_INFO << "path" << path;
    templateDir.setNameFilters(QStringList()<<"*.nii");
    QStringList fileList = templateDir.entryList();
    FUNC_INFO << fileList;
    _fMRIFileNameBox->clear();
    int indexRaw=-1;  int indexMC=-1;
    for (int jList=0; jList<fileList.size(); jList++)
    {
        _fMRIFileNameBox->addItem(fileList.at(jList));
        if ( fileList.at(jList) == "mc.nii")  indexMC = jList;
        if ( fileList.at(jList) == "raw.nii") indexRaw   = jList;
    }
    if ( indexMC >= 0 )
        _fMRIFileNameBox->setCurrentIndex(indexMC);
    else if ( indexRaw >= 0 )
        _fMRIFileNameBox->setCurrentIndex(indexRaw);
    else
        _fMRIFileNameBox->setCurrentIndex(0);

    FUNC_INFO << 2;
    QString fileName = path + "/" + _fMRIFileNameBox->currentText();
    getDimensions(fileName, _dimEPITemplate);

    _fMRIRunItems.resize(folderList.size());
    _fMRIFiles.clear();
    for (int jList=0; jList<folderList.size(); jList++)
    {
        FourDFile file;
        file.name = "epi/" + folderList.at(jList) + "/" + _fMRIFileNameBox->currentText();
        FUNC_INFO << fileName << file.name;
        QString text = getDimensions(file.name, file.dim);
        _fMRIFiles.append(file);
        _fMRIRunItems[jList].setText(folderList.at(jList) + " :    " + text);
        _fMRIRunItems[jList].setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
        _fMRIRunItems[jList].setCheckState(Qt::Checked);
        _fMRIRunItems[jList].setHidden(false);
        _fMRIRunItemsBox->addItem(&_fMRIRunItems[jList]);
    }
    FUNC_INFO << "_fMRIFiles size" << _fMRIFiles.size();
    _fMRITemplateDirectoryBox->setCurrentIndex(0);
    changefMRITemplateDirectory(0);

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
        if ( includeFile )
        {
            FourDFile file = _fMRIFiles[jList];
            allSameDimension &= file.dim.x == _dimEPITemplate.x;
            allSameDimension &= file.dim.y == _dimEPITemplate.y;
            allSameDimension &= file.dim.z == _dimEPITemplate.z;
        }
    }
    _resliceEPIButton->setEnabled(!allSameDimension);
    _motionCorrectEPIButton->setEnabled(allSameDimension);
    _alignEPIButton->setEnabled(allSameDimension);
    FUNC_EXIT << allSameDimension;
}

QString MainWindow::getDimensions(QString fileName, iPoint4D &dim)
{
    ImageIO file;
    if ( !file.readFileHeader(fileName,false) )
    {
        dim = file.getDimensions();
        QString text = QString("%1 x %2 x %3 with %4 time points").arg(dim.x).arg(dim.y).arg(dim.z).arg(dim.t);
        return text;
    }
    else
        return "";
}

void MainWindow::resliceEPI()
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
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)),
            this, SLOT(finishedFMResliceEPI(int, QProcess::ExitStatus)));

    QStringList arguments;
    arguments.append("-O");
    arguments.append("reslice");
    arguments.append("-a");
    QString templateFileName = "./epi/" + _fMRITemplateDirectoryBox->currentText()
            + "/" + _fMRIFileNameBox->currentText();
    arguments.append(templateFileName);
    arguments.append("--preprocess");
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
    QString processName = "/Users/jbm/QtApps/build-FM-Desktop_Qt_5_12_2_clang_64bit-Release/FM.app/Contents/MacOS/FM";
    FUNC_INFO << processName << arguments;
    process->start(processName,arguments);

    FUNC_EXIT;
}
void MainWindow::finishedFMResliceEPI(int exitCode, QProcess::ExitStatus exitStatus )
{
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    _centralWidget->setEnabled(true);
}

void MainWindow::alignEPI()
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
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)),
            this, SLOT(finishedFMAlignEPI(int, QProcess::ExitStatus)));

    QStringList arguments;
    arguments.append("-O");
    arguments.append("align");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->text());
    for (int jFile=0; jFile<_fMRIFiles.size(); jFile++)
    {
        bool includeFile = _fMRIRunItems[jFile].checkState();
        if ( includeFile )
            arguments.append(_fMRIFiles[jFile].name);
    }
    QString processName = "/Users/jbm/QtApps/build-FM-Desktop_Qt_5_12_2_clang_64bit-Release/FM.app/Contents/MacOS/FM";
    FUNC_INFO << processName << arguments;
    process->start(processName,arguments);

    FUNC_EXIT;
}
void MainWindow::finishedFMAlignEPI(int exitCode, QProcess::ExitStatus exitStatus )
{
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    _centralWidget->setEnabled(true);
}

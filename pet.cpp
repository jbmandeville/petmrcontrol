#include <QtWidgets>
#include "petmrMain.h"
#include "ImageIO.h"

void MainWindow::createPETPage()
{
    _petPage = new QWidget();

    auto *petRunLabel = new QLabel("PET run");
    _petDirBox = new QComboBox();
    connect(_petDirBox, SIGNAL(currentIndexChanged(int)),this, SLOT(updatePETDirBox(int)));
    _petDirBox->setToolTip("Sub-directory for PET that contains raw, mc, or align NIFTI");
    _petDirBox->setToolTip("File name: raw, mc, or align NIFTI");
//    connect(_petDirBox, SIGNAL(activated(int)),this, SLOT(updatePETDirBox(int)));

    auto *petFileName = new QLabel("file name");
    _petFileBox = new QComboBox();
    connect(_petFileBox, SIGNAL(activated(int)),this, SLOT(enablePETActionButtons()));

    auto *petRunLayout = new QGridLayout();
    petRunLayout->addWidget(petRunLabel,0,0);
    petRunLayout->addWidget(_petDirBox,0,1);
    petRunLayout->addWidget(petFileName,1,0);
    petRunLayout->addWidget(_petFileBox,1,1);

    auto *displayButton = new QPushButton("display file (fastmap)");
//    connect(displayButton, SIGNAL(pressed()), this, SLOT(displayPET()));

    auto *displayLayout = new QVBoxLayout();
    displayLayout->addWidget(displayButton);

    auto *runsBox = new QGroupBox("PET scan & file name");
    runsBox->setLayout(petRunLayout);
    runsBox->setLayout(displayLayout);

    _petFramesBox = new QListWidget();
    _petFramesBox->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::MinimumExpanding);
    _petFramesBox->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    connect(_petFramesBox, SIGNAL(itemClicked(QListWidgetItem*)),this, SLOT(changedPETFrameSelection(QListWidgetItem*)));

    auto *runLayout = new QVBoxLayout();
    runLayout->addLayout(petRunLayout);
    runLayout->addLayout(displayLayout);

    auto *templateLabel = new QLabel("fMRI template directory (with mcTemplate.nii)");
    auto *fileLabel     = new QLabel("fMRI file name(s) for motion-correction");
    _fMRIForPETTemplate = new QLabel("");
    _fMRIForPETFileName = new QLabel("");

    QString smoothText; smoothText.setNum(_savedSettings.fmSmoothing);
    _smoothingPET = new QLineEdit(smoothText);
    _smoothingPET->setMaximumWidth(150);
    connect(_smoothingPET, SIGNAL(editingFinished()), this, SLOT(changedSmoothingPET()));

    auto *fMRILayout = new QGridLayout();
    fMRILayout->addWidget(templateLabel,0,0);
    fMRILayout->addWidget(_fMRIForPETTemplate,0,1);
    fMRILayout->addWidget(fileLabel,1,0);
    fMRILayout->addWidget(_fMRIForPETFileName,1,1);

    auto *mcLayout = new QVBoxLayout();
    runLayout->addWidget(_petFramesBox);
    mcLayout->addLayout(runLayout);
    mcLayout->addLayout(fMRILayout);

    _mcPETBox = new QGroupBox("Motion-correct PET using EPI");
    _mcPETBox->setLayout(mcLayout);

    _doEverythingPETButton = new QPushButton("Do everything (no interaction required)");
    _doEverythingPETButton->setCheckable(true);
    connect(_doEverythingPETButton, SIGNAL(clicked(bool)), this, SLOT(doEverthingPET()));

    _motionCorrectMatchingMRIButton = new QPushButton("Create matching MRI and motion-correct it");
    connect(_motionCorrectMatchingMRIButton, SIGNAL(pressed()), this, SLOT(motionCorrectMatchingMRI()));

    _motionCorrectPETButton = new QPushButton("Apply motion-correction to PET (raw -->mc)");
    connect(_motionCorrectPETButton,         SIGNAL(pressed()), this, SLOT(applyMotionCorrectionToPET()));

    _alignPETButton   = new QPushButton("Align PET time series (mc --> align)");
    connect(_alignPETButton,                 SIGNAL(pressed()), this, SLOT(alignPET()));

    _analyzeTAC       = new QPushButton("Create kinetic analysis");
    connect(_analyzeTAC,                     SIGNAL(pressed()), this, SLOT(analyzeTAC()));

    auto *smoothingLabel = new QLabel("post-alignment smoothing width");
    _smoothingPET = new QLineEdit("0.");
    _smoothingPET->setMaximumWidth(150);
    connect(_smoothingPET, SIGNAL(editingFinished()), this, SLOT(changedSmoothingPET()));
    auto *alignLayout = new QGridLayout();
    alignLayout->addWidget(smoothingLabel,0,0);
    alignLayout->addWidget(_smoothingPET,0,1);

    auto *actionLayout = new QVBoxLayout();
    actionLayout->addWidget(_doEverythingPETButton);
    actionLayout->addWidget(_motionCorrectMatchingMRIButton);
    actionLayout->addWidget(_motionCorrectPETButton);
    actionLayout->addLayout(alignLayout);
    actionLayout->addWidget(_alignPETButton);
    actionLayout->addWidget(_analyzeTAC);

    auto *actionBox = new QGroupBox("Align PET to template using MRI");
    actionBox->setLayout(actionLayout);

    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(runsBox);
    pageLayout->addWidget(_mcPETBox);
    pageLayout->addWidget(actionBox);
    _petPage->setLayout(pageLayout);
}

void MainWindow::updatePETFileNameBox(QString fileName)
{
    FUNC_ENTER << fileName;
    updatePETFileNameBox();
    int foundIndex=-1;
    for (int jList=0; jList<_petFileBox->count(); jList++)
    {
        if ( !_petFileBox->itemText(jList).compare(fileName) )
            foundIndex = jList;
    }
    if ( foundIndex >= 0 ) _petFileBox->setCurrentIndex(foundIndex);
    FUNC_EXIT << foundIndex;
}


void MainWindow::updatePETFileNameBox()
{
    // For the current selection, show all NIFTI files
    QString path = "./pet/" + _petDirBox->currentText();
    QDir petDir(path);
    petDir.setNameFilters(QStringList()<<"raw.nii"<<"mc.nii"<<"align.nii");
    QStringList fileList = petDir.entryList();

    _petFileBox->clear();
    int indexRaw=-1;  int indexMC=-1;  int indexAlign=-1;
    for (int jList=0; jList<fileList.size(); jList++)
    {
        _petFileBox->addItem(fileList.at(jList));
        if ( fileList.at(jList) == "raw.nii")     indexRaw      = jList;
        if ( fileList.at(jList) == "mc.nii")      indexMC      = jList;
        if ( fileList.at(jList) == "align.nii")   indexAlign   = jList;
    }
    if ( indexAlign >= 0 )
        _petFileBox->setCurrentIndex(indexAlign);
    else if ( indexMC >= 0)
        _petFileBox->setCurrentIndex(indexMC);
    else if ( indexRaw >= 0)
        _petFileBox->setCurrentIndex(indexRaw);
    else if ( _anatomyFileBox->count() > 0 )
        _petFileBox->setCurrentIndex(0);
    // Enable/disable:  // if anatomy file was found, it can be used for either freeSurfer or alignment
    enablePETActionButtons();
}


void MainWindow::changedPETFrameSelection(QListWidgetItem *item)
{
    FUNC_ENTER;
    int iFrame=-1;
    for ( int jItem=0; jItem<_petFrameItems.size(); jItem++)
    {
        if ( item == &_petFrameItems.at(jItem) )
            iFrame = jItem;
    }
    double timePET = _petFile.timeTags.at(iFrame);
    dPoint2D timeFrame = petFrameTime(iFrame);
    qInfo() << QString("PET frame %d has (start, center, end) at times (%1, %2, %3)").arg(timeFrame.lower).arg(timePET).arg(timeFrame.upper);
    if ( _matchingEPI.size() > 0 )
    {
        if ( _matchingEPI[iFrame].size() == 0 )
            qInfo() << "No overlap with EPI: motion-correction parameters will be interpolated for this point";
        else
        {
            int currentRun = _matchingEPI[iFrame].at(0).x;
            qInfo() << "Overlap with run" << _fMRIFilesForPETMC.at(currentRun).name;
            for (int jPair=0; jPair<_matchingEPI[iFrame].size(); jPair++)
            {
                if  ( _matchingEPI[iFrame].at(jPair).x != currentRun )
                {
                    currentRun = _matchingEPI[iFrame].at(jPair).x;
                    qInfo() << "Overlap woith run" << _fMRIFilesForPETMC.at(currentRun).name;
                }
                qInfo() << "time point" << _matchingEPI[iFrame].at(jPair).y;
            }
        }
    }
}

bool MainWindow::petFileExists(QString dirName, QString fileName)
{
    bool fileExists = false;
    QString fullName = "pet/" + dirName + "/" + fileName;
    QFileInfo checkFile(fullName);
    if ( checkFile.exists() && checkFile.isFile() )
         fileExists = true;
    return fileExists;
}

void MainWindow::openedPETPage()
{
    FUNC_ENTER;

    // Make sure anatomy is defined
    openedAnatomyPage();

    QDir const petTopDir("./pet");
    if ( !petTopDir.exists() )
    {
        _petDirBox->clear();
        return;
    }
    QStringList const petFolderList = petTopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    FUNC_INFO << petFolderList;
    _petDirBox->clear();
    for (int jList=0; jList<petFolderList.size(); jList++)
    {
        if ( petFileExists(petFolderList.at(jList),"time-tags.txt") ||
             petFileExists(petFolderList.at(jList),"raw.nii") ||
             petFileExists(petFolderList.at(jList),"align.nii") )
        {
            _petDirBox->addItem(petFolderList.at(jList));
            FUNC_INFO << "add to petDirBox" << petFolderList.at(jList);
        }
    }
    _petDirBox->setCurrentIndex(_petDirBox->count()-1);
    FUNC_INFO << "petDirBox text #1" << _petDirBox->currentText() << _petDirBox->currentIndex();

    QDir const fMRITopDir("./epi");
    if (!fMRITopDir.exists())
        return;
    FUNC_INFO << 1;
    QStringList const epiFolderList = fMRITopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);

    _fMRIFilesForPETMC.clear();
    for (int jList=0; jList<epiFolderList.size(); jList++)
    {
        FourDFile fMRIFile;
        fMRIFile.name = epiFolderList.at(jList);
        QString tagsFileName = fMRITopDir.absolutePath() + "/" + fMRIFile.name + "/time-tags.txt";
        // Because only the time dimension is needed here for the fMRI files, just use the time-tags to avoid file searches
        getTimeTags(tagsFileName,fMRIFile.timeTags,fMRIFile.timeText);
        fMRIFile.dim.t = fMRIFile.timeTags.size();
        if ( fMRIFile.dim.t > 0 )
            _fMRIFilesForPETMC.append(fMRIFile);
    }

    FUNC_INFO << 2 << _fMRIFilesForPETMC.size();
    // find which one is the template
    int iTemplate=-1;
    for (int jList=0; jList<_fMRIFilesForPETMC.size(); jList++)
    {
        QString fileName = fMRITopDir.absolutePath() + "/" + _fMRIFilesForPETMC.at(jList).name + "/mcTemplate.nii";
        QFileInfo checkFile(fileName);
        FUNC_INFO << "check template" << fileName;
        if ( checkFile.exists() && checkFile.isFile() )
            iTemplate = jList;
    }
    iTemplate = qMax(0,iTemplate);
    _fMRIForPETTemplate->setText(_fMRIFilesForPETMC.at(iTemplate).name);

    FUNC_INFO << 3;
    if ( !_fMRIForPETTemplate->text().isEmpty() )
    {
        QString fileName = fMRITopDir.absolutePath() + "/" + _fMRIForPETTemplate->text() + "/reslice.nii";
        QFileInfo checkReslice(fileName);
        if (checkReslice.exists() && checkReslice.isFile())
            _fMRIForPETFileName->setText("reslice.nii");
        else
        {
            fileName = fMRITopDir.absolutePath() + "/" + _fMRIForPETTemplate->text() + "/reslice.nii";
            QFileInfo checkRaw(fileName);
            if (checkRaw.exists() && checkRaw.isFile())
                _fMRIForPETFileName->setText("raw.nii");
            else
                _fMRIForPETTemplate->setText("");  // this should disable action button
        }
    }
    FUNC_INFO << 4;

    FUNC_INFO << "petDirBox text #1" << _petDirBox->currentText() << _petDirBox->currentIndex();
    updatePETDirBox(_petDirBox->currentIndex());

    readSmoothing(0);
    readSmoothing(1);
    readSmoothing(2);

    enablePETActionButtons();

    FUNC_EXIT;
}

void MainWindow::enablePETActionButtons()
{
    if ( !_petFileBox->currentText().compare("raw.nii") )
        _alignFileNameForPETRegistration = "t1/" + _anatomyDirBox->currentText() + "/align.com";
    else if ( !_petFileBox->currentText().compare("align.nii") )
        _alignFileNameForPETRegistration = "pet/" + _petDirBox->currentText() + "/align.com";
    _alignPETButton->setText(QString("Align PET time series (mc --> align) using %1")
                             .arg(_alignFileNameForPETRegistration));

    // enable
    bool bay6 = _radioButtonNHPBay6->isChecked();
    bool timeTagsExist = petFileExists("time-tags.txt");

    if (  bay6 )
        _alignPETButton->setEnabled(petFileExists("raw.nii"));
    else
    {
        _motionCorrectMatchingMRIButton->setEnabled( !_fMRIForPETTemplate->text().isEmpty() && _petDirBox->count() > 0);
        _motionCorrectPETButton->setEnabled(_petDirBox->count() > 0 && petFileExists("mc.dat"));
        _doEverythingPETButton->setEnabled(petFileExists("raw.nii") && anatomyFileExists("align.com") &&
                                           !_fMRIForPETTemplate->text().isEmpty() );
        _alignPETButton->setEnabled(petFileExists("mc.nii")    && anatomyFileExists("align.com"));
    }
    _motionCorrectMatchingMRIButton->setVisible(!bay6);
    _motionCorrectPETButton->setVisible(!bay6);
    _doEverythingPETButton->setVisible(!bay6);
    _mcPETBox->setVisible(!bay6 && timeTagsExist);

    FUNC_INFO << "tests" << !bay6 << timeTagsExist;

    _analyzeTAC->setEnabled(petFileExists("align.nii"));

    // highlight
    if ( _petDirBox->count() > 0 && petFileExists("mc.dat") )
        _motionCorrectMatchingMRIButton->setStyleSheet("background-color:lightYellow;");
    if ( petFileExists("mc.nii") || petFileExists("align.nii") )
        _motionCorrectPETButton->setStyleSheet("background-color:lightYellow;");
    if ( petFileExists("align.nii") )
        _alignPETButton->setStyleSheet("background-color:lightYellow;");
    if ( petFileExists("srtm","timeModel.dat") )
    {
        _analyzeTAC->setStyleSheet("background-color:lightYellow;");
        _doEverythingPETButton->setStyleSheet("background-color:lightYellow;");
    }

}

void MainWindow::updatePETDirBox(int indexInBox)
{
    _petFile.name = "pet/" + _petDirBox->currentText() + "/" + _petFileBox->currentText();
    FUNC_ENTER << indexInBox << _petFile.name;
    if ( petFileExists(_petFileBox->currentText()) )
    {
        QString fileName = "pet/" + _petDirBox->currentText() + "/time-tags.txt";
        getTimeTags(fileName,_petFile.timeTags,_petFile.timeText);
        getDimensions(_petFile.name, _petFile.dim);
        FUNC_INFO << "dim.t" << _petFile.dim.t;
        findPETandFMRIOverlap();
        FUNC_INFO << "sizes" << _petFile.dim.t << _petFile.timeText.size() << _matchingEPI.size();
        _petFrameItems.resize(_petFile.dim.t);
        for (int jFrame=0; jFrame<_petFile.dim.t; jFrame++)
        {
            FUNC_INFO << "jFrame" << jFrame;
            QString text;
            if ( jFrame < _petFile.timeText.size() )
            {
                int nOverlap = _matchingEPI[jFrame].count();
                dPoint2D frameTime = petFrameTime(jFrame);
                double width = (frameTime.upper - frameTime.lower);
                text = QString("%1: start %2, width %3 sec, overlap with %4 fMRI points").arg(jFrame+1).arg(_petFile.timeText.at(jFrame))
                        .arg(width).arg(nOverlap);
            }
            else
                text = QString("%1").arg(jFrame+1);
            QListWidgetItem item;
            item.setText(text);
            _petFrameItems[jFrame] = item;
            _petFramesBox->addItem(&_petFrameItems[jFrame]);
        }
    }
    updatePETFileNameBox();
    FUNC_EXIT;
}
void MainWindow::findPETandFMRIOverlap()
{
    FUNC_ENTER;
    _matchingEPI.resize(_petFile.dim.t);
    for (int jFrame=0; jFrame<_petFile.dim.t; jFrame++)
    {
        dPoint2D timeFrame = petFrameTime(jFrame);
        for (int jFile=0; jFile<_fMRIFilesForPETMC.size(); jFile++)
        {
            FourDFile fmri = _fMRIFilesForPETMC.at(jFile);
            for (int jt=0; jt<fmri.dim.t; jt++)
            {
                double timeEPI = fmri.timeTags.at(jt);
                if ( timeEPI >= timeFrame.lower && timeEPI <= timeFrame.upper )
                {
                    iPoint2D pair; pair.x=jFile;  pair.y=jt;
                    _matchingEPI[jFrame].append(pair);
                }
            }
        }
        // output
        FUNC_INFO << "matching EPI for frame" << jFrame+1;
        for (int jPair=0; jPair<_matchingEPI[jFrame].size(); jPair++)
        {
            int iRun  = _matchingEPI[jFrame].at(jPair).x;
            int iTime = _matchingEPI[jFrame].at(jPair).y;
            FUNC_INFO << "(" << _fMRIFilesForPETMC.at(iRun).name << "," << iTime << ")";
        }
    }
//    writeJipCommandFileForMatchingMRI();
    FUNC_EXIT;
}

void MainWindow::doEverthingPET()
{
    FUNC_ENTER;
    _motionCorrectMatchingMRIButton->setStyleSheet("background-color:white;");
    _motionCorrectPETButton->setStyleSheet("background-color:white;");
    _alignPETButton->setStyleSheet("background-color:white;");
    _analyzeTAC->setStyleSheet("background-color:white;");
    _doEverythingPETButton->setStyleSheet("background-color:white;");

    motionCorrectMatchingMRI();
}

void MainWindow::motionCorrectMatchingMRI()
{
    FUNC_ENTER;
    updatePETFileNameBox("raw.nii");  // not really required
    // First write the jip command file to create the matching MRI volume
    writeJipCommandFileForMatchingMRI();

    auto *process = new QProcess();
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedMotionCorrectMatchingMRI(int, QProcess::ExitStatus)));
    process->setProcessChannelMode(QProcess::MergedChannels);

    QString exe = _scriptDirectory + "createMatchingEPIandMCForPET.csh";
    QStringList arguments;
    // arguments: [EPI MC template dir] [pet dir]
    arguments.append(_fMRIForPETTemplate->text());
    arguments.append(_petDirBox->currentText());
    QString message = "Motion-correct PET-matched EPI volumes; this takes about 10 minutes";
    spawnProcess(process,exe,arguments,message,"Motion-correct matching EPI");
    FUNC_EXIT;
}
void MainWindow::finishedMotionCorrectMatchingMRI(int exitCode, QProcess::ExitStatus exitStatus)
{
    FUNC_ENTER;

    // Is interpolation required?
    finishedProcess();
    if ( _doEverythingPETButton->isChecked() )
    {
        if ( _petDirBox->count() > 0 && petFileExists("mc.dat") )
             _motionCorrectMatchingMRIButton->setStyleSheet("background-color:lightYellow;");
        applyMotionCorrectionToPET();
    }
    else
        enablePETActionButtons();
}
void MainWindow::applyMotionCorrectionToPET()
{
    FUNC_ENTER;

    auto *process = new QProcess();
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedApplyingMCToPET(int, QProcess::ExitStatus)));
    process->setProcessChannelMode(QProcess::MergedChannels);

    QString exe = _scriptDirectory + "applyMCToPET.csh";
    QStringList arguments;
    // arguments: [pet dir] [mc file name]
    arguments.append(_petDirBox->currentText());
    arguments.append("mc.dat");
    QString message = "Apply EPI-based motion-correction to PET; this takes about 10 minutes";
    spawnProcess(process,exe,arguments,message,"Apply motion correction to PET");
    FUNC_EXIT;
}
void MainWindow::finishedApplyingMCToPET(int exitCode, QProcess::ExitStatus exitStatus)
{
    FUNC_ENTER;
    finishedProcess();
    updatePETFileNameBox("mc.nii");  // not really required
    if ( _doEverythingPETButton->isChecked() )
    {
        if ( petFileExists("mc.nii") || petFileExists("align.nii") )
            _motionCorrectPETButton->setStyleSheet("background-color:lightYellow;");
        alignPET();
    }
    else
        enablePETActionButtons();
}

void MainWindow::writeJipCommandFileForMatchingMRI()
{
    FUNC_ENTER;
    QString fileName = "pet/" + _petDirBox->currentText() + "/jip-createMRI.com";
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);

    // Find all required epi runs; all vectors below contains indices into the list of runs _fMRIFilesForPETMC
    for (int jFrame=0; jFrame<_matchingEPI.size(); jFrame++)
    {
        iVector firstPoint;  iVector lastPoint;
        int nEPIFiles = _fMRIFilesForPETMC.size();
        firstPoint.fill(-1,nEPIFiles);
        lastPoint.fill(-1,nEPIFiles);
        for (int jPair=0; jPair<_matchingEPI[jFrame].size(); jPair++)
        {
            int iRun  = _matchingEPI[jFrame].at(jPair).x;  // index into _fMRIFilesForPETMC
            int iTime = _matchingEPI[jFrame].at(jPair).y;
            if ( iTime >= 0 && firstPoint[iRun] < 0 ) firstPoint[iRun] = iTime;
            if ( iTime > lastPoint[iRun] ) lastPoint[iRun] = iTime;
        }

        out << "# create matching MRI volume #" << jFrame+1 << "\n";
        if ( _matchingEPI[jFrame].size() > 0 )
        {
            // Read all required files
            for (int jFile=0; jFile<_fMRIFilesForPETMC.size(); jFile++)
            {
                QString dir = _fMRIFilesForPETMC.at(jFile).name;
                QString name = "../../epi/" + dir + "/" + _fMRIForPETFileName->text();
                if ( firstPoint.at(jFile) >= 0 )
                    out << "read " << name << ">" << firstPoint.at(jFile) << "-" << lastPoint.at(jFile) << " " << dir << "\n";
            }

            // average all points (potentially more than 1 file) to create an average volume
            out << "average volume";
            // average to create volume
            for (int jFile=0; jFile<_fMRIFilesForPETMC.size(); jFile++)
            {
                QString dir = _fMRIFilesForPETMC.at(jFile).name;
                if ( firstPoint.at(jFile) >= 0 )
                    out << " " << dir;
            }
            out << "\n";

            // write the volume
            out << "write matchingMRI.nii volume\n";
            if ( jFrame != _matchingEPI.size()-1 )
                out << "delete\n\n";
        }
    }
    out << "bye\n";
    file.close();
}

bool MainWindow::getPETMCInterpolationRequired()
{
    bool interpolationRequired = false;
    // Find all required epi runs; all vectors below contains indices into the list of runs _fMRIFilesForPETMC
    for (int jFrame=0; jFrame<_matchingEPI.size(); jFrame++)
        interpolationRequired |=  _matchingEPI[jFrame].size() == 0;
}

void MainWindow::alignPET()
{
    FUNC_ENTER;
    auto *process = new QProcess();
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedFMAlignPET(int, QProcess::ExitStatus)));

    QStringList arguments;
    QString mcFileName = "pet/" + _petDirBox->currentText() + "/raw.nii";
    arguments.append(mcFileName);
    arguments.append("-O");
    arguments.append("alignment");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());
    arguments.append("--alignment-file");
    arguments.append(_alignFileNameForPETRegistration);
    arguments.append("--smoothing");
    arguments.append(_smoothingPET->text());
    arguments.append("--output-file");
    arguments.append("align");
    if ( !_petFileBox->currentText().compare("raw.nii") ||
         !_petFileBox->currentText().compare("mc.nii") )
    {
        arguments.append("--preprocess");
        arguments.append("register");
        arguments.append("--quit");
    }

    QString message = "Align/register PET to template space; this takes about 10 minutes";
    spawnProcess(process,_fastmapProcess,arguments,message,"");

    FUNC_EXIT;
}

void MainWindow::finishedFMAlignPET(int exitCode, QProcess::ExitStatus exitStatus )
{
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    finishedProcess();
    updatePETFileNameBox("align.nii");  // not really required
    if ( _doEverythingPETButton->isChecked() )
    {
        if ( petFileExists("align.nii") )
            _alignPETButton->setStyleSheet("background-color:lightYellow;");
        _doEverythingPETButton->setChecked(false);
        analyzeTAC();
    }
    else
        enablePETActionButtons();
}

void MainWindow::analyzeTAC()
{
    FUNC_ENTER;

    // Create a analysis directory on pet ("srtm") and install files
    installSRTMAnalysis();

    auto *process = new QProcess();
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedFMAnalyzeTAC(int, QProcess::ExitStatus)));

    QStringList arguments;
    arguments.append("-t");
    QString timeFile = "timeModel.dat";
    arguments.append(timeFile);
    arguments.append("-O");
    arguments.append("pet");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());

    process->setWorkingDirectory("pet/srtm");
    QString message = "Create a new kinetic analysis; this requires interaction";
    spawnProcess(process,_fastmapProcess,arguments,message,"");

    FUNC_EXIT;
}
void MainWindow::finishedFMAnalyzeTAC(int exitCode, QProcess::ExitStatus exitStatus )
{
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    enablePETActionButtons();
    finishedProcess();
    if ( petFileExists("srtm","timeModel.dat") )
    {
        _analyzeTAC->setStyleSheet("background-color:lightYellow;");
        _doEverythingPETButton->setStyleSheet("background-color:lightYellow;");
    }
}
void MainWindow::installSRTMAnalysis()
{
    // write the pet frames table to the pet directory
    writeFramesTable("pet/frames.table");

    // create the directory if needed
    QString directory = "pet/srtm";
    QDir dir(directory);
    if ( !dir.exists() ) dir.mkpath(directory);

    writeTimeModelFile(directory);
    writeGLMFile(directory);
}
void MainWindow::writeFramesTable(QString fileName)
{
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);

    out << "dt\n";

    for (int jFrame=0; jFrame<_petFile.timeTags.size(); jFrame++)
    {
        dPoint2D timeFrame = petFrameTime(jFrame);
        double width = (timeFrame.upper - timeFrame.lower);
        out << width << "\n";
    }
    file.close();
}
void MainWindow::writeTimeModelFile(QString directoryName)
{
    FUNC_ENTER << directoryName;
    QString fileName = directoryName + "/timeModel.dat";
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);

    out << "time-model SRTM\n";
    out << "conditions B\n";
    out << "scans:\n";
    out << "../"+_petDirBox->currentText()+"/align.nii pet1.glm ../frames.table\n";
    file.close();
}
void MainWindow::writeGLMFile(QString directoryName)
{
    FUNC_ENTER << directoryName;
    QString fileName = directoryName + "/pet1.glm";
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);

    out << "R1 1\n\n";
    out << "k2 a\n\n";
    out << "k2a B\n";
    file.close();
}

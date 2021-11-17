#include <QtWidgets>
#include "petmrMain.h"
#include "ImageIO.h"

void MainWindow::createPETPage()
{
    _petPage = new QWidget();

    auto *petRunLabel = new QLabel("PET run");
    _petDirBox = new QComboBox();
    _petFramesBox = new QListWidget();
    _petFramesBox->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::MinimumExpanding);
    _petFramesBox->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    connect(_petFramesBox, SIGNAL(itemClicked(QListWidgetItem*)),this, SLOT(changedPETFrameSelection(QListWidgetItem*)));
    connect(_petDirBox, SIGNAL(activated(int)),this, SLOT(updatePETRunBox(int)));

    auto *petRunLayout = new QHBoxLayout();
    petRunLayout->addWidget(petRunLabel);
    petRunLayout->addWidget(_petDirBox);

    auto *runLayout = new QVBoxLayout();
    runLayout->addLayout(petRunLayout);
    runLayout->addWidget(_petFramesBox);

    auto *templateLabel = new QLabel("fMRI template directory (with mcTemplate.nii)");
    auto *fileLabel     = new QLabel("fMRI file name(s) for motion-correction");
    _fMRIForPETTemplate = new QLabel("");
    _fMRIForPETFileName = new QLabel("");

    _motionCorrectMatchingMRIButton = new QPushButton("Create matching MRI and motion-correct it");
    _motionCorrectMatchingMRIButton->setEnabled(false);
    connect(_motionCorrectMatchingMRIButton, SIGNAL(pressed()), this, SLOT(motionCorrectMatchingMRI()));

    _motionCorrectPETButton = new QPushButton("Apply motion-correction to PET (raw -->mc)");
    connect(_motionCorrectPETButton,         SIGNAL(pressed()), this, SLOT(applyMotionCorrectionToPET()));

    auto *fMRILayout = new QGridLayout();
    fMRILayout->addWidget(templateLabel,0,0);
    fMRILayout->addWidget(_fMRIForPETTemplate,0,1);
    fMRILayout->addWidget(fileLabel,1,0);
    fMRILayout->addWidget(_fMRIForPETFileName,1,1);

    auto *setupLayout = new QVBoxLayout();
    setupLayout->addLayout(runLayout);
    setupLayout->addLayout(fMRILayout);
    setupLayout->addWidget(_motionCorrectMatchingMRIButton);
    setupLayout->addWidget(_motionCorrectPETButton);

    auto *runsBox = new QGroupBox("Motion-correct PET using EPI");
    runsBox->setLayout(setupLayout);

    _reslicePETButton = new QPushButton("Reslice PET (mc -->reslice)");
    _alignPETButton   = new QPushButton("Align PET time series (reslice --> align)");
    _analyzeTAC       = new QPushButton("Analyze TAC");
    connect(_reslicePETButton,               SIGNAL(pressed()), this, SLOT(reslicePET()));
    connect(_alignPETButton,                 SIGNAL(pressed()), this, SLOT(alignPET()));
    connect(_analyzeTAC,                     SIGNAL(pressed()), this, SLOT(analyzeTAC()));
    _motionCorrectPETButton->setEnabled(false);

    auto *actionLayout = new QVBoxLayout();
    actionLayout->addWidget(_reslicePETButton);
    actionLayout->addWidget(_alignPETButton);
    actionLayout->addWidget(_analyzeTAC);

    auto *actionBox = new QGroupBox("Align PET to template using MRI");
    actionBox->setLayout(actionLayout);

    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(runsBox);
    pageLayout->addWidget(actionBox);
    _petPage->setLayout(pageLayout);
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
//        QString fullFileName = fMRITopDir.absolutePath() + "/" + fMRIFile.name + "/raw.nii";
//        getDimensions(fullFileName, fMRIFile.dim);
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
        QFileInfo checkFile(fileName);
        if (checkFile.exists() && checkFile.isFile())
            _fMRIForPETFileName->setText("reslice.nii");
        else
        {
            QString fileName = fMRITopDir.absolutePath() + "/" + _fMRIForPETTemplate->text() + "/raw.nii";
            QFileInfo checkFile(fileName);
            if (checkFile.exists() && checkFile.isFile())
                _fMRIForPETFileName->setText("raw.nii");
            else
                _fMRIForPETTemplate->setText("");  // this should disable action button
        }
    }
    FUNC_INFO << 4;

    FUNC_INFO << "petDirBox text #1" << _petDirBox->currentText() << _petDirBox->currentIndex();
    updatePETRunBox(_petDirBox->currentIndex());
    enablePETActionButtons();

    FUNC_EXIT;
}

void MainWindow::enablePETActionButtons()
{
    QString templateFileName = "t1/" + _anatomyDirBox->currentText() + "/brain.nii";
    QFileInfo checkBrain(templateFileName);
    if ( checkBrain.exists() && checkBrain.isFile() )
        _anatomyFileNameForPETReslice = templateFileName;
    else
    {
        templateFileName = "t1/" + _anatomyDirBox->currentText() + "/raw.nii";
        QFileInfo checkRaw(templateFileName);
        if ( checkRaw.exists() && checkRaw.isFile() )
            _anatomyFileNameForPETReslice = templateFileName;
    }
    _reslicePETButton->setText(QString("Reslice PET (mc -->reslice) using %1").arg(_anatomyFileNameForPETReslice));\

    _alignFileNameForPETRegistration = "t1/" + _anatomyDirBox->currentText() + "/align.com";
    _alignPETButton->setText(QString("Align PET time series (reslice --> align) using %1")
                             .arg(_alignFileNameForPETRegistration));

    // enable
    _motionCorrectMatchingMRIButton->setEnabled( !_fMRIForPETTemplate->text().isEmpty() && _petDirBox->count() > 0);
    _motionCorrectPETButton->setEnabled(_petDirBox->count() > 0 && petFileExists("mc.dat"));
    _reslicePETButton->setEnabled(petFileExists("mc.nii"));
    _alignPETButton->setEnabled(petFileExists("reslice.nii") && anatomyFileExists("align.com"));

    FUNC_INFO << "***" << petFileExists("reslice.nii") << anatomyFileExists("align.com");

    // highlight
    if ( _petDirBox->count() > 0 && petFileExists("mc.dat") )
        _motionCorrectMatchingMRIButton->setStyleSheet("background-color:lightYellow;");
    if ( petFileExists("mc.nii") || petFileExists("reslice.nii") || petFileExists("align.nii") )
        _motionCorrectPETButton->setStyleSheet("background-color:lightYellow;");
    if ( petFileExists("reslice.nii") )
        _reslicePETButton->setStyleSheet("background-color:lightYellow;");
    if ( petFileExists("align.nii") )
        _alignPETButton->setStyleSheet("background-color:lightYellow;");

}

void MainWindow::updatePETRunBox(int indexInBox)
{
    _petFile.name = "raw.nii";
    FUNC_ENTER << indexInBox << _petFile.name;
    if ( petFileExists(_petFile.name) )
    {
        getDimensions(_petFile.name, _petFile.dim);
        bool timeTagsExists = petFileExists("time-tags.txt");
        FUNC_INFO << "timeTagsExists" << timeTagsExists;
        _petFramesBox->setVisible(timeTagsExists);
        QString fileName = "pet/" + _petDirBox->currentText() + "/time-tags.txt";
        FUNC_INFO << "dim.t" << _petFile.dim.t;
        getTimeTags(fileName,_petFile.timeTags,_petFile.timeText);
        findPETandFMRIOverlap();
        FUNC_INFO << "sizes" << _petFile.dim.t << _petFile.timeText.size() << _matchingEPI.size();
        //        _petFramesBox->clear();
        FUNC_INFO << "here 1";
        _petFrameItems.resize(_petFile.dim.t);
        FUNC_INFO << "here 2";
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
                    _matchingEPI[jFrame].append({jFile,jt});
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

void MainWindow::motionCorrectMatchingMRI()
{
    FUNC_ENTER;
    // First write the jip command file to create the matching MRI volume
    writeJipCommandFileForMatchingMRI();

    auto *process = new QProcess;
    _outputBrowser->setWindowTitle("Motion-correct matching MRI");
    showBrowser(true);
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedMotionCorrectMatchingMRI(int, QProcess::ExitStatus)));
    process->setProcessChannelMode(QProcess::MergedChannels);

    QString exe = _scriptDirectory + "createMatchingEPIandMCForPET.csh";
    QStringList arguments;
    // arguments: [EPI MC template dir] [pet dir]
    arguments.append(_fMRIForPETTemplate->text());
    arguments.append(_petDirBox->currentText());
    qInfo() <<  exe << arguments;
    process->start(exe,arguments);
    _centralWidget->setEnabled(false);
    FUNC_EXIT;
}
void MainWindow::finishedMotionCorrectMatchingMRI(int exitCode, QProcess::ExitStatus exitStatus)
{
    FUNC_ENTER;
    enablePETActionButtons();

    // Is interpolation required?
    _centralWidget->setEnabled(true);
    showBrowser(false);
}
void MainWindow::applyMotionCorrectionToPET()
{
    FUNC_ENTER;

    auto *process = new QProcess;
    _outputBrowser->setWindowTitle("Apply motion correction to PET");
    showBrowser(true);
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedApplyingMCToPET(int, QProcess::ExitStatus)));
    process->setProcessChannelMode(QProcess::MergedChannels);

    QString exe = _scriptDirectory + "applyMCToPET.csh";
    QStringList arguments;
    // arguments: [pet dir] [mc file name]
    arguments.append(_petDirBox->currentText());
    arguments.append("mc.dat");
    qInfo() <<  exe << arguments;
    process->start(exe,arguments);
    _centralWidget->setEnabled(false);
    FUNC_EXIT;
}
void MainWindow::finishedApplyingMCToPET(int exitCode, QProcess::ExitStatus exitStatus)
{
    FUNC_ENTER;
    enablePETActionButtons();
    _centralWidget->setEnabled(true);
    showBrowser(false);
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

void MainWindow::reslicePET()
{
    FUNC_ENTER;
    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedFMReslicePET(int, QProcess::ExitStatus)));

    QStringList arguments;
    arguments.append("-O");
    arguments.append("reslice");
    arguments.append("-a");
    QString templateFileName = _anatomyFileNameForPETReslice;
    arguments.append(templateFileName);
    arguments.append("--preprocess");
    arguments.append("reslice");
    arguments.append("--output-file");
    arguments.append("reslice");
    arguments.append("--quit");
    QString mcFileName = "pet/" + _petDirBox->currentText() + "/mc.nii";
    arguments.append(mcFileName);
    FUNC_INFO << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    FUNC_EXIT;
}
void MainWindow::alignPET()
{
    FUNC_ENTER;
    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedFMAlignPET(int, QProcess::ExitStatus)));

    QStringList arguments;
    QString mcFileName = "pet/" + _petDirBox->currentText() + "/reslice.nii";
    arguments.append(mcFileName);
    arguments.append("-O");
    arguments.append("alignment");
    arguments.append("-T");
    arguments.append(_anatomyTemplateDirectory->currentText());
    // alignment file name
    arguments.append("--preprocess");
    arguments.append("register");
    arguments.append("--output-file");
    arguments.append("align");
    arguments.append("--quit");
    arguments.append("-I");
    arguments.append(_alignFileNameForPETRegistration);   

    qInfo() << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    FUNC_EXIT;
}

void MainWindow::finishedFMReslicePET(int exitCode, QProcess::ExitStatus exitStatus )
{
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    enablePETActionButtons();
    _centralWidget->setEnabled(true);
    showBrowser(false);
}
void MainWindow::finishedFMAlignPET(int exitCode, QProcess::ExitStatus exitStatus )
{
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    enablePETActionButtons();
    _centralWidget->setEnabled(true);
    showBrowser(false);
}

void MainWindow::analyzeTAC()
{
    FUNC_ENTER;

    // Create a analysis directory on pet ("srtm") and install files
    installSRTMAnalysis();

    auto *process = new QProcess;
    _centralWidget->setEnabled(false);
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

    qInfo() << _fastmapProcess << arguments;
    process->start(_fastmapProcess,arguments);

    FUNC_EXIT;
}
void MainWindow::finishedFMAnalyzeTAC(int exitCode, QProcess::ExitStatus exitStatus )
{
    FUNC_INFO << "exit code" << exitCode << "exit status" << exitStatus;
    enablePETActionButtons();
    _centralWidget->setEnabled(true);
    showBrowser(false);
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

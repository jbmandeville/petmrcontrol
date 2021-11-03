#include <QtWidgets>
#include "petmrMain.h"
#include "ImageIO.h"

void MainWindow::createPETPage()
{
    _petPage = new QWidget();

    auto *petRunLabel = new QLabel("PET run");
    _petRunBox = new QComboBox();
    petFramesBox = new QListWidget();
    petFramesBox->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::MinimumExpanding);
    petFramesBox->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    connect(petFramesBox, SIGNAL(itemClicked(QListWidgetItem*)),this, SLOT(changedPETFrameSelection(QListWidgetItem*)));
    connect(_petRunBox, SIGNAL(activated(int)),this, SLOT(updatePETRunBox(int)));

    auto *petRunLayout = new QHBoxLayout();
    petRunLayout->addWidget(petRunLabel);
    petRunLayout->addWidget(_petRunBox);

    auto *runLayout = new QVBoxLayout();
    runLayout->addLayout(petRunLayout);
    runLayout->addWidget(petFramesBox);

    auto *templateLabel = new QLabel("fMRI template directory (with mcTemplate.nii)");
    auto *fileLabel     = new QLabel("file name(s)");
    _fMRIForPETTemplate = new QLabel("");
    _fMRIForPETFileName = new QLabel("");

    auto *fileLayout = new QGridLayout();
    fileLayout->addWidget(templateLabel,0,0);
    fileLayout->addWidget(_fMRIForPETTemplate,0,1);
    fileLayout->addWidget(fileLabel,1,0);
    fileLayout->addWidget(_fMRIForPETFileName,1,1);

    auto *setupLayout = new QVBoxLayout();
    setupLayout->addLayout(runLayout);
    setupLayout->addLayout(fileLayout);

    auto *runsBox = new QGroupBox("PET run to pre-process");
    runsBox->setLayout(setupLayout);

    _motionCorrectMatchingMRIButton = new QPushButton("Create/motion-correct matching MRI");
    _motionCorrectPETButton = new QPushButton("Apply motion-correction to PET");
    connect(_motionCorrectMatchingMRIButton, SIGNAL(pressed()), this, SLOT(motionCorrectMatchingMRI()));
//    connect(_motionCorrectPETButton,         SIGNAL(pressed()), this, SLOT(motionCorrectPET()));

    auto *actionLayout = new QVBoxLayout();
    actionLayout->addWidget(_motionCorrectMatchingMRIButton);
    actionLayout->addWidget(_motionCorrectPETButton);

    auto *actionBox = new QGroupBox("Go do stuff");
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

void MainWindow::openedPETPage()
{
    FUNC_ENTER;

    QDir const petTopDir("./pet");
    if (!petTopDir.exists())
    {
        _petRunBox->clear();
        return;
    }
    QStringList const petFolderList = petTopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    qInfo() << petFolderList;
    _petRunBox->clear();
    for (int jList=0; jList<petFolderList.size(); jList++)
        _petRunBox->addItem(petFolderList.at(jList));
    _petRunBox->setCurrentIndex(_petRunBox->count()-1);

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
        QString fullFileName = fMRITopDir.absolutePath() + "/" + fMRIFile.name + "/raw.nii";
        QString tagsFileName = fMRITopDir.absolutePath() + "/" + fMRIFile.name + "/time-tags.txt";
        getDimensions(fullFileName, fMRIFile.dim);
        getTimeTags(tagsFileName,fMRIFile.timeTags,fMRIFile.timeText);
        if ( fMRIFile.dim.t > 0 )
            _fMRIFilesForPETMC.append(fMRIFile);
    }

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
        }
    }

    updatePETRunBox(_petRunBox->currentIndex());
    FUNC_EXIT;
}

void MainWindow::updatePETRunBox(int indexInBox)
{
    _petFile.name = "pet/" + _petRunBox->currentText() + "/raw.nii";
    QFileInfo checkFile(_petFile.name);
    if (checkFile.exists() && checkFile.isFile())
    {
        getDimensions(_petFile.name, _petFile.dim);
        FUNC_INFO << "dim.t" << _petFile.dim.t;
        QString fileName = "pet/" + _petRunBox->currentText() + "/time-tags.txt";
        getTimeTags(fileName,_petFile.timeTags,_petFile.timeText);
        findPETandFMRIOverlap();
        FUNC_INFO << "sizes" << _petFile.dim.t << _petFile.timeText.size() << _matchingEPI.size();
//        petFramesBox->clear();
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
            petFramesBox->addItem(&_petFrameItems[jFrame]);
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
    // First write the jip command file to create the matching MRI volume
    writeJipCommandFileForMatchingMRI();

    auto *process = new QProcess;
    _outputBrowser->setWindowTitle("Motion-correct matching MRI");
    showBrowser(true);
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedMotionCorrectMatchingMRI(int, QProcess::ExitStatus)));
    process->setProcessChannelMode(QProcess::MergedChannels);

    QString exe = _scriptDirectory + "motionCorrectMatchingEPI.csh";
    QStringList arguments;
    // arguments: [EPI MC template dir] [pet dir]
    arguments.append(_fMRITemplateDirectoryBox->currentText());
    arguments.append(_petRunBox->currentText());
    qInfo() <<  exe << arguments;
    process->start(exe,arguments);
    _centralWidget->setEnabled(false);
}
void MainWindow::finishedMotionCorrectMatchingMRI(int exitCode, QProcess::ExitStatus exitStatus)
{
    _centralWidget->setEnabled(true);
    showBrowser(false);
}
void MainWindow::writeJipCommandFileForMatchingMRI()
{
    QString fileName = "pet/" + _petRunBox->currentText() + "/jip-createMRI.com";
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

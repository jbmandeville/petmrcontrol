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

    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(runsBox);
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
    int currentRun;
    if ( _matchingEPI.size() > 0 )
    {
        currentRun = _matchingEPI[iFrame].at(0).x;
        if ( _matchingEPI[iFrame].size() > 0 )
            qInfo() << "Overlap from run" << _fMRIFilesForPETMC.at(currentRun).name;
        else
            qInfo() << "No overlap with EPI: motion-correction parameters will be interpolated for this point";
    }
    for (int jPair=0; jPair<_matchingEPI[iFrame].size(); jPair++)
    {
        if  ( _matchingEPI[iFrame].at(jPair).x != currentRun )
            qInfo() << "Overlap from run" << _fMRIFilesForPETMC.at(currentRun).name;
        qInfo() << "time point" << _matchingEPI[iFrame].at(jPair).y;
        currentRun = _matchingEPI[iFrame].at(jPair).x;
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

    updatePETRunBox();
    FUNC_EXIT;
}

void MainWindow::updatePETRunBox()
{
    _petFile.name = "pet/" + _petRunBox->currentText() + "/raw.nii";
    QFileInfo checkFile(_petFile.name);
    if (checkFile.exists() && checkFile.isFile())
    {
        getDimensions(_petFile.name, _petFile.dim);
        QString fileName = "pet/" + _petRunBox->currentText() + "/time-tags.txt";
        getTimeTags(fileName,_petFile.timeTags,_petFile.timeText);
        findPETandFMRIOverlap();
        petFramesBox->clear();
        _petFrameItems.resize(_petFile.dim.t);
        for (int jFrame=0; jFrame<_petFile.dim.t; jFrame++)
        {
            QString text;
            if ( jFrame < _petFile.timeText.size() )
            {
                int nOverlap = _matchingEPI[jFrame].count();
                text = QString("%1: start %2, overlap with %3 fMRI points").arg(jFrame+1).arg(_petFile.timeText.at(jFrame))
                        .arg(nOverlap);
            }
            else
                text = QString("%1").arg(jFrame+1);
            QListWidgetItem item;
            item.setText(text);
            _petFrameItems[jFrame] = item;
            petFramesBox->addItem(&_petFrameItems[jFrame]);
        }
    }
}
void MainWindow::findPETandFMRIOverlap()
{
    _matchingEPI.resize(_petFile.dim.t);
    for (int jFrame=0; jFrame<_petFile.dim.t; jFrame++)
    {
        double width;
        if ( jFrame < _petFile.dim.t-1 )
            width = _petFile.timeTags.at(jFrame+1) - _petFile.timeTags.at(jFrame);
        else // last bin has same width as 2nd to last bin
            width = _petFile.timeTags.at(jFrame) - _petFile.timeTags.at(jFrame-1);
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
    }
}

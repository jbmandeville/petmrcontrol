#include <QtWidgets>
#include <QFileInfo>
#include "petmrMain.h"
#include "ImageIO.h"

void MainWindow::createCleanPage()
{
    FUNC_ENTER;
    _cleanPage = new QWidget();

    _cleanScanTypesBox = new QListWidget();
    _cleanScanTypesBox->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::MinimumExpanding);
    _cleanScanTypesBox->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    connect(_cleanScanTypesBox, SIGNAL(itemChanged(QListWidgetItem*)),this, SLOT(changedScanTypeCheckBox(QListWidgetItem*)));

    auto *totalSizeLabel = new QLabel("total size of all sub-directories");
    _totalSizeSubDirs = new QLabel();
    auto *sizeLayout = new QHBoxLayout();
    sizeLayout->addWidget(totalSizeLabel);
    sizeLayout->addWidget(_totalSizeSubDirs);

    auto *typeLayout = new QVBoxLayout();
    typeLayout->addWidget(_cleanScanTypesBox);
    typeLayout->addLayout(sizeLayout);
    auto *typeBox = new QGroupBox("Directories by scan type");
    typeBox->setLayout(typeLayout);

    _cleanDICOMs         = new QPushButton("clean DICOMS",_cleanPage);
    _cleanNII_auxilliary = new QPushButton("clean intermediate NIFTIs",_cleanPage);
    _cleanNII_mc         = new QPushButton("clean pre-alignment files (mc.nii)",_cleanPage);
    _cleanNII_raw        = new QPushButton("clean raw.nii",_cleanPage);
    connect(_cleanDICOMs,         SIGNAL(pressed()), this, SLOT(cleanDICOMFiles()));
    connect(_cleanNII_auxilliary, SIGNAL(pressed()), this, SLOT(cleanAuxNIIFiles()));
    connect(_cleanNII_mc, SIGNAL(pressed()), this, SLOT(cleanMCFiles()));
    connect(_cleanNII_raw,        SIGNAL(pressed()), this, SLOT(cleanRawNIIFiles()));
    _cleanDICOMs->setToolTip("No reason to keep DICOMs after extracting time tags.");
    _cleanNII_auxilliary->setToolTip("These intermediate files can be recovered if you have the raw files.");
    _cleanNII_mc->setToolTip("If you clean these, you can't easily change alignment/smoothing without repeating earlier steps.");
    _cleanNII_raw->setToolTip("These can be recovered by down-loading again, but it takes time.");

    auto *actionLayout = new QVBoxLayout();
    actionLayout->addWidget(_cleanDICOMs);
    actionLayout->addWidget(_cleanNII_auxilliary);
    actionLayout->addWidget(_cleanNII_mc);
    actionLayout->addWidget(_cleanNII_raw);
    auto *actionBox = new QGroupBox("Clean directories (remove files)");
    actionBox->setLayout(actionLayout);

    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(typeBox);
    pageLayout->addWidget(actionBox);

    _cleanPage->setLayout(pageLayout);
}

void MainWindow::changedScanTypeCheckBox(QListWidgetItem *item)
{
    int iSelected=-1;
    for ( int jItem=0; jItem<_cleanScanTypeItems.size(); jItem++)
    {
        if ( item == &_cleanScanTypeItems.at(jItem) )
            iSelected = jItem;
    }
    FUNC_INFO << iSelected;
    updateCleaningList();
}

void MainWindow::setupScanTypes()
{
    FUNC_ENTER;
    for (int jType=0; jType<_scans.size(); jType++)
    {
        downloadScan scan = _scans.at(jType);
        QString categoryName = scan.categoryName;
        if ( scan.existsOnDisk )
        {
            bool found = false;
            for (int jItem=0; jItem<_cleanScanTypeItems.count(); jItem++)
                if ( ! _cleanScanTypeItems.at(jItem).text().compare(categoryName) )
                    found = true;
            FUNC_INFO << "check scan" << scan.categoryName << scan.scanNumberNew << found;
            if ( !found )
            {
                QListWidgetItem item;
                item.setText(categoryName);
                item.setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
                item.setCheckState(Qt::Checked);
                item.setHidden(false);
                _cleanScanTypeItems.append(item);
            }
        }
    }
    for (int jItem=0; jItem<_cleanScanTypeItems.count(); jItem++)
        _cleanScanTypesBox->addItem(&_cleanScanTypeItems[jItem]);
    FUNC_INFO << "new sizes" << _cleanScanTypeItems.count() << _cleanScanTypesBox->count();

    FUNC_EXIT;
}

void MainWindow::updateCleaningList()
{
    for (int jType=0; jType<_scans.size(); jType++)
    {
        downloadScan scan = _scans.at(jType);
        QString categoryName = scan.categoryName;
        if ( scan.existsOnDisk )
        {
            _scans[jType].selectedForCleaning = false;
            for (int jItem=0; jItem<_cleanScanTypeItems.count(); jItem++)
                if ( ! _cleanScanTypeItems.at(jItem).text().compare(categoryName) )
                    _scans[jType].selectedForCleaning = _cleanScanTypeItems.at(jItem).checkState();
        }
    }
    findDICOMs(false);
    findAuxFiles(false);
    findMCFiles(false);
    findRawFiles(false);
    findAllFiles();
}

void MainWindow::openedCleanPage()
{
    FUNC_ENTER;
    setupScanTypes();
    updateCleaningList();

    FUNC_EXIT;
}

void MainWindow::findDICOMs(bool remove)
{
    qint64 totalSizeAll=0;

    for (int jType=0; jType<_scans.size(); jType++)
    {
        downloadScan scan = _scans.at(jType);
        qint64 totalSizeThisCategory=0;
        if ( scan.existsOnDisk && scan.selectedForCleaning )
        {
            QString spec;
            if ( scan.category != category_PET )
                spec = "MR.*";
            else
                spec = "PT.*";
            QString dirName = scan.categoryName + "/" + scan.scanNumberNew;
            QDir dir(dirName);
            QStringList const fileList = dir.entryList( {spec}, QDir::Files | QDir::NoSymLinks);

            qint64 totalSizeThisDir=0;
            for (int jFile=0; jFile<fileList.size(); jFile++)
            {
                QString fileName = dirName + "/" + fileList.at(jFile);
                QFileInfo checkFile(fileName);
                if ( remove )
                {
                    QFile file(fileName);
                    bool success = file.remove();
                    if ( !success )
                        qInfo() << "failed to remove file" << fileName;
                    else
                        totalSizeThisDir += checkFile.size();
                }
                else
                    totalSizeThisDir += checkFile.size();
            }
            totalSizeThisCategory += totalSizeThisDir;
            if ( jType == _scans.size()-1 || (scan.category != _scans.at(jType+1).category) )
                FUNC_INFO << "total size of DICOMS in" << scan.categoryName << "=" << totalSizeThisCategory/(1024*1024) << "Mb";
        }
        totalSizeAll += totalSizeThisCategory;
    }
    QString gigabytes; gigabytes.setNum(static_cast<double>(totalSizeAll)/(1024.*1024.*1024.),'g',3);
    totalSizeAll /= 1024*1024;  // b -> Mb
    if ( remove )
        qInfo() << "Removed" << gigabytes << "Gb";
    else
        _cleanDICOMs->setText(QString("clean DICOMS: disk space = %1 Gb").arg(gigabytes));
    _cleanDICOMs->setEnabled(totalSizeAll > 0);

    FUNC_EXIT;
}

void MainWindow::findAuxFiles(bool remove)
{
    qint64 totalSizeAll=0;

    for (int jType=0; jType<_scans.size(); jType++)
    {
        downloadScan scan = _scans.at(jType);
        qint64 totalSizeThisCategory=0;
        if ( scan.existsOnDisk && scan.selectedForCleaning )
        {
            QStringList spec;
            spec.append("MR*.nii");
            if ( scan.category == category_PET )
            {
                spec.append("matchingMRI-rs.nii");
                spec.append("matchingMRI-mc.nii");
                spec.append("test.nii");
                spec.append("mc.nii");      // mc --> reslice
            }
            else
                spec.append("reslice.nii"); // reslice --> mc

            QString dirName = scan.categoryName + "/" + scan.scanNumberNew;
            QDir dir(dirName);
//            QStringList const fileList = dir.entryList(spec, QDir::Files | QDir::NoSymLinks);
            QStringList const fileList = dir.entryList(spec, QDir::Files);

            qint64 totalSizeThisDir=0;
            for (int jFile=0; jFile<fileList.size(); jFile++)
            {
                QString fileName = dirName + "/" + fileList.at(jFile);
                QFileInfo checkFile(fileName);
                if ( remove )
                {
                    QFile file(fileName);
                    bool success = file.remove();
                    if ( !success )
                        qInfo() << "failed to remove file" << fileName;
                    else
                        totalSizeThisDir += checkFile.size();
                }
                else
                    totalSizeThisDir += checkFile.size();
            }
            totalSizeThisCategory += totalSizeThisDir;
            if ( jType == _scans.size()-1 || (scan.category != _scans.at(jType+1).category) )
                FUNC_INFO << "total size of DICOMS in" << scan.categoryName << "=" << totalSizeThisCategory/(1024*1024) << "Mb";
        }
        totalSizeAll += totalSizeThisCategory;
    }
    QString gigabytes; gigabytes.setNum(static_cast<double>(totalSizeAll)/(1024.*1024.*1024.),'g',3);
    totalSizeAll /= 1024*1024;  // b -> Mb
    if ( remove )
        qInfo() << "Removed" << gigabytes << "Gb";
    else
        _cleanNII_auxilliary->setText(QString("clean intermediate NIFTIs: disk space = %1 Gb").arg(gigabytes));
    _cleanNII_auxilliary->setEnabled(totalSizeAll > 0);
    FUNC_EXIT;
}

void MainWindow::findMCFiles(bool remove)
{
    qint64 totalSizeAll=0;

    for (int jType=0; jType<_scans.size(); jType++)
    {
        downloadScan scan = _scans.at(jType);
        qint64 totalSizeThisCategory=0;
        if ( scan.existsOnDisk && scan.selectedForCleaning )
        {
            QStringList spec;
            spec.append("mc.nii");

            QString dirName = scan.categoryName + "/" + scan.scanNumberNew;
            QDir dir(dirName);
            QStringList const fileList = dir.entryList(spec, QDir::Files | QDir::NoSymLinks);

            qint64 totalSizeThisDir=0;
            for (int jFile=0; jFile<fileList.size(); jFile++)
            {
                QString fileName = dirName + "/" + fileList.at(jFile);
                QFileInfo checkFile(fileName);
                if ( remove )
                {
                    QFile file(fileName);
                    bool success = file.remove();
                    if ( !success )
                        qInfo() << "failed to remove file" << fileName;
                    else
                        totalSizeThisDir += checkFile.size();
                }
                else
                    totalSizeThisDir += checkFile.size();
            }
            totalSizeThisCategory += totalSizeThisDir;
            if ( jType == _scans.size()-1 || (scan.category != _scans.at(jType+1).category) )
                FUNC_INFO << "total size of DICOMS in" << scan.categoryName << "=" << totalSizeThisCategory/(1024*1024) << "Mb";
        }
        totalSizeAll += totalSizeThisCategory;
    }
    QString gigabytes; gigabytes.setNum(static_cast<double>(totalSizeAll)/(1024.*1024.*1024.),'g',3);
    totalSizeAll /= 1024*1024;  // b -> Mb
    if ( remove )
        qInfo() << "Removed" << gigabytes << "Gb";
    else
        _cleanNII_mc->setText(QString("clean pre-alignment files (usually mc.nii): disk space = %1 Gb").arg(gigabytes));
    _cleanNII_mc->setEnabled(totalSizeAll > 0);
    FUNC_EXIT;
}

void MainWindow::findRawFiles(bool remove)
{
    qint64 totalSizeAll=0;

    for (int jType=0; jType<_scans.size(); jType++)
    {
        downloadScan scan = _scans.at(jType);
        qint64 totalSizeThisCategory=0;
        if ( scan.existsOnDisk && scan.selectedForCleaning )
        {
            QStringList spec;
            spec.append("raw.nii");

            QString dirName = scan.categoryName + "/" + scan.scanNumberNew;
            QDir dir(dirName);
            QStringList const fileList = dir.entryList(spec, QDir::Files | QDir::NoSymLinks);

            qint64 totalSizeThisDir=0;
            for (int jFile=0; jFile<fileList.size(); jFile++)
            {
                QString fileName = dirName + "/" + fileList.at(jFile);
                QFileInfo checkFile(fileName);
                if ( remove )
                {
                    QFile file(fileName);
                    bool success = file.remove();
                    if ( !success )
                        qInfo() << "failed to remove file" << fileName;
                    else
                        totalSizeThisDir += checkFile.size();
                }
                else
                    totalSizeThisDir += checkFile.size();
            }
            totalSizeThisCategory += totalSizeThisDir;
            if ( jType == _scans.size()-1 || (scan.category != _scans.at(jType+1).category) )
                FUNC_INFO << "total size of DICOMS in" << scan.categoryName << "=" << totalSizeThisCategory/(1024*1024) << "Mb";
        }
        totalSizeAll += totalSizeThisCategory;
    }
    QString gigabytes; gigabytes.setNum(static_cast<double>(totalSizeAll)/(1024.*1024.*1024.),'g',3);
    totalSizeAll /= 1024*1024;  // b -> Mb
    if ( remove )
        qInfo() << "Removed" << gigabytes << "Gb";
    else
        _cleanNII_raw->setText(QString("clean raw.nii: disk space = %1 Gb").arg(gigabytes));
    _cleanNII_raw->setEnabled(false);
//    _cleanNII_raw->setEnabled(totalSizeAll > 0);
    FUNC_EXIT;
}

void MainWindow::findAllFiles()
{
    FUNC_ENTER;
    QDirIterator it(".", QDirIterator::Subdirectories);
    qint64 total = 0;
    while (it.hasNext())
    {
        FUNC_INFO << "directory" << it.fileName();
        total += it.fileInfo().size();
        it.next();
    }
    QString gigabytes; gigabytes.setNum(static_cast<double>(total)/(1024.*1024.*1024.),'g',3);
    total /= 1024*1024;
    _totalSizeSubDirs->setText(gigabytes + " Gb");
    FUNC_INFO << "gb" << gigabytes;
    FUNC_EXIT << total << "Mb";
}

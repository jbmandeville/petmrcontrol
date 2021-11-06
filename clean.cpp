#include <QtWidgets>
#include <QFileInfo>
#include "petmrMain.h"
#include "ImageIO.h"

void MainWindow::createCleanPage()
{
    _cleanPage = new QWidget();

    auto *pageLayout = new QVBoxLayout();
//    pageLayout->addWidget(runsBox);
//    pageLayout->addWidget(actionBox);
    _petPage->setLayout(pageLayout);
}

void MainWindow::refreshCleanPage()
{
    QDir const fMRITopDir("./epi");
    if (!fMRITopDir.exists())
        return;
    FUNC_INFO << 1;
    QStringList const epiFolderList = fMRITopDir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    FUNC_INFO << "epi folder list" << epiFolderList;

    for (int jDir=0; jDir<epiFolderList.size(); jDir++)
    {
        QString dirName = "epi/" + epiFolderList.at(jDir);
        QDir subDir = dirName;
//        QStringList const niiList = subDir.entryList( {"*.nii"}, QDir::Files);
//        QStringList const niiList = subDir.entryList( {"*.nii"}, QDir::NoSymLinks);
        QStringList const niiList = subDir.entryList( {"*.nii"}, QDir::Files | QDir::NoSymLinks);
        FUNC_INFO << niiList;
        for (int jFile=0; jFile<niiList.size(); jFile++)
        {
            QString fileName = dirName + "/" + niiList.at(jFile);
            QFileInfo checkFile(fileName);
            FUNC_INFO << "size of" << fileName << "is" << checkFile.size()/(1024*1024) << "Mb";
        }
    }

    findDICOMs();
}

void MainWindow::findDICOMs()
{
    qint64 totalSizeAllDICOMs=0;
    for (int jType=0; jType<_scans.size(); jType++)
    {
        downloadScan scan = _scans.at(jType);
        qint64 totalSizeThisCategory=0;
        if ( scan.existsOnDisk )
        {
            QString spec;
            if ( scan.category != category_PET )
                spec = "MR.*";
            else
                spec = "PT.*";
            QString dirName = scan.categoryName + "/" + scan.scanNumberNew;
            QDir dir(dirName);
            //        FUNC_INFO << "check dirName" << dirName;
            QStringList const dicomList = dir.entryList( {spec}, QDir::Files | QDir::NoSymLinks);

            //        FUNC_INFO << "dicomList" << dicomList;
            qint64 totalSizeThisDir=0;
            for (int jFile=0; jFile<dicomList.size(); jFile++)
            {
                QString fileName = dirName + "/" + dicomList.at(jFile);
                QFileInfo checkFile(fileName);
                totalSizeThisDir += checkFile.size();
                //            FUNC_INFO << "fileName" << fileName << "has size" << checkFile.size();
            }
            totalSizeThisCategory += totalSizeThisDir;
            if ( jType == _scans.size()-1 || (scan.category != _scans.at(jType+1).category) )
                FUNC_INFO << "total size of DICOMS in" << scan.categoryName << "=" << totalSizeThisCategory/(1024*1024) << "Mb";
        }
        totalSizeAllDICOMs += totalSizeThisCategory;
    }
    FUNC_INFO << "total size of all DICOMS" << totalSizeAllDICOMs/(1024*1024) << "Mb";
}

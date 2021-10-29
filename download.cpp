#include <QtWidgets>
#include "petmrMain.h"

void MainWindow::createDownloadPage()
{
    _downLoadPage = new QWidget();

    auto *queryLayout = new QGridLayout();
    auto *subjectIDLabel     = new QLabel("Subject ID",_downLoadPage);
    _subjectIDDownload   = new QLineEdit("?");
    queryLayout->addWidget(subjectIDLabel,0,0);
    queryLayout->addWidget(_subjectIDDownload,0,1);

    auto *downloadIDLabel    = new QLabel("Download ID",_downLoadPage);
    _downloadIDBox  = new QComboBox();
    queryLayout->addWidget(downloadIDLabel,1,0);
    queryLayout->addWidget(_downloadIDBox,1,1);
    connect(_downloadIDBox, SIGNAL(activated(int)), this, SLOT(changedDownloadIDBox(int)));

    auto *downloadPathLabel  = new QLabel("Download Path",_downLoadPage);
    _downloadPathBox = new QComboBox();
    queryLayout->addWidget(downloadPathLabel,2,0);
    queryLayout->addWidget(_downloadPathBox,2,1);
    connect(_downloadPathBox, SIGNAL(activated(int)), this, SLOT(changedDownloadPathBox(int)));

    auto *queryDatabaseLabel = new QLabel("Query database",_downLoadPage);
    _querySubjectButton = new QPushButton("Query",_downLoadPage);
    connect(_querySubjectButton, SIGNAL(pressed()), this, SLOT(queryDownloadPaths()));
    queryLayout->addWidget(queryDatabaseLabel,3,0);
    queryLayout->addWidget(_querySubjectButton,3,1);

    _queryDownloadGroupBox = new QGroupBox("Query database and define download path");
    _queryDownloadGroupBox->setLayout(queryLayout);

    _generateScanListButton = new QPushButton("Generate scan list",_downLoadPage);
    connect(_generateScanListButton, SIGNAL(pressed()), this, SLOT(generateScanList()));
    _generateScanListButton->setEnabled(false);

    _scanItemsBox = new QListWidget();
    _scanItemsBox->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::MinimumExpanding);
    _scanItemsBox->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    connect(_scanItemsBox, SIGNAL(itemClicked(QListWidgetItem*)),
            this, SLOT(changedDownloadScanCheckBox(QListWidgetItem*)));

    auto *runsLayout = new QVBoxLayout();
    runsLayout->addWidget(_generateScanListButton);
    runsLayout->addWidget(_scanItemsBox);
    auto *runsBox = new QGroupBox("Download data from server");
    runsBox->setLayout(runsLayout);

    _downloadDataButton = new QPushButton("Download",_downLoadPage);
    connect(_downloadDataButton, SIGNAL(pressed()), this, SLOT(downloadData()));

    auto *downloadLayout = new QVBoxLayout();
    downloadLayout->addWidget(_downloadDataButton);
    auto *downloadBox = new QGroupBox("Download data from server");
    downloadBox->setLayout(downloadLayout);

    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(_queryDownloadGroupBox);
    pageLayout->addWidget(runsBox);
    pageLayout->addWidget(downloadBox);
    _downLoadPage->setLayout(pageLayout);

    QString unpackLog = "unpack.log";
    QFileInfo checkUnpackLog(unpackLog);
    if (checkUnpackLog.exists() && checkUnpackLog.isFile())
        readUnpackLog();

    QString runLog = "scan-list.log";
    QFileInfo checkRunLog(runLog);
    if (checkRunLog.exists() && checkRunLog.isFile())
        readAvailableScanList();
    _downloadDataButton->setEnabled(enableDownloadData());
}

bool MainWindow::enableDownloadData()
{
    bool enable = false;
    for (int jList=0; jList<_scanItems.count(); jList++)
        enable |= _scanItems[jList].checkState();
    return enable;
}

void MainWindow::queryDownloadPaths()
{ // This function uses "findsession" to locate a subject and data path
    QStringList arguments;
    arguments.append(_subjectIDDownload->text());
    QProcess process;
    process.startDetached(_findsessionProcess,arguments);

    process.start(_findsessionProcess,arguments);
    process.waitForFinished(10000);
    QString output = process.readAllStandardOutput();
    QStringList completeList = output.split("\n");
    int iSubjectLine=1;
    _downloadIDBox->clear();
    _downloadPathBox->clear();
    while ( iSubjectLine<completeList.size() )
    {
        QString subjectLine = completeList.at(iSubjectLine);
        QString pathLine    = completeList.at(iSubjectLine+6);
        QStringList list = subjectLine.split(QRegularExpression("\\s+"));
        QString subject  = list.last();
        list = pathLine.split(QRegularExpression("\\s+"));
        QString path  = list.last();
        _downloadIDBox->addItem(subject);
        _downloadPathBox->addItem(path);
        iSubjectLine += 8;
    }

    _downloadIDBox->setCurrentIndex(_downloadIDBox->count()-1);
    _downloadPathBox->setCurrentIndex(_downloadPathBox->count()-1);

    _generateScanListButton->setEnabled(true);
    _downloadDataButton->setEnabled(enableDownloadData());
}

void MainWindow::generateScanList()
{ // given a subject and datapath, this will download everything using "unpacksdcmdir"
    //    unpacksdcmdir -src $DataPath -targ . -unpackerr -scanonly scan-list.log
    auto *process = new QProcess;
    _outputBrowser->setWindowTitle("Query Progress");
    _outputBrowser->show();
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedGeneratingScanList(int, QProcess::ExitStatus)));

    process->setProcessChannelMode(QProcess::MergedChannels);

    QStringList arguments;
    QString exe = _scriptDirectory + "generateScanList.csh";
    arguments.append(_downloadPathBox->currentText());
    qInfo() <<  exe << arguments;
    process->start(exe,arguments);
    _centralWidget->setEnabled(false);
}

void MainWindow::finishedGeneratingScanList(int exitCode, QProcess::ExitStatus exitStatus)
{
    QString runLog = "scan-list.log";
    QFileInfo checkFile(runLog);
    if (checkFile.exists() && checkFile.isFile())
        readAvailableScanList();
    _centralWidget->setEnabled(true);
    _outputBrowser->hide();
}

void MainWindow::readAvailableScanList()
{
    // Read the time model file
    QFile infile("scan-list.log");
    if (!infile.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    QTextStream in_stream(&infile);

    int iLine=0;
    _scans.clear();
    while (!in_stream.atEnd())
    {
        iLine++;
        downloadScan scan;
        QString line = in_stream.readLine();
        QRegExp rx("[,\\s]");// match a comma or a space
        QStringList stringList = line.split(rx, QString::SkipEmptyParts);
//        QStringList stringList = line.split(QRegularExpression("\\s+"));
        int scanNumber = stringList.at(0).toInt();
        if ( scanNumber < 10 )
            scan.scanNumber = "00" + QString("%1").arg(scanNumber);
        else if ( scanNumber < 100 )
            scan.scanNumber = "0" + QString("%1").arg(scanNumber);
        else
            scan.scanNumber = QString("%1").arg(scanNumber);
        scan.sequenceName = stringList.at(1);
        scan.dim.x = stringList.at(3).toInt();
        scan.dim.y = stringList.at(4).toInt();
        scan.dim.z = stringList.at(5).toInt();
        scan.dim.t = stringList.at(6).toInt();
        if ( scan.sequenceName.contains("local",Qt::CaseInsensitive) ||
             scan.sequenceName.contains("scout",Qt::CaseInsensitive) )
        {
            scan.category = category_scout;
            scan.categoryName = "scout";
        }
        else if ( scan.sequenceName.contains("MPRAGE",Qt::CaseInsensitive) ||
                  scan.sequenceName.contains("MEMPR",Qt::CaseInsensitive) )
        {
            scan.category = category_T1;
            scan.categoryName = "t1";
        }
        else if ( scan.sequenceName.contains("epi",Qt::CaseInsensitive)  ||
                  scan.sequenceName.contains("ep2d",Qt::CaseInsensitive) ||
                  scan.sequenceName.contains("BOLD",Qt::CaseInsensitive) ||
                  scan.sequenceName.contains("IRON",Qt::CaseInsensitive) )
        {
            scan.category = category_EPI;
            scan.categoryName = "epi";
        }
        else if ( scan.sequenceName.contains("MRAC",Qt::CaseInsensitive) )
        {
            scan.category = category_MRAC;
            scan.categoryName = "mrac";
        }
        else if ( scan.sequenceName.contains("VIBE",Qt::CaseInsensitive) )
        {
            scan.category = category_VIBE;
            scan.categoryName = "vibe";
        }
        else if ( scan.sequenceName.contains("UTE",Qt::CaseInsensitive) )
        {
            scan.category = category_UTE;
            scan.categoryName = "ute";
        }
        else
        {
            if ( scan.dim.x == 344 && scan.dim.y == 344 && scan.dim.z == 127 )
            {
                scan.category = category_PET;
                scan.categoryName = "pet";
            }
            else
            {
                scan.category = category_UNKNOWN;
                scan.categoryName = "unknown";
            }
        }

        FUNC_INFO << "scan" << scan.sequenceName << "dim" << scan.dim.x << scan.dim.y << scan.dim.z << scan.dim.t;
        if ( scan.dim.x > 0 && scan.dim.y > 0 && scan.dim.z > 0 && scan.dim.t > 0 )
        {
            QString dirname = scan.categoryName + "/" + scan.scanNumber;
            FUNC_INFO << "check dirname" << dirname;
            QFileInfo checkDir(dirname);
            if (checkDir.exists() && checkDir.isDir())
                scan.selectedForDownload = false;
            else
                scan.selectedForDownload = true;
            _scans.append(scan);
        }
    } // while !end
    infile.close();

    // add scans to the table
    _scanItems.resize(_scans.size());
    for (int jList=0; jList<_scans.size(); jList++)
    {
        downloadScan scan = _scans.at(jList);
        QString text = QString("%1 %2 (%3) %4x%5x%6x%7").arg(scan.scanNumber).arg(scan.sequenceName).arg(scan.categoryName)
                .arg(scan.dim.x).arg(scan.dim.y).arg(scan.dim.z).arg(scan.dim.t);
        _scanItems[jList].setText(text);
        _scanItems[jList].setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
        if ( scan.selectedForDownload )
            _scanItems[jList].setCheckState(Qt::Checked);
        else
            _scanItems[jList].setCheckState(Qt::Unchecked);
        _scanItems[jList].setHidden(false);
        _scanItemsBox->addItem(&_scanItems[jList]);
    }
}

void MainWindow::changedDownloadScanCheckBox(QListWidgetItem *item)
{
    FUNC_ENTER;
    int iSelected=-1;
    for ( int jItem=0; jItem<_scanItems.size(); jItem++)
    {
        if ( item == &_scanItems.at(jItem) )
            iSelected = jItem;
    }
    _scans[iSelected].selectedForDownload = _scanItems[iSelected].checkState();
}

void MainWindow::readUnpackLog()
{
    FUNC_ENTER;
    // Read the time model file
    QFile infile("unpack.log");
    if (!infile.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    QTextStream in_stream(&infile);

    int iLine=0;
    while (!in_stream.atEnd())
    {
        iLine++;
        downloadScan scan;
        QString line = in_stream.readLine();
        QStringList stringList = line.split(QRegularExpression("\\s+"));
        FUNC_INFO << stringList << stringList;
        int iPatient=stringList.indexOf("PatientName");
        if ( iPatient >= 0 && iPatient < stringList.size()-1 )
        {
            FUNC_INFO << "add subject" << stringList.at(iPatient+1);
            _downloadIDBox->addItem(stringList.at(iPatient+1));
            _downloadIDBox->setEnabled(false);  // lock the ID box
            _querySubjectButton->setEnabled(false);
        }
        int iSrc=stringList.indexOf("-src");
        FUNC_INFO << stringList << "iSrc" << iSrc;
        if ( iSrc >= 0 && iSrc < stringList.size()-1 )
        {
            FUNC_INFO << "add path" << iSrc+1 << "size" << stringList.size();
            _downloadPathBox->addItem(stringList.at(iSrc+1));
            _downloadPathBox->setEnabled(false);
            _querySubjectButton->setEnabled(false);
            _generateScanListButton->setEnabled(true);
        }
    }
}

void MainWindow::outputConfigurationFile()
{
    QFile file("download-list.dat");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);

    for (int jList=0; jList<_scans.size(); jList++)
    {
        downloadScan scan = _scans.at(jList);
        if ( scan.selectedForDownload )
        {
            if ( scan.category == category_scout )
                out << scan.scanNumber << " scout dicom raw.nii \n";
            else if ( scan.category == category_T1 )
                out << scan.scanNumber << " t1 dicom raw.nii \n";
            else if ( scan.category == category_EPI )
                out << scan.scanNumber << " epi dicom raw.nii \n";
            else if ( scan.category == category_MRAC )
                out << scan.scanNumber << " mrac dicom raw.nii \n";
            else if ( scan.category == category_UTE )
                out << scan.scanNumber << " ute dicom raw.nii \n";
            else if ( scan.category == category_VIBE )
                out << scan.scanNumber << " vibe dicom raw.nii \n";
            else if ( scan.category == category_PET )
                out << scan.scanNumber << " pet dicom raw.nii \n";
            else if ( scan.category == category_UNKNOWN )
                out << scan.scanNumber << " unknown dicom raw.nii \n";
        }
        FUNC_INFO << "scan " << jList << "number" << _scans.at(jList).scanNumber << "name" << _scans.at(jList).sequenceName;
    }
    file.close();
}

void MainWindow::downloadData()
{ // given a subject and datapath, this will download everything using "unpacksdcmdir"
    //    unpacksdcmdir -src $DataPath -targ . -unpackerr -scanonly scan-list.log
    outputConfigurationFile();

    auto *process = new QProcess;
    _outputBrowser->setWindowTitle("Download Progress");
    _outputBrowser->show();
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(enableGUI(int, QProcess::ExitStatus)));

    QStringList arguments;
    QString exe = _scriptDirectory + "unpackData.csh";
    arguments.append(_downloadPathBox->currentText());
    qInfo() <<  exe << arguments;
    process->start(exe,arguments);
    _centralWidget->setEnabled(false);
}

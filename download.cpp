#include <QtWidgets>
#include <QFile>
#include "petmrMain.h"

void MainWindow::createDownloadPage()
{
    _downLoadPage = new QWidget();

    auto *queryLayout = new QGridLayout();
    auto *subjectIDLabel     = new QLabel("Download search string",_downLoadPage);
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
    _querySubjectButton->setToolTip("Find all download paths that contain the given download ID.");

    _queryDownloadGroupBox = new QGroupBox("Query database and define download path");
    _queryDownloadGroupBox->setLayout(queryLayout);

    _generateScanListButton = new QPushButton("Generate scan list",_downLoadPage);
    connect(_generateScanListButton, SIGNAL(pressed()), this, SLOT(generateScanList()));
    _generateScanListButton->setEnabled(false);
    _generateScanListButton->setToolTip("Generate a list of available runs from the given download path.");

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
    _downloadDataButton->setToolTip("Download all scans that are checked in the box above.");

    auto *downloadLayout = new QVBoxLayout();
    downloadLayout->addWidget(_downloadDataButton);
    auto *downloadBox = new QGroupBox("Download data from server");
    downloadBox->setLayout(downloadLayout);

    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(_queryDownloadGroupBox);
    pageLayout->addWidget(runsBox);
    pageLayout->addWidget(downloadBox);
    _downLoadPage->setLayout(pageLayout);

    QString downFile = "download-list.dat";
    QFileInfo checkDownFile(downFile);
    if (checkDownFile.exists() && checkDownFile.isFile())
        _downloadDataButton->setStyleSheet("background-color:lightYellow;");
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
{
    auto *process = new QProcess;
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedGeneratingScanList(int, QProcess::ExitStatus)));

    process->setProcessChannelMode(QProcess::MergedChannels);

    QStringList arguments;
    QString exe = _scriptDirectory + "generateScanList.csh";
    arguments.append(_downloadPathBox->currentText());
    QString message = "Create scan list (scan-list.log); this takes about 15 minutes";
    spawnProcess(process,exe,arguments,message,"Generate scan list");
}

void MainWindow::finishedGeneratingScanList(int exitCode, QProcess::ExitStatus exitStatus)
{
    readAvailableScanList();
    finishedProcess();
    _downloadDataButton->setEnabled(enableDownloadData());
    qInfo() << "finished: generating scan list";
}

void MainWindow::readAvailableScanList()
{
    QString runLog = "scan-list.log";
    QFileInfo checkRunLog(runLog);
    if (! (checkRunLog.exists() && checkRunLog.isFile()) )
        return;

    // Read the time model file
    QFile infile(runLog);
    if (!infile.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    _generateScanListButton->setStyleSheet("background-color:lightYellow;");
    _showNotesAction->setChecked(true);
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
        scan.scanNumber = QString("%1").arg(scanNumber);
        if ( scanNumber < 10 )
            scan.scanNumberNew = "00" + QString("%1").arg(scanNumber);
        else if ( scanNumber < 100 )
            scan.scanNumberNew = "0" + QString("%1").arg(scanNumber);
        else
            scan.scanNumberNew = QString("%1").arg(scanNumber);
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
            QString dirname = scan.categoryName + "/" + scan.scanNumberNew;
            FUNC_INFO << "check dirname" << dirname;
            QFileInfo checkDir(dirname);
            if (checkDir.exists() && checkDir.isDir())
            {
                FUNC_INFO << "exists on disk";
                scan.existsOnDisk        = true;
                scan.selectedForDownload = false;
            }
            else
            {
                FUNC_INFO << "does not exist on disk";
                scan.existsOnDisk        = false;
                scan.selectedForDownload = scan.category != category_UNKNOWN;
            }
            _scans.append(scan);
        }
    } // while !end
    infile.close();

    // add scans to the table
    _scanItems.resize(_scans.size());
    for (int jList=0; jList<_scans.size(); jList++)
    {
        downloadScan scan = _scans.at(jList);
        QString volumes = "volumes";  if ( scan.dim.t == 1 ) volumes = "volume";
        QString text = QString("%1 (%2) %3 %4x%5x%6, %7 %8").arg(scan.scanNumberNew).arg(scan.categoryName).arg(scan.sequenceName)
                .arg(scan.dim.x).arg(scan.dim.y).arg(scan.dim.z).arg(scan.dim.t).arg(volumes);
        _scanItems[jList].setText(text);
        if ( scan.category == category_UNKNOWN )
            _scanItems[jList].setFlags(Qt::ItemIsUserCheckable);
        else
            _scanItems[jList].setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
        if ( scan.selectedForDownload )
            _scanItems[jList].setCheckState(Qt::Checked);
        else
            _scanItems[jList].setCheckState(Qt::Unchecked);
        if ( scan.existsOnDisk )
            _scanItems[jList].setBackgroundColor(Qt::yellow);
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
    _downloadDataButton->setEnabled(enableDownloadData());
}

void MainWindow::readUnpackLog()
{
    FUNC_ENTER;
    QString unpackLog = "unpack.log";
    QFileInfo checkUnpackLog(unpackLog);
    if ( !(checkUnpackLog.exists() && checkUnpackLog.isFile()) )
        return;

    FUNC_INFO << "read unpack -src";
    QString sourceName = readFileTextArgument("unpack.log", "-src");
    if ( !sourceName.isEmpty() )
    {
        FUNC_INFO << "add path" << sourceName;
        _downloadPathBox->addItem(sourceName);
//        _queryDownloadGroupBox->setEnabled(false);
        _queryDownloadGroupBox->setStyleSheet("background-color:lightYellow;");
        _generateScanListButton->setEnabled(true);
    }

    FUNC_INFO << "read unpack PatientName";
    QString patientName = readFileTextArgument("unpack.log", "PatientName");
    FUNC_INFO << "patientName" << patientName << patientName.isEmpty();
    if ( patientName.isEmpty() )
    {
        QString fileName = "t1/" + _anatomyDirBox->currentText() + "/raw.nii-infodump.dat";
        FUNC_INFO << "read" << fileName;
        patientName = readFileTextArgument(fileName, "PatientName");
    }
    if ( !patientName.isEmpty() )
    {
        FUNC_INFO << "add subject" << patientName;
        _downloadIDBox->addItem(patientName);
//        _queryDownloadGroupBox->setEnabled(false);
        _queryDownloadGroupBox->setStyleSheet("background-color:lightYellow;");
    }
}

void MainWindow::outputDownloadList()
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
        FUNC_INFO << "scan " << jList << "number" << _scans.at(jList).scanNumberNew << "name" << _scans.at(jList).sequenceName;
    }
    file.close();
}

void MainWindow::downloadData()
{
    outputDownloadList();

    auto *process = new QProcess();
    QObject::connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::outputToBrowser);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(finishedDownloadData(int, QProcess::ExitStatus)));

    QStringList arguments;
    QString exe = _scriptDirectory + "unpackData.csh";
    arguments.append(_downloadPathBox->currentText());

    QString message = "Download data from server; this takes about 20 minutes";
    spawnProcess(process,exe,arguments,message,"Download Progress");
}

void MainWindow::finishedDownloadData(int exitCode, QProcess::ExitStatus exitStatus)
{
    qInfo() << "finished: downloading data";

    for (int jList=0; jList<_scans.size(); jList++)
        reformatAcquisitionTimes(_scans.at(jList));

    readAvailableScanList();
    finishedProcess();
}

void MainWindow::finishedDownloadDataNew()
{

    for (int jList=0; jList<_scans.size(); jList++)
        reformatAcquisitionTimes(_scans.at(jList));
}

void MainWindow::reformatAcquisitionTimes(downloadScan scan)
{
    FUNC_ENTER;
    QString inputFileName  = scan.categoryName + "/" + scan.scanNumberNew + "/time-tags.tmp";
    QString outputFileName = scan.categoryName + "/" + scan.scanNumberNew + "/time-tags.txt";

    QFile infile(inputFileName);
    if (!infile.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    QTextStream in_stream(&infile);

    QFile outFile(outputFileName);
    if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&outFile);
    out << "Volume secondsTotal time\n";

    in_stream.readLine();  // "done retrieving file list xxx files found"
    in_stream.readLine();  // "Values:"

    int iTime=0;
    while (!in_stream.atEnd())
    {
        QString line = in_stream.readLine();
        QRegExp rx("[\\s]");// match a comma or a space
        QStringList stringList = line.split(rx, QString::SkipEmptyParts);
        if ( stringList.count() > 0 )
        {
            QString text = stringList.at(0);
            FUNC_INFO << text;
            if ( ! text.compare("Mapping:")) break;
            double value   = stringList.at(0).toDouble();
            short hours    = static_cast<short>(value/10000.);
            short minutes  = static_cast<short>((value - 10000*hours)/100.);
            short seconds = static_cast<short>(value - 10000*hours - 100*minutes);
            double secondsTotal = 3600. * hours + 60. * minutes + seconds;
            FUNC_INFO << "time" << hours << minutes << seconds;
            out << QString("%1 %2 %3:%4.%5\n").arg(iTime).arg(secondsTotal).
                   arg(twoDigits(hours)).arg(twoDigits(minutes)).arg(twoDigits(seconds));
            iTime++;
        }
    }
    infile.close();  outFile.close();
    FUNC_EXIT << iTime;
}

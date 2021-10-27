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
    _downloadID  = new QLabel("?");
    queryLayout->addWidget(downloadIDLabel,1,0);
    queryLayout->addWidget(_downloadID,1,1);

    auto *downloadPathLabel  = new QLabel("Download Path",_downLoadPage);
    _downloadPath = new QLineEdit("?");
    queryLayout->addWidget(downloadPathLabel,2,0);
    queryLayout->addWidget(_downloadPath,2,1);

    auto *queryDatabaseLabel = new QLabel("Query database",_downLoadPage);
    _downloadQueryButton = new QPushButton("Query",_downLoadPage);
    connect(_downloadQueryButton, SIGNAL(pressed()), this, SLOT(queryDownloadPaths()));
    queryLayout->addWidget(queryDatabaseLabel,3,0);
    queryLayout->addWidget(_downloadQueryButton,3,1);

    auto *queryBox = new QGroupBox("Query database and define download path");
    queryBox->setLayout(queryLayout);

    _downloadDataButton = new QPushButton("Download",_downLoadPage);
    connect(_downloadDataButton, SIGNAL(pressed()), this, SLOT(downloadData()));
    _downloadDataButton->setEnabled(false);

    auto *downloadLayout = new QVBoxLayout();
    downloadLayout->addWidget(_downloadDataButton);
    auto *downloadBox = new QGroupBox("Download data from server");
    downloadBox->setLayout(downloadLayout);

    auto *pageLayout = new QVBoxLayout();
    pageLayout->addWidget(queryBox);
    pageLayout->addWidget(downloadBox);
    _downLoadPage->setLayout(pageLayout);
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

    int iEnd = output.size();
    int lastEqual = output.lastIndexOf("=");
    output = output.right(iEnd - lastEqual - 1);

    QStringList fullList = output.split("\n");
    QString subjectPlusLabel = fullList.at(1);
    QStringList list = subjectPlusLabel.split(QRegularExpression("\\s+"));
    QString subject = list.at(1);

    QString pathPlusLabel = fullList.at(fullList.size()-2);
    list = pathPlusLabel.split(QRegularExpression("\\s+"));
    QString path = list.at(2);

    _downloadID->setText(subject);
    _downloadPath->setText(path);
    _downloadDataButton->setEnabled(true);

}

void MainWindow::downloadData()
{ // given a subject and datapath, this will download everything using "unpacksdcmdir"
    //    unpacksdcmdir -src $DataPath -targ . -unpackerr -scanonly run-list.log
    auto *process = new QProcess;
    auto *view = new QTextBrowser;
    QObject::connect(process, &QProcess::readyReadStandardOutput, [process,view]()
    {
        auto output=process->readAllStandardOutput();
        view->append(output);
    });
    _centralWidget->setEnabled(false);
    connect(process, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(enableGUI(int, QProcess::ExitStatus)));

    QStringList arguments;
    arguments.append(_downloadPath->text());
    QString exe = _scriptDirectory + "unpack.csh";
    qInfo() <<  exe << arguments;
    process->start(exe,arguments);

    qInfo() << "unpacksdcmdir" << arguments;
}

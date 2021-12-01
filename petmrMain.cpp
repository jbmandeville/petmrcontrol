#include <QtWidgets>
#include "ImageIO.h"
#include "petmrMain.h"

MainWindow::~MainWindow()
{
}

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
    QCoreApplication::setOrganizationName("Martinos");
    QCoreApplication::setOrganizationDomain("http://www.nmr.mgh.harvard.edu");
    QCoreApplication::setApplicationVersion("3.0");
    QCoreApplication::setApplicationName("fastmap");
    QSettings fmSettings;
    _FastmapMSTemplateDirectories = fmSettings.value("templateDirectories","").toStringList();
    FUNC_INFO << "template directories" << _FastmapMSTemplateDirectories;

    QCoreApplication::setOrganizationName("Martinos");
    QCoreApplication::setOrganizationDomain("http://www.nmr.mgh.harvard.edu");
    QCoreApplication::setApplicationVersion("3.0");
    QCoreApplication::setApplicationName("petmrcontrol");

    readQSettings();

    _radioButtonHumanBay7 = new QRadioButton("Human Bay 7");
    _radioButtonNHPBay6   = new QRadioButton("NHP Bay 6");
    _radioButtonNHPBay7   = new QRadioButton("NHP Bay 7");
    connect(_radioButtonHumanBay7, SIGNAL(clicked(bool)), this, SLOT(dataOriginChanged()));
    connect(_radioButtonNHPBay6,   SIGNAL(clicked(bool)), this, SLOT(dataOriginChanged()));
    connect(_radioButtonNHPBay7,   SIGNAL(clicked(bool)), this, SLOT(dataOriginChanged()));
    _radioButtonHumanBay7->setChecked(true);
    auto *radioLayout = new QHBoxLayout();
    radioLayout->addWidget(_radioButtonHumanBay7);
    radioLayout->addWidget(_radioButtonNHPBay6);
    radioLayout->addWidget(_radioButtonNHPBay7);
    auto *radioBox = new QGroupBox("Data origin");
    radioBox->setLayout(radioLayout);
    radioBox->setStyleSheet("background-color:lightBlue;");

    _tabs = new QTabWidget();
    createDownloadPage();
    createAnatomyPage();
    createfMRIPage();
    createPETPage();
    createCleanPage();
    _tabs->addTab(_downLoadPage,  tr("Download"));
    _tabs->addTab(_anatomyPage, tr("Anatomy"));
    _tabs->addTab(_fmriPage, tr("fMRI"));
    _tabs->addTab(_petPage, tr("PET"));
    _tabs->addTab(_cleanPage, tr("clean up"));
    _tabs->setTabToolTip(0,"Query database and download data");
    _tabs->setTabToolTip(1,"Align anatomy; potentially run freeSurfer;\nClick tab to refresh state");
    connect(_tabs, SIGNAL(tabBarClicked(int)), this, SLOT(changedPage(int)));

    _centralWidget = new QWidget(this);
    this->setCentralWidget( _centralWidget );
    auto *mainLayout = new QVBoxLayout( _centralWidget );
    mainLayout->addWidget(radioBox);
    mainLayout->addWidget(_tabs);

    _statusBar = this->statusBar();
    _statusBar->setStyleSheet("color:Darkred");
    mainLayout->addWidget(_statusBar);

    // add a menu
    auto *menuBar = new QMenuBar;
    mainLayout->setMenuBar(menuBar);
    QMenu *mainMenu  = new QMenu(tr("&Menu"), this);
    auto *helpMenu  = new QMenu(tr("Help"), this);
    menuBar->addMenu(mainMenu);
    menuBar->addMenu(helpMenu);

    QAction *aboutAppAct = helpMenu->addAction(tr("About this"));
    connect(aboutAppAct, &QAction::triggered, this, &MainWindow::aboutApp);

    QAction *quitAction = mainMenu->addAction(tr("&Quit"));
    // short-cuts and tooltips
    quitAction->setShortcut(Qt::ControlModifier + Qt::Key_Q);
    connect(quitAction, &QAction::triggered, this, &MainWindow::exitApp);

    QSize iconSizeSmall(24,24);
    const QIcon *showBrowser = new QIcon(":/My-Icons/textOutput.png");
    QToolBar *sideToolBar = addToolBar(tr("tool bar"));
    _browserAction = new QAction(*showBrowser,"browser",this);
    _browserAction->setCheckable(true);
    _browserAction->setChecked(false);
    sideToolBar->addAction(_browserAction);
    sideToolBar->setIconSize(iconSizeSmall);
    addToolBar(Qt::LeftToolBarArea, sideToolBar);
    connect(_browserAction, SIGNAL(toggled(bool)), this, SLOT(showBrowser(bool)));

    _outputBrowser = new QTextBrowser;
    _browserAction->setChecked(false);

    /*
    QSize defaultWindowSize;
    QRect rec = QApplication::desktop()->screenGeometry();
    defaultWindowSize.setWidth(rec.width()/4);
    defaultWindowSize.setHeight(rec.height()/2);
    resize(defaultWindowSize);
    _outputBrowser->resize(defaultWindowSize);
    */

    openedAnatomyPage();
    readUnpackLog();
    readAvailableScanList();
    readSubjectVariables();

    restoreGeometry(_savedSettings.imageWindowGeometry);
    _outputBrowser->restoreGeometry(_savedSettings.browserWindowGeometry);
}

void MainWindow::readSubjectVariables()
{
    FUNC_ENTER;
    QString fileName = "analyze.dat";
    QFileInfo checkFile(fileName);
    if ( !(checkFile.exists() && checkFile.isFile()) )
        return;

    // Read the time model file
    QFile infile(fileName);
    if (!infile.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    QTextStream in_stream(&infile);

    int iLine=0;
    QString templateID;
    while (!in_stream.atEnd())
    {
        iLine++;
        QString line = in_stream.readLine();
        QStringList stringList = line.split(QRegularExpression("\\s+"));
        if ( stringList.isEmpty() ) continue;
        QString variable = stringList.at(0);
        FUNC_INFO << variable << stringList.count();
        if ( !variable.compare("multi-subject-template") && stringList.count() > 1 )
            templateID = stringList.at(1);
        if ( !variable.compare("analysis-type") && stringList.count() > 1 )
        {
            QString type = stringList.at(1);
            FUNC_INFO << type;
            if ( !type.compare("nhp-bay7") )
                _radioButtonNHPBay7->setChecked(true);
            else if ( !type.compare("nhp-bay6") )
                _radioButtonNHPBay6->setChecked(true);
            else
                _radioButtonHumanBay7->setChecked(true);
        }
    }
    infile.close();

    // find the template directory
    int iSelection = 0;
    FUNC_INFO << "save template" << templateID;
    for (int jList=0; jList<_FastmapMSTemplateDirectories.count(); jList+=2)
    {
        if ( !templateID.compare(_FastmapMSTemplateDirectories.at(jList)) )
            iSelection = jList/2;
    }
    _anatomyTemplateDirectory->setCurrentIndex(iSelection);
    FUNC_EXIT;
}
void MainWindow::writeSubjectVariables()
{
    FUNC_ENTER;
    QString fileName = "analyze.dat";
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&file);

    out << "# analysis variables for this session" << "\n";
    out << "multi-subject-template " << _anatomyTemplateDirectory->currentText() << "\n";
    if ( _radioButtonHumanBay7->isChecked() )
        out << "analysis-type human-bay7\n";
    else if ( _radioButtonNHPBay6->isChecked() )
        out << "analysis-type nhp-bay6\n";
    else
        out << "analysis-type nhp-bay7\n";
    file.close();
}

void MainWindow::readSmoothing(int which)
{ // which: 0=anatomy, 1=fMRI, 2=pet
    QString fileName = alignCOMFileName(which);
    FUNC_ENTER << which << fileName;
    QString argument = readFileTextArgument(fileName, "smoothing");
    if ( !argument.isEmpty() )
    {
        FUNC_INFO << "argument" << argument;
        bool ok;
        double smoothing = argument.toDouble(&ok);
        if ( !ok ) smoothing = 0.;
        QString number; number.setNum(smoothing);
        FUNC_INFO << "number" << number;
        if ( which == 0 )
        {
            _smoothingAnatomy->setText(number);
            _smoothingfMRI->setText(number);
            _smoothingPET->setText(number);
        }
        else if ( which == 1 )
        {
            _smoothingfMRI->setText(number);
            _smoothingPET->setText(number);
        }
        else
            _smoothingPET->setText(number);
    }
}

QString MainWindow::readFileTextArgument(QString fileName, QString parameterName)
{
    QFileInfo checkFile(fileName);
    if ( !(checkFile.exists() && checkFile.isFile()) )
        return "";

    QFile infile(fileName);
    if (!infile.open(QIODevice::ReadOnly | QIODevice::Text))
        return "";
    QTextStream in_stream(&infile);

    int iLine=0;
    while (!in_stream.atEnd())
    {
        iLine++;
        QString line = in_stream.readLine();
        if ( line.isEmpty() ) continue;
        FUNC_INFO << "line" << line;
        QStringList stringList = line.split(QRegularExpression("\\s+"));
        FUNC_INFO << stringList << stringList;
        int iParameter=stringList.indexOf(parameterName);
        FUNC_INFO << "iParameter" << iParameter;
        if ( iParameter >= 0 && iParameter < stringList.size()-1 )
        {
            FUNC_INFO << "read iParameter";
            return stringList.at(iParameter+1);
        }
    }
    infile.close();
    FUNC_EXIT;
    return "";
}

QString MainWindow::alignCOMFileName(int which)
{ // which: 0=anatomy, 1=fMRI, 2=pet
    FUNC_ENTER << which << _fMRITemplateDirBox->count();
    QString fileName = "";
    if ( which == 0 && anatomyFileExists("align.com") )
        fileName = "t1/" + _anatomyDirBox->currentText() + "/align.com";
    else if ( which == 1 && epiFileExists("align.com") )
        fileName = "epi/" + _fMRITemplateDirBox->currentText() + "/align.com";
    else if ( which ==2 && petFileExists("align.com") )
        fileName = "pet/" + _petDirBox->currentText() + "/align.com";
    return fileName;
}

void MainWindow::changedSmoothingAnatomy()
{
    bool ok;
    double smoothing = _smoothingAnatomy->text().toDouble(&ok);
    if ( !ok ) smoothing = 0.;
    QString number; number.setNum(smoothing);
    _smoothingAnatomy->setText(number);
    _smoothingfMRI->setText(number);
    _smoothingPET->setText(number);
}

void MainWindow::changedSmoothingfMRI()
{
    bool ok;
    double smoothing = _smoothingfMRI->text().toDouble(&ok);
    if ( !ok ) smoothing = 0.;
    QString number; number.setNum(smoothing);
    _smoothingAnatomy->setText(number);
    _smoothingfMRI->setText(number);
    _smoothingPET->setText(number);
}

void MainWindow::changedSmoothingPET()
{
    bool ok;
    double smoothing = _smoothingPET->text().toDouble(&ok);
    if ( !ok ) smoothing = 0.;
    QString number; number.setNum(smoothing);
    _smoothingAnatomy->setText(number);
    _smoothingfMRI->setText(number);
    _smoothingPET->setText(number);
}

void MainWindow::dataOriginChanged()
{
    _freeSurferGroupBox->setVisible(_radioButtonHumanBay7->isChecked());
    _extractFreeSurferOverlaysButton->setVisible(_radioButtonHumanBay7->isChecked());
    writeSubjectVariables();
}

void MainWindow::changedPage(int index)
{
    if ( index == page_anatomy )
        openedAnatomyPage();
    else if ( index == page_fMRI )
        openedfMRIPage();
    else if ( index == page_PET )
        openedPETPage();
    else if ( index == page_clean )
        openedCleanPage();
}

void MainWindow::aboutApp()
{
    QMessageBox msgBox;
    QString version = qVersion();
    QString text = "Qt Version " + version;
    msgBox.setText(text);
    msgBox.setIcon(QMessageBox::Information);
    msgBox.exec();

    QMessageBox::information(nullptr, QGuiApplication::applicationDisplayName(),
                             QGuiApplication::applicationDisplayName() + ' '
                             + QCoreApplication::applicationVersion() + " , by Joe Mandeville;\n" +
                             "Request bug fixes by email to\njbm@nmr.mgh.harvard.edu\nwith subject line 'simulator bug'.");
}

void MainWindow::exitApp()
{
    writeQSettings();
    QCoreApplication::exit(0);
}

void MainWindow::readQSettings()
{
    QByteArray defaultImageWindowGeometry = saveGeometry();  // defined in ImageWindow constructor
    _savedSettings.imageWindowGeometry    = _savedQSettings.value("imageWindowGeometry",defaultImageWindowGeometry).toByteArray();
    _savedSettings.browserWindowGeometry  = _savedQSettings.value("browserWindowGeometry",defaultImageWindowGeometry).toByteArray();
}

void MainWindow::writeQSettings()
{
    if ( !isMaximized() )
        _savedQSettings.setValue("imageWindowGeometry",saveGeometry());
    if ( !isMaximized() )
        _savedQSettings.setValue("browserWindowGeometry",_outputBrowser->saveGeometry());

    _savedQSettings.sync();
}


QString MainWindow::getDimensions(QString fileName, iPoint4D &dim)
{
    QFileInfo checkFile(fileName);
    if ( checkFile.exists() && checkFile.isFile() )
    {
        ImageIO file;
        if ( !file.readFileHeader(fileName,true) )
        {
            dim = file.getDimensions();
            dPoint4D res = file.getResolution();
            qInfo() << "time step" << res.t;
            double duration = dim.t * res.t / 60.;
            QString text = QString("%1 x %2 x %3 with %4 time points (%5 min)")
                    .arg(dim.x).arg(dim.y).arg(dim.z).arg(dim.t).arg(duration);
            return text;
        }
        else
        {
            dim={0,0,0,0};
            return "";
        }
    }
    else
    {
        dim={0,0,0,0};
        return "";
    }
}

void MainWindow::getTimeTags(QString fileName, dVector &timeTags, sVector &timeText )
{
    FUNC_ENTER << fileName;
    timeTags.clear();  timeText.clear();

    QFile inFile(fileName);
    if (!inFile.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    QTextStream in_stream(&inFile);
    in_stream.readLine(); // header

    int iLine=0;
    while (!in_stream.atEnd())
    {
        iLine++;
        QString line = in_stream.readLine();
        QRegExp rx("[,\\s]");// match a comma or a space
        QStringList stringList = line.split(rx, QString::SkipEmptyParts);
//        FUNC_INFO << stringList;
        if ( stringList.count() >= 3 )
        {
            double value = stringList.at(1).toDouble();
            timeTags.append(value);
            timeText.append(stringList.at(2));
        }
    }
    inFile.close();
}

QString MainWindow::twoDigits(short time)
{
    QString text;
    if ( time < 10 )
        text = QString("0%1").arg(time);
    else
        text = QString("%1").arg(time);
    return text;
}

dPoint2D MainWindow::petFrameTime(int iFrame)
{
    double startPET = _petFile.timeTags.at(iFrame);
    double endPET;
    if ( iFrame < _petFile.timeTags.size()-1 )
        endPET = _petFile.timeTags.at(iFrame+1);
    else
        endPET = startPET + (_petFile.timeTags.at(iFrame)-_petFile.timeTags.at(iFrame-1));
    return {startPET, endPET};

    /*
    double width;
    if ( iFrame < _petFile.dim.t-1 )
        width = _petFile.timeTags.at(iFrame+1) - _petFile.timeTags.at(iFrame);
    else // last bin has same width as 2nd to last bin
        width = _petFile.timeTags.at(iFrame) - _petFile.timeTags.at(iFrame-1);
    double timePET = _petFile.timeTags.at(iFrame);
    double startPET = timePET - width/2.;
    double endPET   = timePET + width/2.;
    return {startPET, endPET};
    */
}

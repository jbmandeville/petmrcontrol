#include <QtWidgets>
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

    _tabs = new QTabWidget();
    createDownloadPage();
    createAnatomyPage();
    createfMRIPage();
    createPETPage();
    _tabs->addTab(_downLoadPage,  tr("Download"));
    _tabs->addTab(_anatomyPage, tr("Anatomy"));
    _tabs->addTab(_fmriPage, tr("fMRI"));
    _tabs->addTab(_petPage, tr("PET"));
    _tabs->setTabToolTip(0,"Query database and download data");
    _tabs->setTabToolTip(1,"Align anatomy; potentially run freeSurfer;\nClick tab to refresh state");
    connect(_tabs, SIGNAL(tabBarClicked(int)), this, SLOT(changedPage(int)));

    _centralWidget = new QWidget(this);
    this->setCentralWidget( _centralWidget );
    auto *mainLayout = new QVBoxLayout( _centralWidget );
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

    QSize defaultWindowSize;
    QRect rec = QApplication::desktop()->screenGeometry();
    defaultWindowSize.setWidth(rec.width()/4);
    defaultWindowSize.setHeight(rec.height()/4);
    resize(defaultWindowSize);
    _outputBrowser->resize(defaultWindowSize);
}

void MainWindow::changedPage(int index)
{
    qInfo() << "enter with index" << index;
    if ( index == page_anatomy )
        openedAnatomyPage();
    else if ( index == page_fMRI )
        openedfMRIPage();
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
    QString defaultDir = "/homes/deltabp/1/users/public/templates/humanMNI305/2mm";
    _lastTemplateDirectory = _savedQSettings.value("lastTemplateDirectory",defaultDir).toString();
    FUNC_INFO << "_lastTemplateDirectory" << _lastTemplateDirectory;
}

void MainWindow::writeQSettings()
{
    _savedQSettings.setValue("lastTemplateDirectory",_anatomyTemplateDirectory->currentText());
    _savedQSettings.sync();
}

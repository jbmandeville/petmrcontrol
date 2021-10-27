#include "petmrMain.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QCoreApplication::setOrganizationName("Martinos");
    QCoreApplication::setOrganizationDomain("http://www.nmr.mgh.harvard.edu");
    QCoreApplication::setApplicationVersion("3.0");
    QCoreApplication::setApplicationName("fastmap");
    QSettings fmSettings;
    QStringList defaultDirs = {"mouse","","rat","","NHP","","human",""};  // pairs of values = (id,path)
    QStringList templateDirectories = fmSettings.value("templateDirectories",defaultDirs).toStringList();
    FUNC_INFO << "template directories" << templateDirectories;

    QCoreApplication::setOrganizationName("Martinos");
    QCoreApplication::setOrganizationDomain("http://www.nmr.mgh.harvard.edu");
    QCoreApplication::setApplicationVersion("3.0");
    QCoreApplication::setApplicationName("petmrcontrol");

    MainWindow win;
    win.show();

    return app.exec();
}
